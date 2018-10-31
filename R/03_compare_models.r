# data manipulation
library(tidyverse)
library(janitor)
library(ggridges)
#devtools::install_github("mjskay/tidybayes")
library(tidybayes)

# parallel processing
library(parallel)
library(future)

# bayes
library(arm)
library(rstanarm)
library(bayesplot)
source('../R/helper03_stan_utility.r')

# misc
library(splines)
library(pROC)
library(ROCR)
source('../R/helper01_process_foo.r')
source('../R/helper02_plot_foo.r')

# important constants
options(mc.cores = parallel::detectCores())
future::plan(strategy = multicore)

# get data in
longi <- read_rds('../data/longitude.rds')
survi <- read_rds('../data/survival.rds')
counti <- read_rds('../data/counting.rds')

# form of data that was analyzed
counti_trans <- prepare_analysis(counti)

# read in model fits
model_key <- c('Past and vary', 
               'Past but no vary',
               'No past but vary', 
               'No past or vary')
disc_fit <- read_rds('../data/disc_fit.rds')

## diagnostics
#np <- map(disc_fit, nuts_params)
#ll <- map(disc_fit, log_posterior)
## mcmc
#rhat <- map(disc_fit, ~ mcmc_rhat(rhat(.x)))
#neff <- map(disc_fit, ~ mcmc_neff(neff_ratio(.x)))
## hmc/nuts
#accept <- map2(np, ll, ~ mcmc_nuts_acceptance(.x, .y))
#diverg <- map2(np, ll, ~ mcmc_nuts_divergence(.x, .y))
#steps <- map2(np, ll, ~ mcmc_nuts_stepsize(.x, .y))
#treed <- map2(np, ll, ~ mcmc_nuts_treedepth(.x, .y))
#energy <- map2(np, ll, ~ mcmc_nuts_energy(.x, .y))

# bayes R2
br2 <- map(disc_fit, bayes_R2)
br2_gg <- reshape2::melt(br2) %>%
  rename(model = L1) %>%
  mutate(model_name = plyr::mapvalues(model, unique(model), model_key),
         model_name = factor(model_name, levels = model_key)) %>%
  ggplot(aes(x = value, y = model_name)) +
  geom_halfeyeh(width = c(0.5, 0.8)) +
  labs(x = 'Bayesian R^2', y = 'Model')
ggsave(filename = '../results/figure/bayes_r2.png', plot = br2_gg,
       width = 6, height = 6)


## estimate of out-of-sample performance
## these should be parallelized; done using future
#cl <- map(disc_fit, ~ future::future(loo(.x)))
#cl <- map(cl, ~ future::value(.x))
#compare_loo_tab <- loo::compare(x = cl)
#rownames(compare_loo_tab) <- plyr::mapvalues(rownames(compare_loo_tab),
#                                             sort(rownames(compare_loo_tab)),
#                                             model_key)
#xtable::print.xtable(xtable::xtable(compare_loo_tab), 
#                     file = '../results/loo_tab_draft.tex')
#
#wl <- map(disc_fit, ~ future::future(waic(.x)))
#wl <- map(wl, ~ future::value(.x))
#compare_waic_tab <- loo::compare(x = wl)
#rownames(compare_waic_tab) <- plyr::mapvalues(rownames(compare_waic_tab),
#                                              sort(rownames(compare_waic_tab)),
#                                              model_key)
#xtable::print.xtable(xtable::xtable(compare_waic_tab), 
#                     file = '../results/waic_tab_draft.tex')
## i should find a way to output these nicely


# posterior probability of observation surviving
pp <- map(disc_fit, ~ future::future(posterior_linpred(.x, 
                                                       transform = TRUE, 
                                                       draws = 100)))
pp_prob <- map(pp, ~ future::value(.x))

# use linpred to get ROC curve
# probability gives finer resolution
#fut_eroc <- map(pp_prob, ~ future::future(apply(.x, 1, function(y) roc(counti_trans$event, y))))
#eroc <- map(fut_eroc, ~ future::value(.x))
eroc <- map(pp_prob, ~ apply(.x, 1, function(y) roc(counti_trans$event, y)))

# extract the parts of the ROC plot
roc_df <- map(eroc, function(x) 
           imap(x, ~ tibble(sim = .y, 
                            fpr = 1 - .x$specificities,
                            tpr = .x$sensitivities))) %>%
  map(., ~ reduce(.x, rbind)) %>%
  imap(., ~ add_column(.x, mod = .y)) %>%
  reduce(., rbind) %>%
  as.tibble(.) %>%
  # long hand way that doesn't need plyr
  mutate(model = case_when(mod == 1 ~ model_key[1],
                           mod == 2 ~ model_key[2],
                           mod == 3 ~ model_key[3],
                           mod == 4 ~ model_key[4])) 

#roc_df_back <- roc_df                  # extra for background


# plot the curves from the probabilities
cur <- roc_df %>%
  ggplot(aes(x = fpr, 
             y = tpr, 
             group = sim)) +
#  geom_line(data = roc_df_back,
#            colour = 'grey10') +
  geom_line(data = roc_df,
            mapping = aes(colour = model),
            alpha = 0.01) +
  facet_wrap(model ~ ., nrow = 2) +
  coord_equal(ratio = 1) +
  labs(x = 'False Positive Rate',
       y = 'True Positive Rate') +
  theme(legend.position = 'none')
ggsave(filename = '../results/figure/roc_curve.png', plot = cur,
       width = 5, height = 8)


# get AUC values for the above
# plot as histogram
auc_hist <- map(eroc, ~ map(.x, function(y) auc(y)[[1]])) %>%
  reshape2::melt(.) %>%
  as.tibble %>%
  rename(model = L1,
         draw = L2) %>%
  mutate(model_name = plyr::mapvalues(model, unique(model), model_key),
         model_name = factor(model_name, levels = model_key)) %>%
  ggplot(aes(x = value, y = model_name)) +
  geom_halfeyeh(.width = c(0.5, 0.8)) +
  labs(x = 'AUC ROC', y = NULL)
ggsave(filename = '../results/figure/auc_hist.png', plot = auc_hist,
       width = 6, height = 6)

# roc as timeseries to see best and worst times
roc_ts <- plot_roc_series(counti_trans, pp_prob, model_key) +
  scale_x_reverse()
ggsave(filename = '../results/figure/roc_ts.png', plot = roc_ts,
