#' Get optimal cut from ROC results
#' 
#' Uses ROCR. I got this from somewhere on the internet.
#'
#' @param perf performance object
#' @param pred prediction object
#' @return vector with sensitivity, specificity, and cutpoint
opt_cut <- function(perf, pred) {
  cut_ind <- mapply(FUN = function(x, y, p) {
                      d <- (x - 0)^2 + (y - 1)^2
                      ind <- which(d == min(d))
                      c(sensitivity = y[[ind]], 
                        specificty = 1 - x[[ind]],
                        cutoff = p[[ind]])
               }, perf@x.values, perf@y.values, pred@cutoffs)
  cut_ind
}


#' Get optimal cutpoints for a lot of ROC results
#'
#' Looped wrapper around opt_cut to get optimal cutpoint for a given posterior predictive distribution
#'
#' @param pp_prob list of posterior predictive draws, where each list element is the PPD for one model.
#' @param target the actual values that are being estimated
get_cutpoints <- function(pp_prob, target) {
  pred <- plyr::alply(pp_prob, 1, function(x) prediction(x, target))
  perf <- map(pred, ~ performance(.x, measure = 'tpr', x.measure = 'fpr'))
  cutpoints <- map2(perf, pred, opt_cut)
}


#' Reassign class membership based on new cutpoints
#'
#' Give probabilities and new cutpoint level, reassigns classes.
#'
#' @param pp_prob list of posterior predictive draws, where each list element is the PPD for one model.
#' @param list_cutpoint list of new cutpoints as estimated from get_cutpoints
#' @return list of reassigned class predictions
cut_newpoint <- function(pp_prob, list_cutpoint) {
  pp_est_new <- list()
  for(jj in seq(length(pp_prob))) {
    mm <- matrix(ncol = ncol(pp_prob[[jj]]), nrow = nrow(pp_prob[[jj]]))
    for(ii in seq(nrow(mm))) {
      mm[ii, ] <- pp_prob[[jj]][ii, ] > list_cutpoint[[jj]]
    }
    pp_est_new[[jj]] <- mm * 1
  }
  pp_est_new
}

#' Posterior roc
#' 
#' This is mostly a convenient function to quickly get ROC estimates for entire posterior distribution (awkward matrix form).
#' 
#' @param x matrix of posterior 
#' @param y matrix with true data. named element event is vector of 0/1
post_roc <- function(x, y) apply(x, 1, function(a) roc(y$event, a))

#' Extract posterior AUC, safely
#'
#' Wrapped with safely to prevent response has only 1 level stuff
#'
#' @param pred matrix of predicted 0/1 values
#' @param counti_fold testing dataset
#' @return list w/ two elements result and error
get_auc_time <- function(pred, counti_fold) {
  fold_time <- split(counti_fold, counti_fold$mybin)

  pred_time <- split(t(pred), counti_fold$mybin) %>%
    map2(., fold_time, ~ matrix(.x, nrow = nrow(.y))) %>%
    map(., t)

  safe_pr <- safely(post_roc)
  auc_time <- map2(pred_time, fold_time, safe_pr) %>%
    map(., function(x) map_dbl(x$result, ~ auc(.x)))
  auc_time
}






# functions below adapted from https://gist.github.com/traversc/1446ebe1dcc2d84bccdca781d4f1fa2a

#' Fast calculation of ROC values
#'
#' calculate roc, giving both true positive rate and false positive rate
#' 
fast_roc <- function(probs, class) {
  class_sorted <- class[order(probs, decreasing=T)]
  TPR <- cumsum(class_sorted) / sum(class)
  FPR <- cumsum(class_sorted == 0) / sum(class == 0)
  return(list(tpr=TPR, fpr=FPR))
}

# Helpful function adapted from: https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
fast_auc <- function(probs, class) {
  x <- probs
  y <- class
  x1 = x[y==1]; n1 = length(x1); 
  x2 = x[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
  return(auc)
}
