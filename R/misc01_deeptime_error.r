library(pacman)

p_load(tidyverse, tidybayes, reprex)
p_load_gh('willgearty/deeptime')

tb <- 
  tibble(time = 0:100) %>%
  group_by_all() %>%
  do(tibble(value = rnorm(100, .$time, 10),
            group = sample(1:4, 100, TRUE)))

gg <- 
  tb %>%
  ggplot() +
  stat_lineribbon(aes(x = time, y = value)) +
  scale_x_reverse() +
  coord_cartesian(xlim = c(0, 100), 
                  ylim = range(tb$value),
                  expand = FALSE) +
  facet_wrap(~ group) +
  scale_fill_brewer()

gggeo_scale(gg)
