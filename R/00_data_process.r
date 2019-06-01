library(pacman)

p_load(tidyverse, janitor, here, furrr,
       sp, dismo, raster, rgdal, XML, maptools, geosphere, GeoRange,
       scales)

# misc
source(here('R', 'helper01_process_foo.r'))
source(here('R', 'helper03_misc_foo.r'))

# parallel set up
plan(multiprocess)


# process timescale information
gts <- read_csv("https://raw.githubusercontent.com/japhir/stratPlot/master/GTS_colours.csv") %>%
  mutate_if(is.character, str_to_lower)
write_csv(gts, path = here('data', 'geologic_time_scale.csv'))

# process the raw data into an anlyzable format

bin_width <- 1
age_max <- 63
form <- 'counti'

# occurrence information
neptune <- list.files(path = here('data'), pattern = 'nannotax', full.names = TRUE) %>%
  purrr::map(., read_tsv) %>%
  purrr::map(., clean_names) %>%
  purrr::reduce(., rbind)

# map
sp <- CRS("+proj=cea +lat_0=0 +lon_0=0 +lat_ts=30 +a=6371228.0 +units=m")
globe.map <- readOGR(here('data', 'ne_10m_coastline.shp')) # from natural earth
proj4string(globe.map) <- sp

# climate data
mgca <- read_tsv(here('data', 'cramer', 'cramer_temp.txt')) %>% 
  clean_names

# longhust codes for atlantic
lnghst_code <- c('FKLD', 'BRAZ', 'BENG', 'GUIN', 'CNRY', 'GUIA', 'NECS', 
                 'ARCT', 'SARC', 'SATL', 'ETRA', 'WTRA', 'CARB', 'NATR', 
                 'NASE', 'NASW', 'MEDI', 'GFST', 'NADR', 'NWCS')


# keep something constant across datasets
partial_raw_to_clean <- partial(raw_to_clean,
                                nano = neptune,
                                form = 'count',
                                bin_width = bin_width, 
                                age_max = age_max,
                                sp = sp,
                                mgca = mgca)
counti <- partial_raw_to_clean()
counti_restrict_time <- partial_raw_to_clean(restrict = TRUE)
counti_restrict_local <- partial_raw_to_clean(prov = lnghst_code)


# try with multiple bin_width -s
partial_short_bin <- partial(raw_to_clean,
                             nano = neptune,
                             form = 'count',
                             age_max = age_max,
                             sp = sp,
                             mgca = mgca)


prep_short <- function(width = bin_width) {

  varname <- paste0("bin_", width)     # tidy eval stuff

  partial_short_bin(bin_width = width)  %>%
  dplyr::select(fullname, mybin, latext, nocc, ncell, maxgcd) %>%
  rename(!!varname := mybin) %>%       # tidy eval stuff
  gather(key = 'key', value = 'value', 
         -fullname, -maxgcd, -latext, -nocc, -ncell) %>%
  separate(key, c('bin', 'width'), sep = '_') %>%
  mutate(width = parse_number(width)) %>%
  dplyr::select(maxgcd, latext, nocc, ncell, width)

}

widths <- c(0.1, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 5)
vary_bin <- future_map_dfr(widths, ~ prep_short(.x))
#counti_01 <- prep_short(0.1) 
#
#counti_025 <- prep_short(0.25) 
#
#counti_05 <- prep_short(0.5)
#
#counti_1 <- prep_short(1)
#
#counti_15 <- prep_short(1.5)
#
#counti_2 <- prep_short(2)
#
#counti_25 <- prep_short(2.5)
#
#counti_3 <- prep_short(3)
#
#counti_5 <- prep_short(5)
#
#vary_bin <- bind_rows(counti_01, 
#                      counti_025, 
#                      counti_05, 
#                      counti_1, 
#                      counti_15, 
#                      counti_2, 
#                      counti_25, 
#                      counti_3, 
#                      counti_5)


# write to file
write_rds(counti, path = here('data', 'counting.rds'))
write_rds(counti_restrict_time, path = here('data', 'counting_restrict_time.rds'))
write_rds(counti_restrict_local, path = here('data', 'counting_restrict_local.rds'))
write_rds(vary_bin, path = here('data', 'counting_vary_binwidth.rds'))
