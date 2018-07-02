library('tidyverse')
library('stringr')
library('eiPack')
library('glue')
library('nnls')
select = dplyr::select

# Set the directory of this file as working directory
if (nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {
  script_folder = dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  frame_files <- lapply(sys.frames(), function(x) x$ofile)
  frame_files <- Filter(Negate(is.null), frame_files)
  script_folder <- dirname(frame_files[[length(frame_files)]])
}

setwd(glue('{script_folder}'))


source('utils.R')


data2018 = read_csv('processed_data/genel25_26_cb2018_meclis2018.csv')


mcmc_estimates_district_level <- function(pattern) {
  fname = (Sys.glob(pattern) %>% sort(decreasing = TRUE))[1]
  print(fname)
  load(fname)
  
  dfx = draws_list[(max(nrow(draws_list)-10, 0)):nrow(draws_list), ] %>% as_tibble %>%
    mutate(n = row_number()) %>%
    gather(key = metric, value = value, matches('lambda')) %>%
    separate(metric, into = c('dummy', 'row', 'col'), sep = '\\.') %>%
    mutate(value = ifelse(is.na(value), 0, value)) %>%
    group_by(row, col) %>%
    summarise(p50 = median(value), 
              mean = mean(value), 
              p05 = quantile(value, 0.05),
              p95 = quantile(value, 0.95),
              p005 = quantile(value, 0.005),
              p995 = quantile(value, 0.995)) %>%
    ungroup
  return(dfx)
}


mcmc_coefs_2015_2018meclis = mcmc_estimates_district_level('kasim2015_haziran2018meclis_districts*') %>%
  mutate(row = ifelse(row == 'row_others', 'OTH_Kasim', row),
         row = ifelse(row == 'row_slack', 'NEW_2018meclis', row),
         col = ifelse(col == 'col_others', 'OTH_2018meclis', col))


# Transition probability matrix from categories specified in the rows
# to categories specified in the columns. For example
# the value in  [ABS_Kasim, AKP_2018meclis] is the probability that a person
# who abstained in 2015 voted for AKP in 2018.
mcmc_coefs_2015_2018meclis %>%
  filter(col != 'col_slack') %>%
  transmute(row, 
            col, 
            value = mean) %>%
  spread(key = col, value = value)


