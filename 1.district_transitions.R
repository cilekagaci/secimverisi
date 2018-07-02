library('tidyverse')
library('stringr')
library('eiPack')
library('glue')
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

# Read election data and keep vote count columns for 2015 November elections and
# 2018 June elections.
data = read_csv('genel25_26_referandum2017_27meclis_cb2018.csv') %>%
  mutate(province.NUTS1 = province.NUTS3 %>% substr(1, 3)) %>%
  select(-matches('Haziran|2018cb|referandum|EVET|HAYIR')) %>%
  # ABS is abstinance.
  mutate(ABS.Kasim = voter.count.Kasim - (total.valid.vote.count.Kasim + invalid.vote.count.Kasim),
         ABS.2018meclis = voter.count.2018meclis - voted.count.2018meclis) 

row_names = c('AKP.Kasim', 
              'CHP.Kasim', 
              'HDP.Kasim', 
              'MHP.Kasim', 
              'ABS.Kasim', 
              'invalid.vote.count.Kasim')
col_names = c('AKP.2018meclis', 
              'CHP.2018meclis', 
              'MHP.2018meclis', 
              'HDP.2018meclis', 
              'IYI.2018meclis', 
              'ABS.2018meclis', 
              'invalid.vote.count.2018meclis')

data %>% 
  run_simulation('kasim2015_haziran2018meclis_districts', 
                 row_names, col_names,
                 major_steps = 2000,
                 burnin = 0,
                 sample = 10,
                 thin = 1000,
                 save_every = 10, 
                 tune = TRUE,
                 starting_step = NULL,
                 row_total_name = 'voter.count.Kasim', 
                 col_total_name = 'voter.count.2018meclis')
