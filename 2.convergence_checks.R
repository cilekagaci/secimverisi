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

plot_district_level('kasim2015_haziran2018meclis_districts*')
