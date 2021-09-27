### Check that values of dataset metadata used in Fig 1 / 2 are correct.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/")

fixed_hackathon_metadata_filt <- readRDS("fixed_hackathon_metadata_filt.rds")
fixed_hackathon_metadata_unfilt <- readRDS("fixed_hackathon_metadata_unfilt.rds")

