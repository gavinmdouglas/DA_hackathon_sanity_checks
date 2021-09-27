### Check that values related to Fig 1 and Fig 2 in the main text are correct.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/")

dataset2name <- read.table("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/mapfiles/dataset_name_mapping.csv",
                           header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

# Reported ranges of % sig
unfilt_sig_percent_fix <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/unfilt_sig_percent_fix.rds")
filt_sig_percent_fix <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/filt_sig_percent_fix.rds")

range(colMeans(filt_sig_percent_fix))
range(colMeans(unfilt_sig_percent_fix))


# Reported ranges of correlation with Aitchison's distance PERMANOVA effect size
unfilt_rhos_df <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/unfilt_rhos_df.rds")
filt_rhos_df <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/filt_rhos_df.rds")

range(unfilt_rhos_df[, "log(Aitch. dist. effect size)"])
range(filt_rhos_df[, "log(Aitch. dist. effect size)"])


# Specific means (and SDs to add into this section) - UNFILTERED
round(sapply(unfilt_sig_percent_fix, sd), 1)
round(sapply(unfilt_sig_percent_fix, mean), 1)

# Specific means (and SDs to add into this section) - FILTERED  
round(sapply(filt_sig_percent_fix, sd), 1)
round(sapply(filt_sig_percent_fix, mean), 1)

# Human HIV (3) mention
unfilt_sig_percent_fix["Human - HIV (3)", ]


# Built - Office and Freshwater - Arctic values
sort(unfilt_sig_percent_fix["Built - Office", ])
sort(unfilt_sig_percent_fix["Freshwater - Arctic", ])


# Wilcoxon (CLR) mention
length(which(unfilt_sig_percent_fix[, "wilcoxonclr"] > 90))


# Human - T1D (1) mention
sort(unfilt_sig_percent_fix["Human - T1D (1)", ])



