rm(list = ls(all.names = TRUE))

setwd("/home/gavin/projects/hackathon/sanity_checks/fig3_feature_chat_intermediate_objects/")

filt_mean_RA <- readRDS("filt_study_tab_all_mean.rds")
unfilt_mean_RA <- readRDS("unfilt_study_tab_all_mean.rds")

filt_mean_RA_clean <- readRDS("mean_RA_sig_data_filt.rds")
unfilt_mean_RA_clean <- readRDS("mean_RA_sig_data_unfilt.rds")


adj_p_all_filt <- readRDS("Adj_p_all_filt.rds")
adj_p_all_unfilt <- readRDS("Adj_p_all_unfilt.rds")

## Unfiltered outliers:

### ALDEx2

aldex2_unfilt_sig_asvs <- rownames(adj_p_all_unfilt[which(adj_p_all_unfilt$aldex2 == 1), ])
median(unfilt_mean_RA[aldex2_unfilt_sig_asvs])

### ANCOM-II
ancom_unfilt_sig_asvs <- rownames(adj_p_all_unfilt[which(adj_p_all_unfilt$ancom == 1), ])
median(unfilt_mean_RA[ancom_unfilt_sig_asvs])


### DESeq2
deseq2_unfilt_sig_asvs <- rownames(adj_p_all_unfilt[which(adj_p_all_unfilt$deseq2 == 1), ])
median(unfilt_mean_RA[deseq2_unfilt_sig_asvs])

median(unfilt_mean_RA_clean[which(unfilt_mean_RA_clean$variable == "deseq2"), "value"], na.rm = TRUE)

## Filtered outliers:

### ALDEx2

## (NOTE THAT THE FILTERED ASVS IN THIS LIST DON'T APPEAR TO BE %'S BUT THAT WAS FIXED IN THE FINAL OBJECT USED FOR PLOTTING)
aldex2_filt_sig_asvs <- rownames(adj_p_all_filt[which(adj_p_all_filt$aldex2 == 1), ])
median(filt_mean_RA[aldex2_filt_sig_asvs])

median(filt_mean_RA_clean[which(filt_mean_RA_clean$variable == "aldex2"), "value"], na.rm = TRUE)

### ANCOM-II
ancom_filt_sig_asvs <- rownames(adj_p_all_filt[which(adj_p_all_filt$ancom == 1), ])
median(filt_mean_RA[ancom_filt_sig_asvs])

median(filt_mean_RA_clean[which(filt_mean_RA_clean$variable == "ancom"), "value"], na.rm = TRUE)
