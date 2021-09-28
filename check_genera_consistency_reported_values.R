### Check values reported in text that are based simply formatted genus tables (because the raw files include many genera that are non-comparable across studies)

rm(list = ls(all.names = TRUE))

diarrhea_combined_overlap <- readRDS("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/consistency_analysis_RDS_out/diarrhea_combined_overlap.rds")

diarrhea_combined_overlap_binary <- readRDS("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/consistency_analysis_RDS_out/diarrhea_outputs_binary_clean_combined.rds")

# ALDEx2
mean(colSums(diarrhea_combined_overlap_binary[, grep("aldex2", colnames(diarrhea_combined_overlap_binary))]))
sd(colSums(diarrhea_combined_overlap_binary[, grep("aldex2", colnames(diarrhea_combined_overlap_binary))]))


# edger
mean(colSums(diarrhea_combined_overlap_binary[, grep("edger", colnames(diarrhea_combined_overlap_binary))]))
sd(colSums(diarrhea_combined_overlap_binary[, grep("edger", colnames(diarrhea_combined_overlap_binary))]))

