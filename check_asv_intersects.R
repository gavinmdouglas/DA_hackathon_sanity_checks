rm(list = ls(all.names = TRUE))

setwd("/home/gavin/projects/hackathon/sanity_checks/fig3_feature_chat_intermediate_objects/")

adj_p_all_filt <- readRDS("Adj_p_all_filt.rds")
adj_p_all_unfilt <- readRDS("Adj_p_all_unfilt.rds")

clean_intersect_unfilt <- readRDS("unfilt_bars_data_join.rds")
clean_intersect_filt <- readRDS("long_filt_melt_join.rds")

# Unfiltered - edgeR and LEfSe

clean_intersect_unfilt[which(clean_intersect_unfilt$Score == 1), ]

edger_unfilt_sig_asvs <- rownames(adj_p_all_unfilt[which(adj_p_all_unfilt$edger == 1), ])
edger_unfilt_sig_intersects <- rowSums(adj_p_all_unfilt[edger_unfilt_sig_asvs,], na.rm = TRUE)
length(which(edger_unfilt_sig_intersects == 1)) / length(edger_unfilt_sig_intersects)

clean_intersect_unfilt[which(clean_intersect_unfilt$variable == "edger" & clean_intersect_unfilt$Score == 1), "value"]



lefse_unfilt_sig_asvs <- rownames(adj_p_all_unfilt[which(adj_p_all_unfilt$lefse == 1), ])
lefse_unfilt_sig_intersects <- rowSums(adj_p_all_unfilt[lefse_unfilt_sig_asvs,], na.rm = TRUE)
length(which(lefse_unfilt_sig_intersects == 1)) / length(lefse_unfilt_sig_intersects)

clean_intersect_unfilt[which(clean_intersect_unfilt$variable == "lefse" & clean_intersect_unfilt$Score == 1), "value"]



# Percentage of ASVs called by more than 12 tools between filtered and unfiltered

clean_intersect_unfilt$Score_num <- as.numeric(as.character(clean_intersect_unfilt$Score))

clean_intersect_unfilt_above_12 <- clean_intersect_unfilt[which(clean_intersect_unfilt$Score_num > 12), ]

clean_intersect_unfilt_above_12_summed <- aggregate(value ~ variable, data = clean_intersect_unfilt_above_12, FUN = sum)

mean(clean_intersect_unfilt_above_12_summed$value)
sd(clean_intersect_unfilt_above_12_summed$value)


clean_intersect_filt$Score_num <- as.numeric(as.character(clean_intersect_filt$Score))

clean_intersect_filt_above_12 <- clean_intersect_filt[which(clean_intersect_filt$Score_num > 12), ]

clean_intersect_filt_above_12_summed <- aggregate(value ~ variable, data = clean_intersect_filt_above_12, FUN = sum)

mean(clean_intersect_filt_above_12_summed$value)
sd(clean_intersect_filt_above_12_summed$value)



