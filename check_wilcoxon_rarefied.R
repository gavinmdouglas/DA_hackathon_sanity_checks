rm(list = ls(all.names = TRUE))

source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

datasets_to_test <- c("BISCUIT", "MALL", "Exercise", "ob_ross", "wood_plastic_kesy",
                      "Ji_WTP_DS", "cdi_vincent")

wilcoxon_summary_out <- data.frame(matrix(NA, nrow = length(datasets_to_test), ncol = 2))
colnames(wilcoxon_summary_out) <- c("filt", "unfilt")
rownames(wilcoxon_summary_out) <- datasets_to_test

for(dataset in datasets_to_test) {
  
  dataset_filt_abun_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                      dataset,
                                      "/Fix_Results_0.1/fixed_rare_tables/",
                                      dataset,
                                      "_ASVs_table.tsv",
                                      sep = "")
  
  dataset_group_tab_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                      dataset,
                                      "/",
                                      dataset,
                                      "_metadata.tsv",
                                      sep = "")
  
  dataset_filt_abun_wilcoxon_sig <- wilcoxon_2group_sig_features(dataset_filt_abun_filename, dataset_group_tab_filename, convert_CLR = FALSE)
  
  dataset_filt_wilcoxon_orig <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                 dataset,
                                                 "/Fix_Results_0.1/Wilcoxon_rare_out/Wil_rare_results.tsv",
                                                 sep = ""),
                                           header=TRUE, sep="\t", check.names=FALSE, row.names = 1)
  
  dataset_filt_wilcoxon_orig_sig <- rownames(dataset_filt_wilcoxon_orig)[which(p.adjust(dataset_filt_wilcoxon_orig$x, "BH") < 0.05)]
  
  wilcoxon_summary_out[dataset, "filt"] <- identical(sort(dataset_filt_abun_wilcoxon_sig), sort(dataset_filt_wilcoxon_orig_sig))
  
  
  dataset_unfilt_abun_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                        dataset,
                                        "/No_filt_Results/fixed_rare_tables/",
                                        dataset,
                                        "_ASVs_table.tsv",
                                        sep = "")
  
  
  dataset_unfilt_abun_wilcoxon_sig <- wilcoxon_2group_sig_features(dataset_unfilt_abun_filename, dataset_group_tab_filename, convert_CLR = FALSE)
  
  dataset_unfilt_wilcoxon_orig <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                   dataset,
                                                   "/No_filt_Results/Wilcoxon_rare_out/Wil_rare_results.tsv",
                                                   sep = ""),
                                             header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  
  dataset_unfilt_wilcoxon_orig_sig <- rownames(dataset_unfilt_wilcoxon_orig)[which(p.adjust(dataset_unfilt_wilcoxon_orig$x, "BH") < 0.05)]
  
  wilcoxon_summary_out[dataset, "unfilt"] <- identical(sort(dataset_unfilt_abun_wilcoxon_sig), sort(dataset_unfilt_wilcoxon_orig_sig))
  
}
