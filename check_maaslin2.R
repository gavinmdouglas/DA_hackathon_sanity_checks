rm(list = ls(all.names = TRUE))

source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

library("Maaslin2")

datasets_to_test <- c("BISCUIT", "MALL", "Exercise", "ob_ross", "wood_plastic_kesy",
                      "Ji_WTP_DS", "cdi_vincent")

maaslin2_summary_out <- data.frame(matrix(NA, nrow = length(datasets_to_test), ncol = 2))
colnames(maaslin2_summary_out) <- c("filt", "unfilt")
rownames(maaslin2_summary_out) <- datasets_to_test

for(dataset in datasets_to_test) {

  dataset_filt_abun_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                      dataset,
                                      "/Fix_Results_0.1/fixed_non_rare_tables/",
                                      dataset,
                                      "_ASVs_table.tsv",
                                      sep = "")
  
  dataset_group_tab_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                      dataset,
                                      "/",
                                      dataset,
                                      "_metadata.tsv",
                                      sep = "")
  
  
  if(! file.exists(dataset_group_tab_filename)) {
    dataset_group_tab_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                        dataset,
                                        "/",
                                        dataset,
                                        "_metadata.csv",
                                        sep = "")
  }
                    
  dataset_filt_abun_maaslin2_sig <- run_maaslin2(dataset_filt_abun_filename, dataset_group_tab_filename, ncores = 40)
  
  dataset_filt_maaslin2_orig <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                dataset,
                                                "/Fix_Results_0.1/Maaslin2_out/significant_results.tsv",
                                                sep = ""),
                                          header=TRUE, sep="\t", check.names=FALSE)
  
  dataset_filt_maaslin2_orig_sig <- as.character(dataset_filt_maaslin2_orig$feature[which(dataset_filt_maaslin2_orig$qval < 0.05)])
  
  maaslin2_summary_out[dataset, "filt"] <- identical(sort(dataset_filt_abun_maaslin2_sig), sort(dataset_filt_maaslin2_orig_sig))
  
          
        dataset_unfilt_abun_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                              dataset,
                                              "/No_filt_Results/fixed_non_rare_tables/",
                                              dataset,
                                              "_ASVs_table.tsv",
                                              sep = "")
        
          
        dataset_unfilt_abun_maaslin2_sig <- run_maaslin2(dataset_unfilt_abun_filename, dataset_group_tab_filename, ncores = 40)
          
        dataset_unfilt_maaslin2_orig <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                        dataset,
                                                        "/No_filt_Results/Maaslin2_out/significant_results.tsv",
                                                        sep = ""),
                                                  header = TRUE, sep = "\t", check.names = FALSE)
          
        dataset_unfilt_maaslin2_orig_sig <- as.character(dataset_unfilt_maaslin2_orig$feature[which(dataset_unfilt_maaslin2_orig$qval < 0.05)])
  
        maaslin2_summary_out[dataset, "unfilt"] <- identical(sort(dataset_unfilt_abun_maaslin2_sig), sort(dataset_unfilt_maaslin2_orig_sig))
        
}
