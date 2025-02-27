### Re-run corncob on a few datasets

rm(list = ls(all.names = TRUE))

source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

library("corncob")
library("phyloseq")
library("parallel")

datasets_to_test <- c("BISCUIT", "MALL", "Exercise", "ob_ross", "wood_plastic_kesy",
                      "Ji_WTP_DS", "cdi_vincent")

corncob_summary_out <- data.frame(matrix(NA, nrow = length(datasets_to_test), ncol = 2))
colnames(corncob_summary_out) <- c("filt", "unfilt")
rownames(corncob_summary_out) <- datasets_to_test

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
  
  dataset_filt_abun_corncob_sig <- run_corncob(asv = dataset_filt_abun_filename, meta = dataset_group_tab_filename)
  
  dataset_filt_corncob_orig <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                dataset,
                                                "/Fix_Results_0.1/Corncob_out/Corncob_results.tsv",
                                                sep = ""),
                                          header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  
  dataset_filt_corncob_orig_sig <- rownames(dataset_filt_corncob_orig)[which(dataset_filt_corncob_orig$x < 0.05)]
  
  corncob_summary_out[dataset, "filt"] <- identical(sort(dataset_filt_abun_corncob_sig), sort(dataset_filt_corncob_orig_sig))
  
          
        dataset_unfilt_abun_filename <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                              dataset,
                                              "/No_filt_Results/fixed_non_rare_tables/",
                                              dataset,
                                              "_ASVs_table.tsv",
                                              sep = "")
        
          
        dataset_unfilt_abun_corncob_sig <- run_corncob(asv = dataset_unfilt_abun_filename, meta = dataset_group_tab_filename)
          
        dataset_unfilt_corncob_orig <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                        dataset,
                                                        "/No_filt_Results/Corncob_out/Corncob_results.tsv",
                                                        sep = ""),
                                                  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
          
        dataset_unfilt_corncob_orig_sig <- rownames(dataset_unfilt_corncob_orig)[which(dataset_unfilt_corncob_orig$x < 0.05)]
          
        corncob_summary_out[dataset, "unfilt"] <- identical(sort(dataset_unfilt_abun_corncob_sig), sort(dataset_unfilt_corncob_orig_sig))
        
}
