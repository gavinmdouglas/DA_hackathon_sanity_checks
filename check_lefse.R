### Check LefSe results for all datasets
### Also evaluate the effect of removing ASVs missing in all features on the results

rm(list=ls(all.names=TRUE))

source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

outdir_path <- "/home/gavin/projects/hackathon/sanity_checks/lefse_testing"

datasets <- list.files("/home/jacob/projects/Hackathon/Studies/")

datasets <- datasets[-which(datasets %in% c("Simulation", "ob_zupancic"))]

invisible(mclapply(datasets, function(dataset) { format_for_lefse(dataset=dataset, outdir=outdir_path) }, mc.cores=40))

# Then on console navigated to ~/projects/hackathon/sanity_checks/lefse_testing
# Created lefse-formatted files with "lefse-format_input.py" and ran lefse with "run_lefse.py"
# (run over all input files with parallel)
# parallel -j 40 'lefse-format_input.py {} {.}.lefse -c 1 -u 2 -o 1000000' ::: *tsv
# parallel -j 40 'run_lefse.py {} {.}.out.tsv' ::: *lefse


tmp1 <- read.table("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/orig_inputs/lefse_format_file.tsv", header = TRUE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)
tmp2 <- read.table("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/BISCUIT_filt_lefse_input.tsv", header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)

tmp1 <- tmp1[-1, ]
colnames(tmp2) <- as.character(tmp2["id", ])
tmp2 <- tmp2[-c(1, 2), ]

tmp1_reorder <- tmp1[, colnames(tmp2)]

tmp1_reorder_num <- tmp1_reorder
tmp2_num <- tmp2

for(sample_name in colnames(tmp1_reorder)) {
  tmp1_reorder_num[, sample_name] <- as.numeric(tmp1_reorder[, sample_name])
  tmp2_num[, sample_name] <- as.numeric(tmp2[, sample_name])
}

tmp1_results <- read.table("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/orig_inputs/lefse_format_file_edit.lefse.results", header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)
#tmp1_results <- tmp1_results[-which(is.na(tmp1_results$V4)), ]
#tmp1_results$V5 <- as.numeric(tmp1_results$V5)

tmp2_results <- read.table("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/orig_inputs/tmp_biscuit_results.tsv", header = FALSE, sep ="\t", row.names = 1, stringsAsFactors = FALSE)
#tmp2_results <- tmp2_results[-which(is.na(tmp2_results$V4)), ]
#tmp2_results$V5 <- as.numeric(tmp2_results$V5)

tmp1_results <- tmp1_results[, -2]
tmp2_results <- tmp2_results[, -2]

tmp1_results$V5[which(tmp1_results$V5 == "-")] <- NA
tmp2_results$V5[which(tmp2_results$V5 == "-")] <- NA

tmp1_results$V5 <- as.numeric(tmp1_results$V5)
tmp2_results$V5 <- as.numeric(tmp2_results$V5)

tmp1_results <- tmp1_results[rownames(tmp2_results), ]
all.equal(tmp1_results, tmp2_results)

overlapping_rows <- rownames(tmp2_results)[which(rownames(tmp2_results) %in% rownames(tmp1_results))]
tmp1_results <- tmp1_results[overlapping_rows, ]
tmp2_results <- tmp2_results[overlapping_rows, ]


# After running lefse, I then read in the outputs and compared them with the original results.

lefse_out <- list()

for(dataset in datasets) {
  
  lefse_out[[dataset]] <- list()
  
  lefse_out[[dataset]][["orig_filt"]] <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                          dataset,
                                                          "/Fix_Results_0.1/Lefse_out/Lefse_results.tsv", sep=""),
                                                    sep="\t", row.names=1, stringsAsFactors = FALSE)
  
  lefse_out[[dataset]][["orig_unfilt"]] <- read.table(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                                            dataset,
                                                            "/No_filt_Results/Lefse_out/Lefse_results.tsv", sep=""),
                                                      sep="\t", row.names=1, stringsAsFactors = FALSE)
  
  
  # filt new
  filt_new_out <- paste("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/", dataset, "_filt_lefse_input.out.tsv", sep="")
  
  lefse_out[[dataset]][["new_filt"]] <- read.table(filt_new_out, sep="\t", row.names=1, stringsAsFactors = FALSE)
  
  filt_nomiss_new_out <- paste("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/", dataset, "_filt_nomiss_lefse_input.out.tsv", sep="")
  
  if(file.exists(filt_nomiss_new_out)) {
    lefse_out[[dataset]][["new_filt_nomiss"]] <- read.table(filt_nomiss_new_out, sep="\t", row.names=1, stringsAsFactors = FALSE)
  }
  
  # unfilt new
  unfilt_new_out <- paste("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/", dataset, "_unfilt_lefse_input.out.tsv", sep="")
  
  if(file.exists(unfilt_new_out)) {
    lefse_out[[dataset]][["new_unfilt"]] <- read.table(unfilt_new_out, sep="\t", row.names=1, stringsAsFactors = FALSE)
  }
  
  unfilt_nomiss_new_out <- paste("/home/gavin/projects/hackathon/sanity_checks/lefse_testing/", dataset, "_unfilt_nomiss_lefse_input.out.tsv", sep="")
  
  if(file.exists(unfilt_nomiss_new_out)) {
    lefse_out[[dataset]][["new_ununfilt_nomiss"]] <- read.table(unfilt_nomiss_new_out, sep="\t", row.names=1, stringsAsFactors = FALSE)
  }

}



# Get % of how many sig ASVs differ between old and new files.
lefse_out_perdiff <- lapply(lefse_out, compare_lefse_outputs)
lefse_out_perdiff_df <- do.call(rbind.data.frame, lefse_out_perdiff)

colnames(lefse_out_perdiff_df) <- c("new_filt", "new_filt_nomiss", "new_unfilt", "new_unfilt_nomiss")
rownames(lefse_out_perdiff_df) <- names(lefse_out)

write.table(x = lefse_out_perdiff_df,
            file="/home/gavin/projects/hackathon/sanity_checks/lefse_percent_diff_tab.tsv",
            col.names = NA, row.names = TRUE, quote=FALSE, sep="\t")
