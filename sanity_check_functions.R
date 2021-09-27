### This code was written independently of the code used for the manuscript.
### The idea was to verify a subset of results with independent code to catch possible errors


library(edgeR)
library(limma)
library(parallel)

parse_two_meta_groups <- function(meta_filename) {
  
  # If file doesn't exist then try adding ".csv" to end instead.
  if(! file.exists(meta_filename)) {
    meta_filename <- sub(".tsv", ".csv", meta_filename)
    
    if(! file.exists(meta_filename)) {
      meta_filename <- sub("_metadata.csv", "_meta.tsv", meta_filename)
    
      if(! file.exists(meta_filename)) { stop("file not found (with either .tsv or .csv extension or as _meta.tsv") }
      
    }
  }
    
  meta_tab <- read.table(meta_filename, header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE)
  comparison_groups <- table(meta_tab[, 1])
  
  if(length(comparison_groups) != 2) { stop("# of groups != 2") }
  
  return(list(group1=rownames(meta_tab)[which(meta_tab[, 1] == names(comparison_groups)[1])],
              group2=rownames(meta_tab)[which(meta_tab[, 1] == names(comparison_groups)[2])]))
}

limma_voom_two_group <- function(table, group1, group2, norm_approach) {
  
  group1_subset <- group1[which(group1 %in% colnames(table))]
  group2_subset <- group2[which(group2 %in% colnames(table))]
  
  # Create metadata df.
  metadata_tab <- data.frame(matrix(NA, nrow=(length(group1_subset) + length(group2_subset)), ncol=1))
  rownames(metadata_tab) <- c(group1_subset, group2_subset)
  colnames(metadata_tab) <- c("group")
  metadata_tab[group1_subset, "group"] <- "group1"
  metadata_tab[group2_subset, "group"] <- "group2"
  metadata_tab$group <- as.factor(metadata_tab$group)
  
  table <- table[, c(group1_subset, group2_subset)]
  
  counts <- floor(as.matrix(table))

  dge <- DGEList(counts = counts)
  
  ### Check if upper quartile method works for selecting reference
  upper_quartile_norm_test <- calcNormFactors(dge, method="upperquartile")
  
  summary_upper_quartile <- summary(upper_quartile_norm_test$samples$norm.factors)[3]
  if(is.na(summary_upper_quartile) | is.infinite(summary_upper_quartile)){
    message("Upper Quartile reference selection failed will use find sample with largest sqrt(read_depth) to use as reference")
    Ref_col <- which.max(colSums(sqrt(counts)))
    dge_norm <- calcNormFactors(dge, method = norm_approach, refColumn = Ref_col)
  }else{
    dge_norm <- calcNormFactors(dge, method=norm_approach)
  }
  
  model_matrix <- model.matrix(as.formula("~ group"), metadata_tab)
  
  voom_out <- voom(dge_norm, model_matrix, plot=FALSE)
  
  voom_out_fit <- lmFit(voom_out, model_matrix)
  voom_out_fit_eBayes <- eBayes(voom_out_fit)
  
  return(topTable(voom_out_fit_eBayes, coef = 2, number = nrow(dge_norm), sort.by="none"))
  
}


format_for_lefse <- function(dataset, outdir) {
  
  abun_rare_filt_path <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                               dataset, "/Fix_Results_0.1/fixed_rare_tables/", dataset, "_ASVs_table.tsv", sep="")
  abun_rare_filt_in <- read.table(abun_rare_filt_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  
  abun_rare_unfilt_path <- paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                 dataset, "/", dataset, "_ASVs_table_rare.tsv", sep="")
  
  # Check to see what format unfiltered table is in.
  abun_rare_unfilt_path_con <- file(abun_rare_unfilt_path)
  abun_rare_unfilt_line1 <- readLines(abun_rare_unfilt_path_con, n=1)
  close(abun_rare_unfilt_path_con)
  
  if(grepl("Constructed from biom file", abun_rare_unfilt_line1)) {
    abun_rare_unfilt_in <- read.table(abun_rare_unfilt_path, sep="\t", skip=1, header=T, row.names = 1,
                                      comment.char = "", quote="", check.names = FALSE)
  } else {
    abun_rare_unfilt_in <- read.table(abun_rare_unfilt_path, sep="\t", header=T, row.names = 1,
                                      comment.char = "", quote="", check.names = FALSE)
  }
  
  groups_in <- parse_two_meta_groups(paste("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/",
                                           dataset, "/", dataset, "_metadata.tsv", sep=""))
  
  group1_filt <- groups_in$group1[which(groups_in$group1 %in% colnames(abun_rare_filt_in))]
  group2_filt <- groups_in$group2[which(groups_in$group2 %in% colnames(abun_rare_filt_in))]
  
  group1_unfilt <- groups_in$group1[which(groups_in$group1 %in% colnames(abun_rare_unfilt_in))]
  group2_unfilt <- groups_in$group2[which(groups_in$group2 %in% colnames(abun_rare_unfilt_in))]
  
  abun_rare_filt_in <- abun_rare_filt_in[, c(group1_filt, group2_filt)]
  abun_rare_unfilt_in <- abun_rare_unfilt_in[, c(group1_unfilt, group2_unfilt)]
  
  #abun_rare_filt_in <- data.frame(sweep(abun_rare_filt_in, 2, colSums(abun_rare_filt_in), '/'), check.names = FALSE) * 100
  #abun_rare_unfilt_in <- data.frame(sweep(abun_rare_unfilt_in, 2, colSums(abun_rare_unfilt_in), '/'), check.names = FALSE) * 100
  
  filt_header <- data.frame(matrix(NA, nrow=2, ncol=length(c(group1_filt, group2_filt)) + 1))
  filt_header[1, ] <- c("comparison", rep("group1", length(group1_filt)), rep("group2", length(group2_filt)))
  filt_header[2, ] <- c("id", colnames(abun_rare_filt_in))
  filt_outfile <- paste(outdir, "/", dataset, "_filt_lefse_input.tsv", sep="")
  write.table(x = filt_header, file = filt_outfile, quote = FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
  write.table(x = abun_rare_filt_in, file = filt_outfile, quote = FALSE, col.names=FALSE, row.names=TRUE, sep="\t", append=TRUE)
  
  unfilt_header <- data.frame(matrix(NA, nrow=2, ncol=length(c(group1_unfilt, group2_unfilt)) + 1))
  unfilt_header[1, ] <- c("comparison", rep("group1", length(group1_unfilt)), rep("group2", length(group2_unfilt)))
  unfilt_header[2, ] <- c("id", colnames(abun_rare_unfilt_in))
  unfilt_outfile <- paste(outdir, "/", dataset, "_unfilt_lefse_input.tsv", sep="")
  write.table(x = unfilt_header, file = unfilt_outfile, quote = FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
  write.table(x = abun_rare_unfilt_in, file = unfilt_outfile, quote = FALSE, col.names=FALSE, row.names=TRUE, sep="\t", append=TRUE)
  
  
  if(min(rowSums(abun_rare_filt_in)) == 0) {
    abun_rare_filt_in_nomiss <- abun_rare_filt_in[-which(rowSums(abun_rare_filt_in) == 0), ]
    #abun_rare_filt_in_nomiss <- data.frame(sweep(abun_rare_filt_in_nomiss, 2, colSums(abun_rare_filt_in_nomiss), '/'), check.names = FALSE) * 100
    filt_nomiss_outfile <- paste(outdir, "/", dataset, "_filt_nomiss_lefse_input.tsv", sep="")
    write.table(x = filt_header, file = filt_nomiss_outfile, quote = FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
    write.table(x = abun_rare_filt_in_nomiss, file = filt_nomiss_outfile, quote = FALSE, col.names=FALSE, row.names=TRUE, sep="\t", append=TRUE)
  }
  
  if(min(rowSums(abun_rare_unfilt_in)) == 0) {
    abun_rare_unfilt_in_nomiss <- abun_rare_unfilt_in[-which(rowSums(abun_rare_unfilt_in) == 0), ]
    #abun_rare_unfilt_in_nomiss <- data.frame(sweep(abun_rare_unfilt_in_nomiss, 2, colSums(abun_rare_unfilt_in_nomiss), '/'), check.names = FALSE) * 100
    unfilt_nomiss_outfile <- paste(outdir, "/", dataset, "_unfilt_nomiss_lefse_input.tsv", sep="")
    write.table(x = unfilt_header, file = unfilt_nomiss_outfile, quote = FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
    write.table(x = abun_rare_unfilt_in_nomiss, file = unfilt_nomiss_outfile, quote = FALSE, col.names=FALSE, row.names=TRUE, sep="\t", append=TRUE)
  }
  
}


compare_lefse_outputs <- function(out_tabs) {
  
  orig_filt_sig <- gsub("f_", "", rownames(out_tabs$orig_filt)[which(out_tabs$orig_filt$V5 == "-")])
  orig_unfilt_sig <- gsub("f_", "", rownames(out_tabs$orig_unfilt)[which(out_tabs$orig_unfilt$V5 == "-")])
  
  orig_filt_sig_count <- length(orig_filt_sig)
  orig_unfilt_sig_count <- length(orig_unfilt_sig)
  
  percent_diff <- c()
  
  if(exists('new_filt', where=out_tabs)) {
    new_filt_sig <- gsub("f_", "", rownames(out_tabs$new_filt)[which(out_tabs$new_filt$V5 == "-")])
    percent_diff <- c(percent_diff, ((length(which(! new_filt_sig %in% orig_filt_sig)) + length(which(! orig_filt_sig %in% new_filt_sig))) / orig_filt_sig_count) * 100)
  } else {
    percent_diff <- c(percent_diff, NA)
  }
  
  if(exists('new_filt_nomiss', where=out_tabs)) {
    new_filt_nomiss_sig <- gsub("f_", "", rownames(out_tabs$new_filt_nomiss)[which(out_tabs$new_filt_nomiss$V5 == "-")])
    percent_diff <- c(percent_diff, ((length(which(! new_filt_nomiss_sig %in% orig_filt_sig)) + length(which(! orig_filt_sig %in% new_filt_nomiss_sig))) / orig_filt_sig_count) * 100)
  } else {
    percent_diff <- c(percent_diff, NA)
  }
  
  if(exists('new_unfilt', where=out_tabs)) {
    new_unfilt_sig <- gsub("f_", "", rownames(out_tabs$new_unfilt)[which(out_tabs$new_unfilt$V5 == "-")])
    percent_diff <- c(percent_diff, ((length(which(! new_unfilt_sig %in% orig_unfilt_sig)) + length(which(! orig_unfilt_sig %in% new_unfilt_sig))) / orig_unfilt_sig_count) * 100)
  } else {
    percent_diff <- c(percent_diff, NA)
  }
  
  if(exists('new_unfilt_nomiss', where=out_tabs)) {
    new_unfilt_nomiss_sig <- gsub("f_", "", rownames(out_tabs$new_unfilt_nomiss)[which(out_tabs$new_unfilt_nomiss$V5 == "-")])
    percent_diff <- c(percent_diff, ((length(which(! new_unfilt_nomiss_sig %in% orig_unfilt_sig)) + length(which(! orig_unfilt_sig %in% new_unfilt_nomiss_sig))) / orig_unfilt_sig_count) * 100)
  } else {
    percent_diff <- c(percent_diff, NA)
  }

  return(percent_diff)
  
}



run_corncob <- function(asv, meta) {
  ##read in the asv and meta tables
  ASV <- read.table(asv, header=TRUE, check.names = FALSE, row.names=1, sep="\t")
  meta <- read.table(meta, header=TRUE, row.names=1, check.names=FALSE, comment.char = "", quote = "", sep="\t")
  
  ##take the overlapping sample ids in the metadata and asv tables
  overlapping_ids <- colnames(ASV)[which(colnames(ASV) %in% rownames(meta))]
  ASV_ordered <- ASV[, overlapping_ids]
  meta_ordered <- meta[overlapping_ids,, drop=FALSE]
  
  ##use the ASV table and metadata table to create phyloseq otu table and sample data, then combine into a phyloseq object
  otu <- otu_table(ASV_ordered, taxa_are_rows = TRUE)
  sam_data <- sample_data(meta_ordered, errorIfNULL = TRUE)
  phylo <- merge_phyloseq(otu,sam_data)
  
  sample_var <- colnames(meta_ordered)[1]
  
  ##creates the formula for the differential test based on the sample_var provided in the input to the function
  my_formula <- as.formula(paste("~", sample_var, sep=" ", collapse=""))
  
  da_analysis <- differentialTest(formula = my_formula,
                                  phi.formula = my_formula,
                                  formula_null = ~ 1,
                                  phi.formula_null = my_formula,
                                  test= "Wald", boot = FALSE,
                                  data = phylo,
                                  fdr_cutoff = 0.05)
  
  
  return(da_analysis$significant_taxa)
}


run_deseq2 <- function(asv, meta){
  counts <- read.table(asv, sep="\t", row.names = 1, comment.char="", quote="", header=T, check.names=F)
  metadata <- read.table(meta, sep="\t", row.names = 1, header=T, comment.char="", check.names = F)
  col_data <- which( colnames(metadata)=="comparison" )
  colnames(metadata)[col_data] <- "grouping"
  
  metadata$Dup <- metadata$grouping
  
  ##take the overlapping sample ids in the metadata and asv tables
  overlapping_ids <- colnames(counts)[which(colnames(counts) %in% rownames(metadata))]
  ASV_ordered <- counts[, overlapping_ids]
  meta_ordered <- metadata[overlapping_ids,, drop=FALSE]
  #rownames(metadata) <- paste0("X",rownames(metadata))
  
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = meta_ordered,
                                design = ~ grouping)
  mid_dds <- DESeq(dds, sfType = "poscounts")
  
  mid_res <- results(mid_dds)
  
  resultsNames(mid_dds)
  mid_res = mid_res[order(mid_res$padj, na.last=NA), ]
  alpha = 0.05
  sighits <- mid_res[(mid_res$padj < alpha),]
  sighits
  sig_hits <- rownames(sighits)
  return(sig_hits)
  
}


### Taken from phyloseq authors at: https://joey711.github.io/phyloseq-extensions/edgeR.html
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}


run_edgeR <- function(ASV_table_path, groupings_path) {
  
  ASV_table <- read.table(ASV_table_path, sep="\t", row.names = 1, comment.char="", quote="", header=T, check.names=F)
  groupings <- read.table(groupings_path, sep="\t", row.names = 1, header=T, comment.char="", check.names = F)
  
  intersecting_ids <- colnames(ASV_table)[which(colnames(ASV_table) %in% rownames(groupings))]
  ASV_table <- ASV_table[, intersecting_ids]
  groupings <- groupings[intersecting_ids, , drop = FALSE]
  
  OTU <- phyloseq::otu_table(ASV_table, taxa_are_rows = T)
  sampledata <- phyloseq::sample_data(groupings, errorIfNULL = T)
  phylo <- phyloseq::merge_phyloseq(OTU, sampledata)
  
  test <- phyloseq_to_edgeR(physeq = phylo, group=colnames(groupings)[1])
  
  et = exactTest(test)
  
  tt = topTags(et, n=nrow(test$table), adjust.method="fdr", sort.by="PValue")
  outtable <- tt@.Data[[1]]
  
  return(rownames(outtable)[which(outtable$FDR < 0.05)])
}


run_metagenomeSeq <- function(ASV_table_path, grouping_path) {

  con <- file(ASV_table_path)
  file_1_line1 <- readLines(con,n=1)
  close(con)
  
  if(grepl("Constructed from biom file", file_1_line1)){
    ASV_table <- read.table(ASV_table_path, sep="\t", skip=1, header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }else{
    ASV_table <- read.table(ASV_table_path, sep="\t", header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }
  
  groupings <- read.table(grouping_path, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
  
  
  rows_to_keep <- intersect(colnames(ASV_table), rownames(groupings))
  groupings <- groupings[rows_to_keep,,drop=F]
  ASV_table <- ASV_table[,rows_to_keep]
  
  data_list <- list()
  data_list[["counts"]] <- ASV_table
  data_list[["taxa"]] <- rownames(ASV_table)
  
  pheno <- AnnotatedDataFrame(groupings)
  
  counts <- AnnotatedDataFrame(ASV_table)
  feature_data <- data.frame("ASV"=rownames(ASV_table),
                             "ASV2"=rownames(ASV_table))
  feature_data <- AnnotatedDataFrame(feature_data)
  rownames(feature_data) <- feature_data@data$ASV
  
  
  test_obj <- newMRexperiment(counts = data_list$counts, phenoData = pheno, featureData = feature_data)
  
  p <- cumNormStat(test_obj, pFlag = T)
  
  test_obj_norm <- cumNorm(test_obj, p=p)
  
  fromula <- as.formula(paste(~1, colnames(groupings)[1], sep=" + "))
  pd <- pData(test_obj_norm)
  mod <- model.matrix(fromula, data=pd)
  regres <- fitFeatureModel(test_obj_norm, mod)
  
  res_table <- MRfulltable(regres, number = length(rownames(ASV_table)))
  
  return(rownames(res_table)[which(res_table$adjPvalues < 0.05)])

}


run_maaslin2 <- function(ASV_table_path, grouping_path, tmp_dir = "/home/gavin/tmp", ncores=1) {  
  
  con <- file(ASV_table_path)
  file_1_line1 <- readLines(con,n=1)
  close(con)
  
  if(grepl("Constructed from biom file", file_1_line1)){
    ASV_table <- read.table(ASV_table_path, sep="\t", skip=1, header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }else{
    ASV_table <- read.table(ASV_table_path, sep="\t", header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }
  
  ASV_table <- data.frame(t(ASV_table), check.names = FALSE)
  
  groupings <- read.table(grouping_path, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
  
  intersecting_samples <- rownames(ASV_table)[which(rownames(ASV_table) %in% rownames(groupings))]
  
  ASV_table <- ASV_table[intersecting_samples, , drop = FALSE]
  groupings <- groupings[intersecting_samples, , drop = FALSE]
  
  orig_table_col_order <- colnames(ASV_table)
  ASV_table$ID <- rownames(ASV_table)
  ASV_table <- ASV_table[, c("ID", orig_table_col_order)]
  
  orig_metadata_col_order <- colnames(groupings)
  groupings$ID <- rownames(groupings)
  groupings <- groupings[, c("ID", orig_metadata_col_order)]
  
  tmp_input_tab_path <- paste(tmp_dir, "maaslin2_input_tab.tsv", sep = "/")
  tmp_metadata_tab_path <- paste(tmp_dir, "maaslin2_metadata_tab.tsv", sep = "/")
  tmp_output_path <- paste(tmp_dir, "maaslin2_output", sep = "/")
  
  write.table(x = ASV_table, file = tmp_input_tab_path, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(x = groupings, file = tmp_metadata_tab_path, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  maaslin2_output <- Maaslin2(tmp_input_tab_path,
                              tmp_metadata_tab_path,
                              tmp_output_path,
                              min_abundance = 0.0,
                              min_prevalence = 0.1,
                              normalization = "TSS",
                              transform = "AST",
                              analysis_method = "LM",
                              max_significance = 0.25,
                              fixed_effects = c('comparison'),
                              standardize = FALSE,
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE,
                              cores = ncores)
  
  unlink(c(tmp_input_tab_path, tmp_metadata_tab_path, tmp_output_path))
  
  return(gsub("^X", "", maaslin2_output$results$feature[which(maaslin2_output$results$qval < 0.05)]))
  
}


clr_transform_by_col <- function(in_df) {
  # Perform CLR transformation on table assuming samples are columns.
  return(data.frame(apply(in_df, 2, function(x){log(x) - mean(log(x))}), check.names=FALSE))
}

wilcoxon_2group_sig_features <- function(intable_path, metadata_path, convert_CLR=FALSE) {
  
  con <- file(intable_path)
  file_1_line1 <- readLines(con,n=1)
  close(con)
  
  if(grepl("Constructed from biom file", file_1_line1)){
    intable <- read.table(intable_path, sep="\t", skip=1, header=T, row.names = 1, 
                          comment.char = "", quote="", check.names = F)
  }else{
    intable <- read.table(intable_path, sep="\t", header=T, row.names = 1, 
                          comment.char = "", quote="", check.names = F)
  }
  
  group_in <- parse_two_meta_groups(metadata_path)
  
  group1_samples <- group_in$group1
  group2_samples <- group_in$group2
  
  group1_samples <- group1_samples[which(group1_samples %in% colnames(intable))]
  group2_samples <- group2_samples[which(group2_samples %in% colnames(intable))]
  
  group1_intable <- intable[, group1_samples]
  group2_intable <- intable[, group2_samples]
  if(convert_CLR) {
    group1_intable <- data.frame(t(clr_transform_by_col(group1_intable + 1)), check.names = FALSE)
    group2_intable <- data.frame(t(clr_transform_by_col(group2_intable + 1)), check.names = FALSE)
  } else {
    group1_intable <- data.frame(t(group1_intable), check.names = FALSE)
    group2_intable <- data.frame(t(group2_intable), check.names = FALSE)
  }
  
  wilcox_p <- c()
  for(feature in colnames(group1_intable)) {
    raw_wilcoxon_out <- wilcox.test(group1_intable[, feature], group2_intable[, feature], exact = FALSE)
    wilcox_p <- c(wilcox_p, raw_wilcoxon_out$p.value)
  }
  
  wilcox_fdr <- p.adjust(wilcox_p, "BH")
  
  return(colnames(group1_intable)[which(wilcox_fdr < 0.05)])
}


ttest_2group_sig_features <- function(intable_path, metadata_path) {
  
  con <- file(intable_path)
  file_1_line1 <- readLines(con,n=1)
  close(con)
  
  if(grepl("Constructed from biom file", file_1_line1)){
    intable <- read.table(intable_path, sep="\t", skip=1, header=T, row.names = 1, 
                          comment.char = "", quote="", check.names = F)
  }else{
    intable <- read.table(intable_path, sep="\t", header=T, row.names = 1, 
                          comment.char = "", quote="", check.names = F)
  }
  
  group_in <- parse_two_meta_groups(metadata_path)
  
  group1_samples <- group_in$group1
  group2_samples <- group_in$group2
  
  group1_samples <- group1_samples[which(group1_samples %in% colnames(intable))]
  group2_samples <- group2_samples[which(group2_samples %in% colnames(intable))]
  
  group1_intable <- intable[, group1_samples]
  group2_intable <- intable[, group2_samples]
  
  group1_intable <- data.frame(t(group1_intable), check.names = FALSE)
  group2_intable <- data.frame(t(group2_intable), check.names = FALSE)
  
  t.test_p <- c()
  for(feature in colnames(group1_intable)) {
    raw_ttest_out <- t.test(group1_intable[, feature], group2_intable[, feature])
    t.test_p <- c(t.test_p, raw_ttest_out$p.value)
  }
  
  t.test_fdr <- p.adjust(t.test_p, "BH")
  
  return(colnames(group1_intable)[which(t.test_fdr < 0.05)])
}
