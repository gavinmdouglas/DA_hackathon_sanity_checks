### Check that correlations (and p-values) in Fig 2 are correct.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/")

tool2name <- read.table("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/mapfiles/tool_name_mapping.csv",
                           header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

fixed_hackathon_metadata_filt <- readRDS("fixed_hackathon_metadata_filt.rds")
fixed_hackathon_metadata_unfilt <- readRDS("fixed_hackathon_metadata_unfilt.rds")

unfilt_rhos_df <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/unfilt_rhos_df.rds")
filt_rhos_df <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/filt_rhos_df.rds")

unfilt_pvals_df <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/unfilt_pvals_df.rds")
filt_pvals_df <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/filt_pvals_df.rds")

unfilt_sig_percent_fix <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/unfilt_sig_percent_fix.rds")
filt_sig_percent_fix <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/filt_sig_percent_fix.rds")


### aldex2 % sig vs Sparsity - filt 
aldex2_vs_Sparsity_spearman_filt <- cor.test(filt_sig_percent_fix$aldex2, fixed_hackathon_metadata_filt[rownames(filt_sig_percent_fix), "Sparsity"], method = "spearman", exact = FALSE)
aldex2_vs_Sparsity_spearman_filt$estimate == filt_rhos_df[tool2name["aldex2", "clean"], "Sparsity"]
aldex2_vs_Sparsity_spearman_filt$p.value == filt_pvals_df[tool2name["aldex2", "clean"], "Sparsity"]

# corncob vs log(Read depth range) - filt
corncob_spearman_filt <- cor.test(filt_sig_percent_fix$corncob, fixed_hackathon_metadata_filt[rownames(filt_sig_percent_fix), "log(Read depth range)"], method = "spearman", exact = FALSE)
corncob_spearman_filt$estimate == filt_rhos_df[tool2name["corncob", "clean"], "log(Read depth range)"]
corncob_spearman_filt$p.value == filt_pvals_df[tool2name["corncob", "clean"], "log(Read depth range)"]


# maaslin2 vs richness - filt
maaslin2_spearman_filt <- cor.test(filt_sig_percent_fix$maaslin2, fixed_hackathon_metadata_filt[rownames(filt_sig_percent_fix), "Richness"], method = "spearman", exact = FALSE)
maaslin2_spearman_filt$estimate == filt_rhos_df[tool2name["maaslin2", "clean"], "Richness"]
maaslin2_spearman_filt$p.value == filt_pvals_df[tool2name["maaslin2", "clean"], "Richness"]


# lefse vs richness - filt
lefse_spearman_filt <- cor.test(filt_sig_percent_fix$lefse, fixed_hackathon_metadata_filt[rownames(filt_sig_percent_fix), "log(Sample size)"], method = "spearman", exact = FALSE)
lefse_spearman_filt$estimate == filt_rhos_df[tool2name["lefse", "clean"],"log(Sample size)"]
lefse_spearman_filt$p.value == filt_pvals_df[tool2name["lefse", "clean"], "log(Sample size)"]







### aldex2 % sig vs Sparsity - unfilt 
aldex2_vs_Sparsity_spearman_unfilt <- cor.test(unfilt_sig_percent_fix$aldex2, fixed_hackathon_metadata_unfilt[rownames(unfilt_sig_percent_fix), "Sparsity"], method = "spearman", exact = FALSE)
aldex2_vs_Sparsity_spearman_unfilt$estimate == unfilt_rhos_df[tool2name["aldex2", "clean"], "Sparsity"]
aldex2_vs_Sparsity_spearman_unfilt$p.value == unfilt_pvals_df[tool2name["aldex2", "clean"], "Sparsity"]

# corncob vs log(Read depth range) - unfilt
corncob_spearman_unfilt <- cor.test(unfilt_sig_percent_fix$corncob, fixed_hackathon_metadata_unfilt[rownames(unfilt_sig_percent_fix), "log(Read depth range)"], method = "spearman", exact = FALSE)
corncob_spearman_unfilt$estimate == unfilt_rhos_df[tool2name["corncob", "clean"], "log(Read depth range)"]
corncob_spearman_unfilt$p.value == unfilt_pvals_df[tool2name["corncob", "clean"], "log(Read depth range)"]


# maaslin2 vs richness - unfilt
maaslin2_spearman_unfilt <- cor.test(unfilt_sig_percent_fix$maaslin2, fixed_hackathon_metadata_unfilt[rownames(unfilt_sig_percent_fix), "Richness"], method = "spearman", exact = FALSE)
maaslin2_spearman_unfilt$estimate == unfilt_rhos_df[tool2name["maaslin2", "clean"], "Richness"]
maaslin2_spearman_unfilt$p.value == unfilt_pvals_df[tool2name["maaslin2", "clean"], "Richness"]


# lefse vs richness - unfilt
lefse_spearman_unfilt <- cor.test(unfilt_sig_percent_fix$lefse, fixed_hackathon_metadata_unfilt[rownames(unfilt_sig_percent_fix), "log(Sample size)"], method = "spearman", exact = FALSE)
lefse_spearman_unfilt$estimate == unfilt_rhos_df[tool2name["lefse", "clean"],"log(Sample size)"]
lefse_spearman_unfilt$p.value == unfilt_pvals_df[tool2name["lefse", "clean"], "log(Sample size)"]
