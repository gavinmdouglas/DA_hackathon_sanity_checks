### Check ANCOM-II results on a few small datasets.

rm(list = ls(all.names = TRUE))

deps = c("exactRankTests", "nlme", "dplyr", "ggplot2", "compositions")
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(dep)
  }
  library(dep, character.only = TRUE)
}

source("/home/gavin/projects/hackathon/sanity_checks/ancom_2.1_code.R")

# BISCUIT unfiltered
BISCUIT_unfilt_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/BISCUIT/No_filt_Results/fixed_non_rare_tables/BISCUIT_ASVs_table.tsv",
                                  header=TRUE, sep="\t", row.names=1, check.names=FALSE)

BISCUIT_group_tab <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/BISCUIT/BISCUIT_metadata.tsv",
                                header=TRUE, sep="\t", row.names=1, check.names=FALSE)

BISCUIT_group_tab <- BISCUIT_group_tab[which(rownames(BISCUIT_group_tab) %in% colnames(BISCUIT_unfilt_abun)), , drop=FALSE]

BISCUIT_group_tab$Sample <- rownames(BISCUIT_group_tab)

BISCUIT_unfilt_abun_prepo <- feature_table_pre_process(feature_table = BISCUIT_unfilt_abun, meta_data = BISCUIT_group_tab, sample_var = 'Sample', 
                                                       group_var = NULL, out_cut = 0.05, zero_cut = 0.90,
                                                       lib_cut = 1000, neg_lb=FALSE)

#run ancom
p_adj_method = "BH"
alpha=0.05
adj_formula=NULL
rand_formula=NULL
BISCUIT_ANCOM_rerun <- ANCOM(feature_table = BISCUIT_unfilt_abun_prepo$feature_table, meta_data = BISCUIT_unfilt_abun_prepo$meta_data,
                             struc_zero = BISCUIT_unfilt_abun_prepo$structure_zeros, main_var = "comparison", p_adj_method = p_adj_method,
                             alpha=alpha, adj_formula = adj_formula, rand_formula = rand_formula)

BISCUIT_ANCOM_orig <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/BISCUIT/No_filt_Results/ANCOM_out/Ancom_res.tsv",
                                  header=TRUE, sep="\t", row.names=1, check.names=FALSE)

all.equal(BISCUIT_ANCOM_rerun$out, BISCUIT_ANCOM_orig)



# MALL filtered
MALL_filt_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/MALL/Fix_Results_0.1//fixed_non_rare_tables/MALL_ASVs_table.tsv",
                                  header=TRUE, sep="\t", row.names=1, check.names=FALSE)

MALL_group_tab <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/MALL/MALL_metadata.tsv",
                                header=TRUE, sep="\t", row.names=1, check.names=FALSE)

MALL_group_tab <- MALL_group_tab[which(rownames(MALL_group_tab) %in% colnames(MALL_filt_abun)), , drop=FALSE]

MALL_group_tab$Sample <- rownames(MALL_group_tab)

MALL_filt_abun_prepo <- feature_table_pre_process(feature_table = MALL_filt_abun, meta_data = MALL_group_tab, sample_var = 'Sample', 
                                                       group_var = NULL, out_cut = 0.05, zero_cut = 0.90,
                                                       lib_cut = 1000, neg_lb=FALSE)

#run ancom
p_adj_method = "BH"
alpha=0.05
adj_formula=NULL
rand_formula=NULL


MALL_ANCOM_rerun <- ANCOM(feature_table = MALL_filt_abun_prepo$feature_table, meta_data = MALL_filt_abun_prepo$meta_data,
                             struc_zero = MALL_filt_abun_prepo$structure_zeros, main_var = "comparison", p_adj_method = p_adj_method,
                             alpha=alpha, adj_formula = adj_formula, rand_formula = rand_formula)


MALL_ANCOM_orig <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/MALL/Fix_Results_0.1/ANCOM_out/Ancom_res.tsv",
                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)

all.equal(MALL_ANCOM_rerun$out, MALL_ANCOM_orig)

