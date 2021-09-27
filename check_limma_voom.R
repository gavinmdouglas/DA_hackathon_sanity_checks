### Check a subset of limma-voom results by hand ###

# Start with two major outliers
# HIV_nog and limma_voom_TMMwsp (filtered)
# crc_zeller and limma_voom_TMMwsp (filtered)

# Then also look at three other cases:
# t1d_mej and limm_voom (unfiltered)
# MALL and limma_voom (filtered)
# glass and limma_voom (unfiltered)


################ HIV_nog and limma_voom_TMMwsp (filtered) ##################
rm(list=ls(all.names=TRUE))

options(warn=2)

source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

HIV_nog_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian/Fix_Results_0.1/fixed_non_rare_tables/hiv_noguerajulian_ASVs_table.tsv",
                         header=TRUE, sep="\t", row.names=1, check.names=FALSE)

HIV_nog_groups <- parse_two_meta_groups("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian/hiv_noguerajulian_metadata.tsv")

HIV_nog_limma_voom_TMM_new_out <- limma_voom_two_group(table=HIV_nog_abun,
                                                   group1=HIV_nog_groups$group1,
                                                   group2=HIV_nog_groups$group2,
                                                   norm_approach = "TMM")

HIV_nog_limma_voom_TMMwsp_new_out <- limma_voom_two_group(table=HIV_nog_abun,
                                                          group1=HIV_nog_groups$group1,
                                                          group2=HIV_nog_groups$group2,
                                                          norm_approach = "TMMwsp")


HIV_nog_limma_voom_TMMwsp_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian/Fix_Results_0.1/Limma_voom_TMMwsp/limma_voom_tmmwsp_res.tsv",
                                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)

HIV_nog_limma_voom_TMM_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian/Fix_Results_0.1/limma_voom_tmm_out/limma_voom_tmm_res.tsv",
                                              header=TRUE, sep="\t", row.names=1, check.names=FALSE)

HIV_nog_limma_voom_TMMwsp_orig_out_sig <- rownames(HIV_nog_limma_voom_TMMwsp_orig_out)[which(HIV_nog_limma_voom_TMMwsp_orig_out$adj.P.Val < 0.05)]
HIV_nog_limma_voom_TMM_orig_out_sig <- rownames(HIV_nog_limma_voom_TMM_orig_out)[which(HIV_nog_limma_voom_TMM_orig_out$adj.P.Val < 0.05)]


HIV_nog_limma_voom_TMMwsp_new_out_sig <- rownames(HIV_nog_limma_voom_TMMwsp_new_out)[which(HIV_nog_limma_voom_TMMwsp_new_out$adj.P.Val < 0.05)]
HIV_nog_limma_voom_TMM_new_out_sig <- rownames(HIV_nog_limma_voom_TMM_new_out)[which(HIV_nog_limma_voom_TMM_new_out$adj.P.Val < 0.05)]

# Key check to see if results are the same:
identical(HIV_nog_limma_voom_TMMwsp_new_out_sig, HIV_nog_limma_voom_TMMwsp_orig_out_sig)
identical(HIV_nog_limma_voom_TMM_new_out_sig, HIV_nog_limma_voom_TMM_orig_out_sig)

################ crc_zeller and limma_voom_TMMwsp (filtered) ##################

rm(list=ls(all.names=TRUE))
source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

crc_zeller_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/crc_zeller/Fix_Results_0.1/fixed_non_rare_tables/crc_zeller_ASVs_table.tsv",
                           header=TRUE, sep="\t", row.names=1, check.names=FALSE)

crc_zeller_groups <- parse_two_meta_groups("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/crc_zeller/crc_zeller_metadata.tsv")

crc_zeller_limma_voom_TMM_new_out <- limma_voom_two_group(table=crc_zeller_abun,
                                                       group1=crc_zeller_groups$group1,
                                                       group2=crc_zeller_groups$group2,
                                                       norm_approach = "TMM")

crc_zeller_limma_voom_TMMwsp_new_out <- limma_voom_two_group(table=crc_zeller_abun,
                                                          group1=crc_zeller_groups$group1,
                                                          group2=crc_zeller_groups$group2,
                                                          norm_approach = "TMMwsp")


crc_zeller_limma_voom_TMMwsp_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/crc_zeller/Fix_Results_0.1/Limma_voom_TMMwsp/limma_voom_tmmwsp_res.tsv",
                                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)

crc_zeller_limma_voom_TMM_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/crc_zeller/Fix_Results_0.1/limma_voom_tmm_out/limma_voom_tmm_res.tsv",
                                              header=TRUE, sep="\t", row.names=1, check.names=FALSE)

crc_zeller_limma_voom_TMMwsp_orig_out_sig <- rownames(crc_zeller_limma_voom_TMMwsp_orig_out)[which(crc_zeller_limma_voom_TMMwsp_orig_out$adj.P.Val < 0.05)]
crc_zeller_limma_voom_TMM_orig_out_sig <- rownames(crc_zeller_limma_voom_TMM_orig_out)[which(crc_zeller_limma_voom_TMM_orig_out$adj.P.Val < 0.05)]


crc_zeller_limma_voom_TMMwsp_new_out_sig <- rownames(crc_zeller_limma_voom_TMMwsp_new_out)[which(crc_zeller_limma_voom_TMMwsp_new_out$adj.P.Val < 0.05)]
crc_zeller_limma_voom_TMM_new_out_sig <- rownames(crc_zeller_limma_voom_TMM_new_out)[which(crc_zeller_limma_voom_TMM_new_out$adj.P.Val < 0.05)]

# Key check to see if results are the same:
identical(crc_zeller_limma_voom_TMMwsp_new_out_sig, crc_zeller_limma_voom_TMMwsp_orig_out_sig)
identical(crc_zeller_limma_voom_TMM_new_out_sig, crc_zeller_limma_voom_TMM_orig_out_sig)




################ t1d_mej and limma_voom (UNFILTERED) ##################

rm(list=ls(all.names=TRUE))
source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

t1d_mejialeon_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/t1d_mejialeon/t1d_mejialeon_ASVs_table.tsv",
                              header=TRUE, sep="\t", row.names=1, check.names=FALSE, comment.char="")

t1d_mejialeon_groups <- parse_two_meta_groups("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/t1d_mejialeon/t1d_mejialeon_metadata.tsv")

t1d_mejialeon_limma_voom_TMM_new_out <- limma_voom_two_group(table=t1d_mejialeon_abun,
                                                          group1=t1d_mejialeon_groups$group1,
                                                          group2=t1d_mejialeon_groups$group2,
                                                          norm_approach = "TMM")

t1d_mejialeon_limma_voom_TMMwsp_new_out <- limma_voom_two_group(table=t1d_mejialeon_abun,
                                                             group1=t1d_mejialeon_groups$group1,
                                                             group2=t1d_mejialeon_groups$group2,
                                                             norm_approach = "TMMwsp")


t1d_mejialeon_limma_voom_TMMwsp_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/t1d_mejialeon/No_filt_Results/Limma_voom_TMMwsp/limma_voom_tmmwsp_res.tsv",
                                                    header=TRUE, sep="\t", row.names=1, check.names=FALSE)

t1d_mejialeon_limma_voom_TMM_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/t1d_mejialeon/No_filt_Results/limma_voom_tmm_out/limma_voom_tmm_res.tsv",
                                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)

t1d_mejialeon_limma_voom_TMMwsp_orig_out_sig <- rownames(t1d_mejialeon_limma_voom_TMMwsp_orig_out)[which(t1d_mejialeon_limma_voom_TMMwsp_orig_out$adj.P.Val < 0.05)]
t1d_mejialeon_limma_voom_TMM_orig_out_sig <- rownames(t1d_mejialeon_limma_voom_TMM_orig_out)[which(t1d_mejialeon_limma_voom_TMM_orig_out$adj.P.Val < 0.05)]


t1d_mejialeon_limma_voom_TMMwsp_new_out_sig <- rownames(t1d_mejialeon_limma_voom_TMMwsp_new_out)[which(t1d_mejialeon_limma_voom_TMMwsp_new_out$adj.P.Val < 0.05)]
t1d_mejialeon_limma_voom_TMM_new_out_sig <- rownames(t1d_mejialeon_limma_voom_TMM_new_out)[which(t1d_mejialeon_limma_voom_TMM_new_out$adj.P.Val < 0.05)]

# Key check to see if results are the same:
identical(t1d_mejialeon_limma_voom_TMMwsp_new_out_sig, t1d_mejialeon_limma_voom_TMMwsp_orig_out_sig)
identical(t1d_mejialeon_limma_voom_TMM_new_out_sig, t1d_mejialeon_limma_voom_TMM_orig_out_sig)



################ MALL and limma_voom (filtered) ##################

rm(list=ls(all.names=TRUE))
source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

MALL_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/MALL/Fix_Results_0.1/fixed_non_rare_tables/MALL_ASVs_table.tsv",
                              header=TRUE, sep="\t", row.names=1, check.names=FALSE)

MALL_groups <- parse_two_meta_groups("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/MALL/MALL_metadata.tsv")

MALL_limma_voom_TMM_new_out <- limma_voom_two_group(table=MALL_abun,
                                                          group1=MALL_groups$group1,
                                                          group2=MALL_groups$group2,
                                                          norm_approach = "TMM")

MALL_limma_voom_TMMwsp_new_out <- limma_voom_two_group(table=MALL_abun,
                                                             group1=MALL_groups$group1,
                                                             group2=MALL_groups$group2,
                                                             norm_approach = "TMMwsp")


MALL_limma_voom_TMMwsp_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/MALL/Fix_Results_0.1/Limma_voom_TMMwsp/limma_voom_tmmwsp_res.tsv",
                                                    header=TRUE, sep="\t", row.names=1, check.names=FALSE)

MALL_limma_voom_TMM_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/MALL/Fix_Results_0.1/limma_voom_tmm_out/limma_voom_tmm_res.tsv",
                                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)

MALL_limma_voom_TMMwsp_orig_out_sig <- rownames(MALL_limma_voom_TMMwsp_orig_out)[which(MALL_limma_voom_TMMwsp_orig_out$adj.P.Val < 0.05)]
MALL_limma_voom_TMM_orig_out_sig <- rownames(MALL_limma_voom_TMM_orig_out)[which(MALL_limma_voom_TMM_orig_out$adj.P.Val < 0.05)]


MALL_limma_voom_TMMwsp_new_out_sig <- rownames(MALL_limma_voom_TMMwsp_new_out)[which(MALL_limma_voom_TMMwsp_new_out$adj.P.Val < 0.05)]
MALL_limma_voom_TMM_new_out_sig <- rownames(MALL_limma_voom_TMM_new_out)[which(MALL_limma_voom_TMM_new_out$adj.P.Val < 0.05)]

# Key check to see if results are the same:
identical(MALL_limma_voom_TMMwsp_new_out_sig, MALL_limma_voom_TMMwsp_orig_out_sig)
identical(MALL_limma_voom_TMM_new_out_sig, MALL_limma_voom_TMM_orig_out_sig)


################ glass_plastic_oberbeckmann and limma_voom (UNFILTERED) ##################

rm(list=ls(all.names=TRUE))
source("/home/gavin/projects/hackathon/sanity_checks/sanity_check_functions.R")

glass_plastic_oberbeckmann_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/glass_plastic_oberbeckmann/glass_plastic_oberbeckmann_ASVs_table.tsv",
                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE, comment.char="")

glass_plastic_oberbeckmann_groups <- parse_two_meta_groups("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/glass_plastic_oberbeckmann/glass_plastic_oberbeckmann_metadata.tsv")

glass_plastic_oberbeckmann_limma_voom_TMM_new_out <- limma_voom_two_group(table=glass_plastic_oberbeckmann_abun,
                                                             group1=glass_plastic_oberbeckmann_groups$group1,
                                                             group2=glass_plastic_oberbeckmann_groups$group2,
                                                             norm_approach = "TMM")

glass_plastic_oberbeckmann_limma_voom_TMMwsp_new_out <- limma_voom_two_group(table=glass_plastic_oberbeckmann_abun,
                                                                group1=glass_plastic_oberbeckmann_groups$group1,
                                                                group2=glass_plastic_oberbeckmann_groups$group2,
                                                                norm_approach = "TMMwsp")


glass_plastic_oberbeckmann_limma_voom_TMMwsp_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/glass_plastic_oberbeckmann/No_filt_Results/Limma_voom_TMMwsp/limma_voom_tmmwsp_res.tsv",
                                                       header=TRUE, sep="\t", row.names=1, check.names=FALSE)

glass_plastic_oberbeckmann_limma_voom_TMM_orig_out <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/glass_plastic_oberbeckmann/No_filt_Results/limma_voom_tmm_out/limma_voom_tmm_res.tsv",
                                                    header=TRUE, sep="\t", row.names=1, check.names=FALSE)

glass_plastic_oberbeckmann_limma_voom_TMMwsp_orig_out_sig <- rownames(glass_plastic_oberbeckmann_limma_voom_TMMwsp_orig_out)[which(glass_plastic_oberbeckmann_limma_voom_TMMwsp_orig_out$adj.P.Val < 0.05)]
glass_plastic_oberbeckmann_limma_voom_TMM_orig_out_sig <- rownames(glass_plastic_oberbeckmann_limma_voom_TMM_orig_out)[which(glass_plastic_oberbeckmann_limma_voom_TMM_orig_out$adj.P.Val < 0.05)]


glass_plastic_oberbeckmann_limma_voom_TMMwsp_new_out_sig <- rownames(glass_plastic_oberbeckmann_limma_voom_TMMwsp_new_out)[which(glass_plastic_oberbeckmann_limma_voom_TMMwsp_new_out$adj.P.Val < 0.05)]
glass_plastic_oberbeckmann_limma_voom_TMM_new_out_sig <- rownames(glass_plastic_oberbeckmann_limma_voom_TMM_new_out)[which(glass_plastic_oberbeckmann_limma_voom_TMM_new_out$adj.P.Val < 0.05)]

# Key check to see if results are the same:
identical(glass_plastic_oberbeckmann_limma_voom_TMMwsp_new_out_sig, glass_plastic_oberbeckmann_limma_voom_TMMwsp_orig_out_sig)
identical(glass_plastic_oberbeckmann_limma_voom_TMM_new_out_sig, glass_plastic_oberbeckmann_limma_voom_TMM_orig_out_sig)

