### Check some of the values displayed in Figure 1

rm(list = ls(all.names = TRUE))

dataset2name <- read.table("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/mapfiles/dataset_name_mapping.csv",
                           header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

tool2name <- read.table("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/mapfiles/tool_name_mapping.csv",
                           header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)


unfilt_sig_percent_fix <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/unfilt_sig_percent_fix.rds")
filt_sig_percent_fix <- readRDS(file = "/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/filt_sig_percent_fix.rds")

# Check five random dataset / tool combinations.
# Check them separately as unfilt and filt

# "Ji_WTP_DS"         "hiv_lozupone"      "ArcticFireSoils"   "hiv_noguerajulian" "Exercise"
# "wilcoxonrare"  "metagenomeSeq" "lefse"         "maaslin2"      "edger"   

# Ji_WTP_DS and wilcoxonrare - unfilt

Ji_WTP_DS_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Ji_WTP_DS//No_filt_Results/fixed_rare_tables/Ji_WTP_DS_ASVs_table.tsv",
                                        header = TRUE, sep = "\t", row.names = 1)

Ji_WTP_DS_unfilt_wilcoxonrare <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Ji_WTP_DS/No_filt_Results/Wilcoxon_rare_out/Wil_rare_results.tsv",
                                            header = TRUE, sep = "\t", row.names = 1)

Ji_WTP_DS_unfilt_wilcoxonrare_sig <- length(which(p.adjust(Ji_WTP_DS_unfilt_wilcoxonrare$x, "fdr") < 0.05))

# Check that this value is shown in the figure
dataset2name["Ji_WTP_DS", "clean"]
Ji_WTP_DS_unfilt_wilcoxonrare_sig

# Check that percent is correct too:
(Ji_WTP_DS_unfilt_wilcoxonrare_sig / nrow(Ji_WTP_DS_unfilt_ASV_abun)) * 100 == unfilt_sig_percent_fix[dataset2name["Ji_WTP_DS", "clean"], "wilcoxonrare"]




# hiv_lozupone and metagenomeSeq - unfilt

hiv_lozupone_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_lozupone//No_filt_Results/fixed_non_rare_tables/hiv_lozupone_ASVs_table.tsv",
                                        header = TRUE, sep = "\t", row.names = 1)

hiv_lozupone_unfilt_metagenomeSeq <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_lozupone/No_filt_Results/metagenomeSeq_out/mgSeq_res.tsv",
                                            header = TRUE, sep = "\t", row.names = 1)

hiv_lozupone_unfilt_metagenomeSeq_sig <- length(which(hiv_lozupone_unfilt_metagenomeSeq$adjPvalues < 0.05))

# Check that this value is shown in the figure
dataset2name["hiv_lozupone", "clean"]
hiv_lozupone_unfilt_metagenomeSeq_sig

# Check that percent is correct too:
(hiv_lozupone_unfilt_metagenomeSeq_sig / nrow(hiv_lozupone_unfilt_ASV_abun)) * 100 == unfilt_sig_percent_fix[dataset2name["hiv_lozupone", "clean"], "metagenomeSeq"]




# ArcticFireSoils and lefse - unfilt

ArcticFireSoils_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFireSoils//No_filt_Results/fixed_rare_tables/ArcticFireSoils_ASVs_table.tsv",
                                           header = TRUE, sep = "\t", row.names = 1)

ArcticFireSoils_unfilt_lefse <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFireSoils/No_filt_Results/Lefse_out/Lefse_results.tsv",
                                                header = FALSE, sep = "\t", row.names = 1)

ArcticFireSoils_unfilt_lefse <- ArcticFireSoils_unfilt_lefse[-which(ArcticFireSoils_unfilt_lefse$V5 == "-"), ]

ArcticFireSoils_unfilt_lefse$V5 <- as.numeric(as.character(ArcticFireSoils_unfilt_lefse$V5))

ArcticFireSoils_unfilt_lefse_sig <- length(which(ArcticFireSoils_unfilt_lefse$V5 < 0.05))

# Check that this value is shown in the figure
dataset2name["ArcticFireSoils", "clean"]
ArcticFireSoils_unfilt_lefse_sig

# Check that percent is correct too:
(ArcticFireSoils_unfilt_lefse_sig / nrow(ArcticFireSoils_unfilt_ASV_abun)) * 100 == unfilt_sig_percent_fix[dataset2name["ArcticFireSoils", "clean"], "lefse"]




# hiv_noguerajulian and maaslin2 - unfilt

hiv_noguerajulian_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian//No_filt_Results/fixed_non_rare_tables//hiv_noguerajulian_ASVs_table.tsv",
                                              header = TRUE, sep = "\t", row.names = 1)

hiv_noguerajulian_unfilt_maaslin2 <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian/No_filt_Results/Maaslin2_out/significant_results.tsv",
                                           header = TRUE, sep = "\t")

hiv_noguerajulian_unfilt_maaslin2_sig <- length(which(hiv_noguerajulian_unfilt_maaslin2$qval < 0.05))

# Check that this value is shown in the figure
dataset2name["hiv_noguerajulian", "clean"]
hiv_noguerajulian_unfilt_maaslin2_sig

# Check that percent is correct too:
(hiv_noguerajulian_unfilt_maaslin2_sig / nrow(hiv_noguerajulian_unfilt_ASV_abun)) * 100 == unfilt_sig_percent_fix[dataset2name["hiv_noguerajulian", "clean"], "maaslin2"]





# Exercise and edger - unfilt

Exercise_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Exercise//No_filt_Results/fixed_non_rare_tables/Exercise_ASVs_table.tsv",
                                                header = TRUE, sep = "\t", row.names = 1)

Exercise_unfilt_edger <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Exercise/No_filt_Results/edgeR_out/edgeR_res.tsv",
                                                header = TRUE, sep = "\t")

Exercise_unfilt_edger_sig <- length(which(Exercise_unfilt_edger$FDR < 0.05))

# Check that this value is shown in the figure
dataset2name["Exercise", "clean"]
Exercise_unfilt_edger_sig

# Check that percent is correct too:
(Exercise_unfilt_edger_sig / nrow(Exercise_unfilt_ASV_abun)) * 100 == unfilt_sig_percent_fix[dataset2name["Exercise", "clean"], "edger"]



# Do a limma voom test just for fun

ArcticFreshwaters_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFreshwaters/No_filt_Results/fixed_non_rare_tables/ArcticFreshwaters_ASVs_table.tsv",
                                       header = TRUE, sep = "\t", row.names = 1)

ArcticFreshwaters_unfilt_limma_voom_TMMwsp <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFreshwaters/No_filt_Results/Limma_voom_TMMwsp/limma_voom_tmmwsp_res.tsv",
                                    header = TRUE, sep = "\t")

ArcticFreshwaters_unfilt_limma_voom_TMMwsp_sig <- length(which(ArcticFreshwaters_unfilt_limma_voom_TMMwsp$adj.P.Val < 0.05))

# Check that this value is shown in the figure
dataset2name["ArcticFreshwaters", "clean"]
ArcticFreshwaters_unfilt_limma_voom_TMMwsp_sig

# Check that percent is correct too:
(ArcticFreshwaters_unfilt_limma_voom_TMMwsp_sig / nrow(ArcticFreshwaters_unfilt_ASV_abun)) * 100 == unfilt_sig_percent_fix[dataset2name["ArcticFreshwaters", "clean"], "limma_voom_TMMwsp"]












########### As above, but look at filtered data instead.


# Ji_WTP_DS and wilcoxonrare - filt

Ji_WTP_DS_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Ji_WTP_DS//Fix_Results_0.1//fixed_rare_tables/Ji_WTP_DS_ASVs_table.tsv",
                                        header = TRUE, sep = "\t", row.names = 1)

Ji_WTP_DS_filt_wilcoxonrare <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Ji_WTP_DS/Fix_Results_0.1/Wilcoxon_rare_out/Wil_rare_results.tsv",
                                            header = TRUE, sep = "\t", row.names = 1)

Ji_WTP_DS_filt_wilcoxonrare_sig <- length(which(p.adjust(Ji_WTP_DS_filt_wilcoxonrare$x, "fdr") < 0.05))

# Check that this value is shown in the figure
dataset2name["Ji_WTP_DS", "clean"]
Ji_WTP_DS_filt_wilcoxonrare_sig

# Check that percent is correct too:
(Ji_WTP_DS_filt_wilcoxonrare_sig / nrow(Ji_WTP_DS_filt_ASV_abun)) * 100 == filt_sig_percent_fix[dataset2name["Ji_WTP_DS", "clean"], "wilcoxonrare"]




# hiv_lozupone and metagenomeSeq - filt

hiv_lozupone_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_lozupone//Fix_Results_0.1/fixed_non_rare_tables/hiv_lozupone_ASVs_table.tsv",
                                           header = TRUE, sep = "\t", row.names = 1)

hiv_lozupone_filt_metagenomeSeq <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_lozupone/Fix_Results_0.1/metagenomeSeq_out/mgSeq_res.tsv",
                                                header = TRUE, sep = "\t", row.names = 1)

hiv_lozupone_filt_metagenomeSeq_sig <- length(which(hiv_lozupone_filt_metagenomeSeq$adjPvalues < 0.05))

# Check that this value is shown in the figure
dataset2name["hiv_lozupone", "clean"]
hiv_lozupone_filt_metagenomeSeq_sig

# Check that percent is correct too:
(hiv_lozupone_filt_metagenomeSeq_sig / nrow(hiv_lozupone_filt_ASV_abun)) * 100 == filt_sig_percent_fix[dataset2name["hiv_lozupone", "clean"], "metagenomeSeq"]




# ArcticFireSoils and lefse - filt

ArcticFireSoils_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFireSoils//Fix_Results_0.1/fixed_rare_tables/ArcticFireSoils_ASVs_table.tsv",
                                              header = TRUE, sep = "\t", row.names = 1)

ArcticFireSoils_filt_lefse <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFireSoils/Fix_Results_0.1/Lefse_out/Lefse_results.tsv",
                                           header = FALSE, sep = "\t", row.names = 1)

ArcticFireSoils_filt_lefse <- ArcticFireSoils_filt_lefse[-which(ArcticFireSoils_filt_lefse$V5 == "-"), ]

ArcticFireSoils_filt_lefse$V5 <- as.numeric(as.character(ArcticFireSoils_filt_lefse$V5))

ArcticFireSoils_filt_lefse_sig <- length(which(ArcticFireSoils_filt_lefse$V5 < 0.05))

# Check that this value is shown in the figure
dataset2name["ArcticFireSoils", "clean"]
ArcticFireSoils_filt_lefse_sig

# Check that percent is correct too:
(ArcticFireSoils_filt_lefse_sig / nrow(ArcticFireSoils_filt_ASV_abun)) * 100 == filt_sig_percent_fix[dataset2name["ArcticFireSoils", "clean"], "lefse"]




# hiv_noguerajulian and maaslin2 - filt

hiv_noguerajulian_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian//Fix_Results_0.1/fixed_non_rare_tables//hiv_noguerajulian_ASVs_table.tsv",
                                                header = TRUE, sep = "\t", row.names = 1)

hiv_noguerajulian_filt_maaslin2 <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/hiv_noguerajulian/Fix_Results_0.1/Maaslin2_out/significant_results.tsv",
                                                header = TRUE, sep = "\t")

hiv_noguerajulian_filt_maaslin2_sig <- length(which(hiv_noguerajulian_filt_maaslin2$qval < 0.05))

# Check that this value is shown in the figure
dataset2name["hiv_noguerajulian", "clean"]
hiv_noguerajulian_filt_maaslin2_sig

# Check that percent is correct too:
(hiv_noguerajulian_filt_maaslin2_sig / nrow(hiv_noguerajulian_filt_ASV_abun)) * 100 == filt_sig_percent_fix[dataset2name["hiv_noguerajulian", "clean"], "maaslin2"]





# Exercise and edger - filt

Exercise_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Exercise//Fix_Results_0.1/fixed_non_rare_tables/Exercise_ASVs_table.tsv",
                                       header = TRUE, sep = "\t", row.names = 1)

Exercise_filt_edger <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Exercise/Fix_Results_0.1/edgeR_out/edgeR_res.tsv",
                                    header = TRUE, sep = "\t")

Exercise_filt_edger_sig <- length(which(Exercise_filt_edger$FDR < 0.05))

# Check that this value is shown in the figure
dataset2name["Exercise", "clean"]
Exercise_filt_edger_sig

# Check that percent is correct too:
(Exercise_filt_edger_sig / nrow(Exercise_filt_ASV_abun)) * 100 == filt_sig_percent_fix[dataset2name["Exercise", "clean"], "edger"]



# Do a limma voom test just for fun

ArcticFreshwaters_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFreshwaters/Fix_Results_0.1/fixed_non_rare_tables/ArcticFreshwaters_ASVs_table.tsv",
                                                header = TRUE, sep = "\t", row.names = 1)

ArcticFreshwaters_filt_limma_voom_TMMwsp <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ArcticFreshwaters/Fix_Results_0.1/Limma_voom_TMMwsp/limma_voom_tmmwsp_res.tsv",
                                                         header = TRUE, sep = "\t")

ArcticFreshwaters_filt_limma_voom_TMMwsp_sig <- length(which(ArcticFreshwaters_filt_limma_voom_TMMwsp$adj.P.Val < 0.05))

# Check that this value is shown in the figure
dataset2name["ArcticFreshwaters", "clean"]
ArcticFreshwaters_filt_limma_voom_TMMwsp_sig

# Check that percent is correct too:
(ArcticFreshwaters_filt_limma_voom_TMMwsp_sig / nrow(ArcticFreshwaters_filt_ASV_abun)) * 100 == filt_sig_percent_fix[dataset2name["ArcticFreshwaters", "clean"], "limma_voom_TMMwsp"]


