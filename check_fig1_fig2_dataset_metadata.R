### Check that values of dataset metadata used in Fig 1 / 2 are correct.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/projects/hackathon/sanity_checks/fig1_fig2_intermediate_objects/")

dataset2name <- read.table("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/mapfiles/dataset_name_mapping.csv",
                           header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

fixed_hackathon_metadata_filt <- readRDS("fixed_hackathon_metadata_filt.rds")
fixed_hackathon_metadata_unfilt <- readRDS("fixed_hackathon_metadata_unfilt.rds")

permanova_results <- read.table("/home/gavin/github_repos/hackathon/Comparison_of_DA_microbiome_methods/Misc_datafiles/aitchison_permanova_results.tsv",
                                header = TRUE, sep = "\t", row.names = 1)

### Blueberry unfilt

Blueberry_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Blueberry/No_filt_Results/fixed_non_rare_tables/Blueberry_ASVs_table.tsv",
                                        header=TRUE, sep ="\t", row.names = 1)

# Check sparsity
fixed_hackathon_metadata_unfilt["Soil - Blueberry", "Sparsity"] == length(which(Blueberry_unfilt_ASV_abun == 0)) / (dim(Blueberry_unfilt_ASV_abun)[1] * dim(Blueberry_unfilt_ASV_abun)[2])

# Check median richness
fixed_hackathon_metadata_unfilt["Soil - Blueberry", "Richness"] == median(colSums(Blueberry_unfilt_ASV_abun > 0))

# Check log(Sample size)
fixed_hackathon_metadata_unfilt["Soil - Blueberry", "log(Sample size)"] == log(ncol(Blueberry_unfilt_ASV_abun))

# Check read depth range
fixed_hackathon_metadata_unfilt["Soil - Blueberry", "log(Read depth range)"] == log(max(colSums(Blueberry_unfilt_ASV_abun)) - min(colSums(Blueberry_unfilt_ASV_abun)))

# Check read depth CV
fixed_hackathon_metadata_unfilt["Soil - Blueberry", "Read depth variation"]  == sd(colSums(Blueberry_unfilt_ASV_abun)) / mean(colSums(Blueberry_unfilt_ASV_abun))

# Check permanova R-squared
fixed_hackathon_metadata_unfilt["Soil - Blueberry", "log(Aitch. dist. effect size)"] == log(permanova_results["Blueberry", "unfilt_nonrare_R2"])



### glass_plastic_oberbeckmann unfilt

glass_plastic_oberbeckmann_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/glass_plastic_oberbeckmann/No_filt_Results/fixed_non_rare_tables/glass_plastic_oberbeckmann_ASVs_table.tsv",
                                        header=TRUE, sep ="\t", row.names = 1)

# Check sparsity
fixed_hackathon_metadata_unfilt[dataset2name["glass_plastic_oberbeckmann", "clean"], "Sparsity"] == length(which(glass_plastic_oberbeckmann_unfilt_ASV_abun == 0)) / (dim(glass_plastic_oberbeckmann_unfilt_ASV_abun)[1] * dim(glass_plastic_oberbeckmann_unfilt_ASV_abun)[2])

# Check median richness
fixed_hackathon_metadata_unfilt[dataset2name["glass_plastic_oberbeckmann", "clean"], "Richness"] == median(colSums(glass_plastic_oberbeckmann_unfilt_ASV_abun > 0))

# Check log(Sample size)
fixed_hackathon_metadata_unfilt[dataset2name["glass_plastic_oberbeckmann", "clean"], "log(Sample size)"] == log(ncol(glass_plastic_oberbeckmann_unfilt_ASV_abun))

# Check read depth range
fixed_hackathon_metadata_unfilt[dataset2name["glass_plastic_oberbeckmann", "clean"], "log(Read depth range)"] == log(max(colSums(glass_plastic_oberbeckmann_unfilt_ASV_abun)) - min(colSums(glass_plastic_oberbeckmann_unfilt_ASV_abun)))

# Check read depth CV
fixed_hackathon_metadata_unfilt[dataset2name["glass_plastic_oberbeckmann", "clean"], "Read depth variation"]  == sd(colSums(glass_plastic_oberbeckmann_unfilt_ASV_abun)) / mean(colSums(glass_plastic_oberbeckmann_unfilt_ASV_abun))

# Check permanova R-squared
fixed_hackathon_metadata_unfilt[dataset2name["glass_plastic_oberbeckmann", "clean"], "log(Aitch. dist. effect size)"] == log(permanova_results["glass_plastic_oberbeckmann", "unfilt_nonrare_R2"])



### ob_turnbaugh unfilt

ob_turnbaugh_unfilt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ob_turnbaugh/No_filt_Results/fixed_non_rare_tables/ob_turnbaugh_ASVs_table.tsv",
                                                         header=TRUE, sep ="\t", row.names = 1)

# Check sparsity
fixed_hackathon_metadata_unfilt[dataset2name["ob_turnbaugh", "clean"], "Sparsity"] == length(which(ob_turnbaugh_unfilt_ASV_abun == 0)) / (dim(ob_turnbaugh_unfilt_ASV_abun)[1] * dim(ob_turnbaugh_unfilt_ASV_abun)[2])

# Check median richness
fixed_hackathon_metadata_unfilt[dataset2name["ob_turnbaugh", "clean"], "Richness"] == median(colSums(ob_turnbaugh_unfilt_ASV_abun > 0))

# Check log(Sample size)
fixed_hackathon_metadata_unfilt[dataset2name["ob_turnbaugh", "clean"], "log(Sample size)"] == log(ncol(ob_turnbaugh_unfilt_ASV_abun))

# Check read depth range
fixed_hackathon_metadata_unfilt[dataset2name["ob_turnbaugh", "clean"], "log(Read depth range)"] == log(max(colSums(ob_turnbaugh_unfilt_ASV_abun)) - min(colSums(ob_turnbaugh_unfilt_ASV_abun)))

# Check read depth CV
fixed_hackathon_metadata_unfilt[dataset2name["ob_turnbaugh", "clean"], "Read depth variation"]  == sd(colSums(ob_turnbaugh_unfilt_ASV_abun)) / mean(colSums(ob_turnbaugh_unfilt_ASV_abun))

# Check permanova R-squared
fixed_hackathon_metadata_unfilt[dataset2name["ob_turnbaugh", "clean"], "log(Aitch. dist. effect size)"] == log(permanova_results["ob_turnbaugh", "unfilt_nonrare_R2"])















### Blueberry filt

Blueberry_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/Blueberry/Fix_Results_0.1/fixed_non_rare_tables/Blueberry_ASVs_table.tsv",
                                        header=TRUE, sep ="\t", row.names = 1)

# Check sparsity
fixed_hackathon_metadata_filt["Soil - Blueberry", "Sparsity"] == length(which(Blueberry_filt_ASV_abun == 0)) / (dim(Blueberry_filt_ASV_abun)[1] * dim(Blueberry_filt_ASV_abun)[2])

# Check median richness
fixed_hackathon_metadata_filt["Soil - Blueberry", "Richness"] == median(colSums(Blueberry_filt_ASV_abun > 0))

# Check log(Sample size)
fixed_hackathon_metadata_filt["Soil - Blueberry", "log(Sample size)"] == log(ncol(Blueberry_filt_ASV_abun))

# Check read depth range
fixed_hackathon_metadata_filt["Soil - Blueberry", "log(Read depth range)"] == log(max(colSums(Blueberry_filt_ASV_abun)) - min(colSums(Blueberry_filt_ASV_abun)))

# Check read depth CV
fixed_hackathon_metadata_filt["Soil - Blueberry", "Read depth variation"]  == sd(colSums(Blueberry_filt_ASV_abun)) / mean(colSums(Blueberry_filt_ASV_abun))

# Check permanova R-squared
fixed_hackathon_metadata_filt["Soil - Blueberry", "log(Aitch. dist. effect size)"] == log(permanova_results["Blueberry", "filt_nonrare_R2"])



### glass_plastic_oberbeckmann filt

glass_plastic_oberbeckmann_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/glass_plastic_oberbeckmann/Fix_Results_0.1/fixed_non_rare_tables/glass_plastic_oberbeckmann_ASVs_table.tsv",
                                                         header=TRUE, sep ="\t", row.names = 1)

# Check sparsity
fixed_hackathon_metadata_filt[dataset2name["glass_plastic_oberbeckmann", "clean"], "Sparsity"] == length(which(glass_plastic_oberbeckmann_filt_ASV_abun == 0)) / (dim(glass_plastic_oberbeckmann_filt_ASV_abun)[1] * dim(glass_plastic_oberbeckmann_filt_ASV_abun)[2])

# Check median richness
fixed_hackathon_metadata_filt[dataset2name["glass_plastic_oberbeckmann", "clean"], "Richness"] == median(colSums(glass_plastic_oberbeckmann_filt_ASV_abun > 0))

# Check log(Sample size)
fixed_hackathon_metadata_filt[dataset2name["glass_plastic_oberbeckmann", "clean"], "log(Sample size)"] == log(ncol(glass_plastic_oberbeckmann_filt_ASV_abun))

# Check read depth range
fixed_hackathon_metadata_filt[dataset2name["glass_plastic_oberbeckmann", "clean"], "log(Read depth range)"] == log(max(colSums(glass_plastic_oberbeckmann_filt_ASV_abun)) - min(colSums(glass_plastic_oberbeckmann_filt_ASV_abun)))

# Check read depth CV
fixed_hackathon_metadata_filt[dataset2name["glass_plastic_oberbeckmann", "clean"], "Read depth variation"]  == sd(colSums(glass_plastic_oberbeckmann_filt_ASV_abun)) / mean(colSums(glass_plastic_oberbeckmann_filt_ASV_abun))

# Check permanova R-squared
fixed_hackathon_metadata_filt[dataset2name["glass_plastic_oberbeckmann", "clean"], "log(Aitch. dist. effect size)"] == log(permanova_results["glass_plastic_oberbeckmann", "filt_nonrare_R2"])



### ob_turnbaugh filt

ob_turnbaugh_filt_ASV_abun <- read.table("/home/jacob/projects/HACKATHON_ANCOM_FIX_21_03_13/Hackathon/Studies/ob_turnbaugh/Fix_Results_0.1/fixed_non_rare_tables/ob_turnbaugh_ASVs_table.tsv",
                                           header=TRUE, sep ="\t", row.names = 1)

# Check sparsity
fixed_hackathon_metadata_filt[dataset2name["ob_turnbaugh", "clean"], "Sparsity"] == length(which(ob_turnbaugh_filt_ASV_abun == 0)) / (dim(ob_turnbaugh_filt_ASV_abun)[1] * dim(ob_turnbaugh_filt_ASV_abun)[2])

# Check median richness
fixed_hackathon_metadata_filt[dataset2name["ob_turnbaugh", "clean"], "Richness"] == median(colSums(ob_turnbaugh_filt_ASV_abun > 0))

# Check log(Sample size)
fixed_hackathon_metadata_filt[dataset2name["ob_turnbaugh", "clean"], "log(Sample size)"] == log(ncol(ob_turnbaugh_filt_ASV_abun))

# Check read depth range
fixed_hackathon_metadata_filt[dataset2name["ob_turnbaugh", "clean"], "log(Read depth range)"] == log(max(colSums(ob_turnbaugh_filt_ASV_abun)) - min(colSums(ob_turnbaugh_filt_ASV_abun)))

# Check read depth CV
fixed_hackathon_metadata_filt[dataset2name["ob_turnbaugh", "clean"], "Read depth variation"]  == sd(colSums(ob_turnbaugh_filt_ASV_abun)) / mean(colSums(ob_turnbaugh_filt_ASV_abun))

# Check permanova R-squared
fixed_hackathon_metadata_filt[dataset2name["ob_turnbaugh", "clean"], "log(Aitch. dist. effect size)"] == log(permanova_results["ob_turnbaugh", "filt_nonrare_R2"])

