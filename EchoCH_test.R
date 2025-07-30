## Mona Arabzadeh -- July 29th 2025
## Test file for EchoCH: quantitative framework for Clonal Hematopoiesis (CH) Statistical Analysis

library(devtools)
devtools::install_github("marabzadeh/EchoCH", force = T)

library(EchoCH)

# Example usage
pathRead <- "~/inps/"
CH.df <- readxl::read_xlsx( paste0(pathRead, "TableS2.xlsx" ) )
head(CH.df)

#CH.df a data.frame with CH mutations with the following columns::
#'                                              patient_id (combination of number and letter), gene, treatment, amino_acid_change,
#'                                              depth, vaf, temporal_order, pre/post_treatment,
#'                                             days_from_treatment_start_to_sample_collection )


# Example usage

# Growth analysis per-variant per-patient for each time-point
PerVarGrowth.df <- EchoCH::GrowthAnalysisPerVariant(CH.df)

# Growth analysis per-gene per-treatment
PerGeneGrowth.df <- EchoCH::GrowthAnalysisPerGene(PerVarGrowth.df, time_threshold = 540)

# Growth analysis per-treatment
PerTreatmentGrowth.df <- EchoCH::GrowthAnalysisPerTreatment(PerVarGrowth.df, time_threshold = 540)


# Calculate Neff for each patient with the Fitness likelihood for each variant in that patient
Neff_Pr_List <- EchoCH::NeffAnalysis(CH.df, mutNumber_threshold = 1, time_threshold = 540)
NeffList.df <- Neff_Pr_List[[1]]

#Calculate Fitness likelihood for each variant in each patient
prList.df <- EchoCH::FitnessLikelihoodAnalysis(CH.df, Neff_Pr_List, SigOrderOfMagnitude_threshold = 5, time_threshold = 540)
prList.df <- prList.df [ prList.df$days == "<=540", ]


# Write outputs
pathWrite <- "~/outs/"
EchoCH::WriteFunction (input.df=PerVarGrowth.df, path = pathWrite, name = "PerVarGrowth" )
EchoCH::WriteFunction (input.df=PerGeneGrowth.df, path = pathWrite, name = "PerGeneGrowth" )
EchoCH::WriteFunction (input.df=PerTreatmentGrowth.df, path = pathWrite, name = "PerTreatmentGrowth" )
EchoCH::WriteFunction (input.df=NeffList.df, path = pathWrite, name = "NeffList" )
EchoCH::WriteFunction (input.df=prList.df, path = pathWrite, name = "prList" )
