# EchoCH 
Here we present a quantitative framework for Clonal Hematopoiesis (CH) Analysis 

<img width="500" height="300" alt="image" src="https://github.com/user-attachments/assets/e603bc7e-538e-4f4f-b17b-ed74d2658164" />


### Installation

Install **devtools** using:
```
install.packages("devtools",dependencies = T)
```
Once **devtools** is installed, run the following to install **EchoCH**:
```
library(devtools)
devtools::install_github("marabzadeh/EchoCH", force = T)
```
## Instructions for Use

EchoCH_test.R includes guides to how to use the method and read the output. 
For the input file we need the list of curated CH mutations with the requested file format. 
Output files are being produced accordingly. 

```
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
```

## Authors
* **Mona Arabzadeh**
* More details on the method and the data that has been analyzed in the manuscript under preparation. 

## Acknowledgments
* [Khiabanian Lab](https://khiabanian-lab.org)
* [Gillis Lab](https://www.moffitt.org/research-science/researchers/nancy-gillis-johnson/)
