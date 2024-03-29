# GBM-Radiogenomics
Biological pathways underlying prognostic radiomics phenotypes from paired MRI and RNA sequencing in glioblastoma, accepted by Radiology


--------------------------------------------------------------------------------
R Script Introduction

Run radiomics_main.R for radiomics analysis, including feature selection, signature building, signature validation, Kaplan-Meier analysis, calibration curve calculation, and decision curve analysis. Run radiogenomics_main.R for radiogenomics analysis, including WGCNA analysis, pathway enrichment analysis, radiomics-correlated driving pathway identification, radiomics-correlation key genes identification, prognostic significance test of key genes, and replication analysis on external test set. All data required to run the R script should be stored in the .csv files in the “input” folder. The detailed introductions for each .csv file are listed in the following.

Train.csv: the data for patients of the radiomics training subset. The column names from left to right indicate the patient number, survival status (0 or 1), overall survival (in days or months), sex (0 or 1), age (years), KPS, IDH1 mutation status (0 or 1), radiation therapy yes or no (0 or 1), chemotherapy yes or no (0 or 1), extent of surgical resection complete or incomplete (0 or 1), and the radiomics feature values. The column orders from “Status” to “SurgicalResection” should not be changed. X1…Xn are radiomics features. Users can modify the number of features.

Test.csv: the data for patients of the radiomics validation subset.

gmt folder: databases of annotated pathways. Do not change them.

externalTestData.csv: RNA-seq data of the external test set.

externalTestRadiomics.csv: overall survival, survival status, radiomics signature, and each of the selected features constituting the signature for patients in the external test set. In our study there were 13 features used to build the signature. Users can modify them based on their own study data.

module.csv: IDs of the finally selected pathways in gene modules. Class is not generated by the program but should be determined by the user according to the biological meaning of each pathway.

radioGenTrainData.csv: RNA-seq data of the radiogenomics training set.

radioGenTrainRadiomics.csv: overall survival, survival status, radiomics signature, and each of the selected features constituting the signature for patients in the radiogenomics training set.

tcgaTestclinical.csv: overall survival and survival status for patients in the TCGA test set.

tcgaTestData.csv: RNA-seq data of the TCGA test set.
 
