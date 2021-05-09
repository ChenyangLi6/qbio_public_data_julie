if(!requireNamespace("TCGAbiolinks", quietly = TRUE))
  install.packages("TCGAbiolinks")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
BiocManager::install("maftools")
library(maftools)

if(!requireNamespace("SummarizedExperiment", quietly = TRUE))
  install.packages("SummarizedExperiment")
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")
library(survival)
library(survminer)
library(arsenal)
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))
colnames(clinic)[1] <- "Tumor_Sample_Barcode"
my_maf <- data.table::fread("GDCdata/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.csv")
maf_dataframe = read.maf(maf = my_maf, clinicalData = clinic, isTCGA = TRUE)

young_samples <- maf_dataframe@clinical.data[age_category == "Young", Tumor_Sample_Barcode]
mutation_pairs_young <- subsetMaf(maf_dataframe, tsb=young_samples )
old_samples <- maf_dataframe@clinical.data[age_category == "Old", Tumor_Sample_Barcode]
mutation_pairs_old <- subsetMaf(maf_dataframe, tsb=old_samples )
pdf( "Oncoplots_result.pdf")
coOncoplot(mutation_pairs_old, mutation_pairs_young, m1Name="Old", m2Name="Young",removeNonMutated = FALSE,gene_mar= 3.3,titleFontSize = 0.8)
coOncoplot(mutation_pairs_young, mutation_pairs_old, m1Name="Young", m2Name="Old", removeNonMutated = FALSE,gene_mar= 0.8,titleFontSize = 0.8,genes = c("TP53", "PIK3CA","GATA3","TTN","APOB","HYDIN","RYR2"))
dev.off()