if(!requireNamespace("TCGAbiolinks", quietly = TRUE))
  install.packages("TCGAbiolinks")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

if(!requireNamespace("SummarizedExperiment", quietly = TRUE))
  install.packages("SummarizedExperiment")
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

library(TCGAbiolinks)
library(SummarizedExperiment)

barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
                     "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
                     "TCGA-A2-A0T3-01A-====21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)
patient_data$tumor_stage[is.na(patient_data$tumor_stage)] <- "missing"
patient_data<-patient_data[!(patient_data$tumor_stage=="stage iv"),]
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
patient_data$age_category = ifelse(patient_ages < 40, "Young", ifelse(patient_ages >= 60, "Old", "Mid"))

htseq_counts <- assays(sum_exp)$"HTSeq - Counts"

patient_data$MAP3K1_counts = htseq_counts["ENSG00000228650",]

patient_data$MAP3K1_counts_log = sapply(htseq_counts["ENSG00000228650",], log10)

boxplot(MAP3K1_counts_log~age_category, main = "Boxplot of HTSeq - Counts for MAP3K1 by Age Category", data = patient_data)
pdf("box_plots.pdf")
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$MAP3K1_counts)
dev.off()