#jonathan code

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!requireNamespace("devtools")) BiocManager::install(c("devtools"))
if(!requireNamespace("robustbase"))BiocManager::install(c("robustbase"))
library(devtools)
library(robustbase)
if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
if(!requireNamespace("survival"))install.packages(c("survival"))

#load required libraries
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(survival)

#obtain MAF data 
mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="mutect2")
maf_dataframe = read.maf(mutation)

#query and obtain HT-SEQ data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query) #only need this line of code ONCE to download the data
sum_exp <- GDCprepare(query)


patient_data <- colData(sum_exp)
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
#here we are creating a NEW COLUMN in patient_data called "age_category"
#NOTE: This will NOT be added to colData(sum_exp). Instead it will only be added to the patient_data data table.

#create new column with age categories for all patient samples
patient_data$age_category = ifelse(patient_ages < 40, "Young", ifelse(patient_ages >= 60, "Old", "Mid"))
#The ifelse() form is: ifelse( condition, action when condition is true, action when condition is false ). Here we have two ifelse() embedded together

patient_data <- patient_data[, c("barcode", "age_category", "paper_age_at_initial_pathologic_diagnosis")]


#get shortened patient barcodes so that we can compare with
short_maf <- substr(maf_dataframe@clinical.data$Tumor_Sample_Barcode, 1,12)

#create a new column in maf_dataframe
maf_dataframe@clinical.data$short_barcodes <- short_maf

#extract age category information for each barcode that we have
maf_ages <- patient_data[short_maf, "age_category"]

maf_dataframe@clinical.data$Ages <- maf_ages

#Extract codes for each age group
young_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Young",]

mid_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Mid",]

old_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Old",]


#count how many patient samples do not have age associated with them
NA_codes <- maf_dataframe@clinical.data[is.na(maf_dataframe@clinical.data$Ages),]
print(dim(NA_codes))

#create maf subsets for each age group
young_maf <- subsetMaf(maf_dataframe, tsb = young_codes$Tumor_Sample_Barcode)

mid_maf <- subsetMaf(maf_dataframe, tsb = mid_codes$Tumor_Sample_Barcode)

old_maf <- subsetMaf(maf_dataframe, tsb = old_codes$Tumor_Sample_Barcode)


#get oncoplots
oncoplot(maf = maf_dataframe, draw_titv = TRUE, top = 3)
oncoplot(maf = young_maf, draw_titv = TRUE, top = 3)
oncoplot(maf = mid_maf, draw_titv = TRUE, top = 3)
oncoplot(maf = old_maf, draw_titv = TRUE, top = 3)


#get boxplots for top 3 genes for each age group

htseq_counts <- assays(sum_exp)$"HTSeq - Counts"

# #Genes of interest: TP53, PIK3CA, TTN, GATA3
# #external_gene_name does NOT need quotation marks because there are no spaces in the name of the column
# obtain IDs for genes using Boolean indexing
TP53_mask <- rowData(sum_exp)$external_gene_name == "TP53"
TP53_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TP53_mask]

PIK3CA_mask <- rowData(sum_exp)$external_gene_name == "PIK3CA"
PIK3CA_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[PIK3CA_mask]

TTN_mask <- rowData(sum_exp)$external_gene_name == "TTN"
TTN_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[TTN_mask]

GATA3_mask <- rowData(sum_exp)$external_gene_name == "GATA3"
GATA3_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[GATA3_mask]

# #get HTseq counts for desired genes
patient_data$TP53_counts_log = sapply(htseq_counts[TP53_ENSG_name,], log10)
patient_data$PIK3CA_counts_log = sapply(htseq_counts[PIK3CA_ENSG_name,], log10)
patient_data$TTN_counts_log = sapply(htseq_counts[TTN_ENSG_name,], log10)
patient_data$GATA3_counts_log = sapply(htseq_counts[GATA3_ENSG_name,], log10)


#exclude mid data
patient_data <-na.omit(patient_data)
patient_data <- patient_data[patient_data$age_category != "Mid",]

# 
# # Boxplots of HTseq counts by age
jpeg("TP53_log_counts_by_age.jpg")
boxplot(TP53_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for TP53 by Age Category",
        xlab = "Age Category", ylab = "TP53 Log Counts")
dev.off()

jpeg("PIK3CA_log_counts_by_age.jpg")
boxplot(PIK3CA_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for PIK3CA by Age Category",
        xlab = "Age Category", ylab = "PIK3CA Log Counts")
dev.off()

jpeg("TTN_log_counts_by_age.jpg")
boxplot(TTN_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for TTN by Age Category",
        xlab = "Age Category", ylab = "TTN Log Counts")
dev.off()

jpeg("GATA3_log_counts_by_age.jpg")
boxplot(GATA3_counts_log~age_category, data = patient_data, main = "Boxplot of HTSeq - Counts for GATA3 by Age Category",
        xlab = "Age Category", ylab = "GATA3 Log Counts")
dev.off()

#table of mean HTseq counts
mean_counts <- c(mean(patient_data$PIK3CA_counts_log[patient_data$age_category == "Young"]),
                 mean(patient_data$PIK3CA_counts_log[patient_data$age_category == "Old"]),
                 mean(patient_data$TP53_counts_log[patient_data$age_category == "Young"]),
                 mean(patient_data$TP53_counts_log[patient_data$age_category == "Old"]),
                 mean(patient_data$TTN_counts_log[patient_data$age_category == "Young"]),
                 mean(patient_data$TTN_counts_log[patient_data$age_category == "Old"]),
                 mean(patient_data$GATA3_counts_log[patient_data$age_category == "Young"]),
                 mean(patient_data$GATA3_counts_log[patient_data$age_category == "Old"])
)

mean_counts_matrix <- matrix(mean_counts, nrow = 2, ncol = 4, byrow = F)
colnames(mean_counts_matrix) <- c("PIK3CA", 'TP53', 'TTN', 'GATA3')
rownames(mean_counts_matrix) <- c("Young", "Old")
mean_counts_matrix


#t-tests to test differnces in gene expression between young and old patients
t.test(patient_data$PIK3CA_counts_log[patient_data$age_category == "Young"],
       patient_data$PIK3CA_counts_log[patient_data$age_category == "Old"])

t.test(patient_data$TP53_counts_log[patient_data$age_category == "Young"],
       patient_data$TP53_counts_log[patient_data$age_category == "Old"])

t.test(patient_data$TTN_counts_log[patient_data$age_category == "Young"],
       patient_data$TTN_counts_log[patient_data$age_category == "Old"])

t.test(patient_data$GATA3_counts_log[patient_data$age_category == "Young"],
       patient_data$GATA3_counts_log[patient_data$age_category == "Old"])