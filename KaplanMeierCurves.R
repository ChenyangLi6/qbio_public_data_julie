#Kaplan Meier plots for the Genes of Interest

#Mutated Genes in Young: GATA3
#Mutated Genes in Old: TTN

#DESeq2 Young: TRH
#DESeq2 Old: CARTPT


# Loading packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)

# Barcodes
# barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
#                     "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
#                    "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
# barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
#                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

# Accessing RNAseq data 
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
# barcode = barcodes_rnaseq
# GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

# Accessing Clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type = "xml")
# barcode=barcodes_clinic
# GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

# Getting counts for specific genes, replacing column names with shorter barcodes (to match clinical)
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient

# Matching clinical data sample order to RNAseq sample order
row_order <- match(colnames(htseq_counts), clinic$bcr_patient_barcode) 
# match function takes two parameters: first is the order you wish to match something to, and the second is the dataframe
# you wish to alter the order of.
clinic_ordered  <- clinic[row_order, ]

# Get rid of nonmatching samples in clinical and htseq
matching <- which(clinic_ordered$bcr_patient_barcode %in% colnames(htseq_counts))
# which function basically only takes the values of clinic_ordered$bcr_patient_barcode that are also found in the colnames of htseq
clinic_matched <- clinic_ordered[matching, ]

# Adding age to clinical data
age_clinical = clinic_matched$age_at_initial_pathologic_diagnosis
clinic_matched$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))
# Getting rid of middle-aged patients from clinical
clinic_matched <- subset(clinic_matched, clinic_matched$age_category=="Old" | clinic_matched$age_category=="Young")

# Accessing counts data for TTN, categorize expression, and add to clinical data
TTN_mask <- rowData(sum_exp)$external_gene_name == "TTN"
TTN_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ TTN_mask ]
TTN_counts <- htseq_counts[TTN_ENSG_name, clinic_matched$bcr_patient_barcode]
TTN_quartiles <- quantile(TTN_counts) # Categorizing the expression level based on quartile analysis
TTN_expression_level <- ifelse(TTN_counts > TTN_quartiles[4], "High", ifelse(TTN_counts < TTN_quartiles[2], "Low", "Mid"))
clinic_matched$TTN_expression = TTN_expression_level

# Accessing counts data for GATA3, categorize expression, and add to clinical data
GATA3_mask <- rowData(sum_exp)$external_gene_name == "GATA3"
GATA3_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ GATA3_mask ]
GATA3_counts <- htseq_counts[GATA3_ENSG_name, clinic_matched$bcr_patient_barcode]
GATA3_quartiles <- quantile(GATA3_counts) # Categorizing the expression level based on quartile analysis
GATA3_expression_level <- ifelse(GATA3_counts > GATA3_quartiles[4], "High", ifelse(GATA3_counts < GATA3_quartiles[2], "Low", "Mid"))
clinic_matched$GATA3_expression = GATA3_expression_level

# Accessing counts data for ARHGAP36, categorize expression, and add to clinical data
TRH_mask <- rowData(sum_exp)$external_gene_name == "TRH"
TRH_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ TRH_mask ]
TRH_counts <- htseq_counts[TRH_ENSG_name, clinic_matched$bcr_patient_barcode]
TRH_quartiles <- quantile(TRH_counts) # Categorizing the expression level based on quartile analysis
TRH_expression_level <- ifelse(TRH_counts > TRH_quartiles[4], "High", ifelse(TRH_counts < TRH_quartiles[2], "Low", "Mid"))
clinic_matched$TRH_expression = TRH_expression_level

# Accessing counts data for CARTPT, categorize expression, and add to clinical data
CARTPT_mask <- rowData(sum_exp)$external_gene_name == "CARTPT"
CARTPT_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ CARTPT_mask ]
CARTPT_counts <- htseq_counts[CARTPT_ENSG_name, clinic_matched$bcr_patient_barcode]
CARTPT_quartiles <- quantile(CARTPT_counts) # Categorizing the expression level based on quartile analysis
CARTPT_expression_level <- ifelse(CARTPT_counts > CARTPT_quartiles[4], "High", ifelse(CARTPT_counts < CARTPT_quartiles[2], "Low", "Mid"))
clinic_matched$CARTPT_expression = CARTPT_expression_level

# Taking only the patients with high expression of the genes of interest
clinic_TTN_high <- subset(clinic_matched, TTN_expression_level=="High")
clinic_GATA3_high <- subset(clinic_matched, GATA3_expression_level=="High")
clinic_TRH_high <- subset(clinic_matched, TRH_expression_level=="High")
clinic_CARTPT_high <- subset(clinic_matched, CARTPT_expression_level=="High")

# Plotting the Kaplan Meier survival graphs
TCGAanalyze_survival( clinic_TTN_high, "age_category", main="Kaplan-Meier Survival Curves for Patients with High TTN Expression", filename="TTN_survival.pdf")
TCGAanalyze_survival( clinic_GATA3_high, "age_category", main="Kaplan-Meier Survival Curves for Patients with High GATA3 Expression", filename="GATA3_survival.pdf")
TCGAanalyze_survival( clinic_TRH_high, "age_category", main="Kaplan-Meier Survival Curves for Patients with High TRH Expression", filename="TRH_survival.pdf")
TCGAanalyze_survival( clinic_CARTPT_high, "age_category", main="Kaplan-Meier Survival Curves for Patients with High CARTPT Expression", filename="CARTPT_survival.pdf")