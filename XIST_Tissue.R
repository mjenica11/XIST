#!/usr/bin/env Rscript

# Perform linear regression on the gene count data (version 7).
# Fit linear models: 
# Reponse variable: Median X chromosome expression and various sets of X-linked genes 
# Predictor variable: XIST
# Inlude expressed genes on the X chromosome only 
setwd("~/XIST/")

# Constants
COUNTS <- "~/XIST/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct" # TPM normalized
METRICS <- "~/XIST/Files/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"
PHENOTYPES <- "~/XIST/Files/GTEX_v7_Annotations_SubjectPhenotypesDS.txt"
GENCODE <- "gencode.v19.genes.v7.patched_contigs.gff3"
GENE_LST <- "~/XIST/Files/X_Genes_Status.json"

# Results
LM_FEM <- "~/XIST/Tissue/Xpressed/Median/Female_Tissue_Correlations.csv"
LM_MALE <- "~/XIST/Tissue/Xpressed/Median/Male_Tissue_Correlations.csv"

# Session data 
DATA <- "~/XIST/Tissue/Xpressed/Median/Xpressed_012720.RData"

# Load libraries
library(readr) 
library(refGenome) # Parse .gff file
library(plyr) # Load plyr before dplyer
library(dplyr) 
library(data.table) 
library(stringr)
library(broom)
library(rjson)
library(robustbase)

# Read in files
Metrics <- read_tsv(METRICS) # Contains tissue sample info
Phenotypes <- read_tsv(PHENOTYPES) # Contains sex info
Gene_Cts <- fread(COUNTS) # df of counts
Gene_Lst <- fromJSON(file=GENE_LST)

# Create ensemblGenome object for storing Ensembl genomic annotation data
ENS <- ensemblGenome()

# Read in .gff annotation file as ensemblGenome object
read.gtf(ENS, GENCODE)

# ______________________________________________________________________________________________________________________
# Get list of genes on X chromosome
# ______________________________________________________________________________________________________________________
# Get annotations for X chromosome
X_Annot <- extractSeqids(ENS, 'X')

# List genes on X chromosome
X_Genes = X_Annot@ev$gtf[, c("gene_id","gene_name")]

# Remove numbers after decimal in gene ID to get the gene IDs, not the exon IDs
X_Genes$gene_id <- sub("\\.\\d+$", "", X_Genes$gene_id)

# Remove duplicate genes
# Only difference bw transcript and gene IDs is the decimal after the gene ID
X_Genes <- X_Genes[!duplicated(X_Genes$gene_id), ]

# Drop XIST from list of X chromosome genes
# Don't want to include XIST in median X chromsome count values
X_Genes <- X_Genes[!X_Genes$gene_name == 'XIST',]
# ______________________________________________________________________________________________________________________
# Organize samples into list of dfs by tissue type
# ______________________________________________________________________________________________________________________
# Drop columns in Metrics that aren't needed 
Metrics <- Metrics %>% select(Sample, Note)

# Rename column
colnames(Metrics)[2] <- "Tissue"

# Get list of individual GTEx IDs
Individual_IDs <- unique(str_extract(Metrics$Sample, "GTEX-[0-9A-Z]+"))

# Remove any missing values
which(is.na(Individual_IDs)) # 738; last item in list
Individual_IDs <- Individual_IDs[!is.na(Individual_IDs)]

# For each tissue, make a df of samples from that tissue and store in list.
Types <- unique(Metrics$Tissue)

Tissue_Lst <- list()
for (i in Types){
  Tissue_Lst[[i]] <- Metrics[Metrics$Tissue == i,] 
}

# Add column containing sample ID to each df in list
for (i in seq_along(Tissue_Lst)){
  Tissue_Lst[[i]]$ID <- str_extract(Tissue_Lst[[i]]$Sample, "GTEX-[0-9A-Z]+")
}

# Check for missing values
sapply(Tissue_Lst, function(x) sum(is.na(x))) # Cells - Leukemia cell line (CML) 

# Drop cell lines
Cells <- c("Cells - Leukemia cell line (CML)", "Cells - EBV-transformed lymphocytes", "Cells - Transformed fibroblasts")
Tissue_Lst <- Tissue_Lst[!names(Tissue_Lst) %in% Cells]

# Get list of female IDs
Fem.IDs <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Get list of sample sex for each tissue
Sex_Lst <- lapply(Tissue_Lst, function(x) {
  x <- with(x['ID'], ifelse(ID %in% Fem.IDs, "Female", "Male"))
})

# Add column containing sample sex to each df in list 
Tissue_Lst <- mapply(cbind, Tissue_Lst, "Sex"=Sex_Lst, SIMPLIFY=F)

# For each individual, make a data frame of samples that comes from the same 
# person and store in list.
Ind_Tissues <- list()
for (i in Individual_IDs){
  Ind_Tissues[[i]] <- Metrics[Metrics$Sample %like% i, ]
}

# Get list of sample replicates
# i.e. people who have multiple samples for the same tissue type
Dup_Ind_Tissues <- lapply(Ind_Tissues, function(x) {
  x[duplicated(x[,2]), ]
})
names(Dup_Ind_Tissues) <- names(Ind_Tissues)

Sample_Replicates <- lapply(Dup_Ind_Tissues, function(x) {
  if (as.character(x['Sample']) != "character(0)") {
    as.character(x['Sample'])
  } else {}
})
Sample_Replicates <- unlist(Sample_Replicates)

# Remove sample replicates
Tissue_Lst <- lapply(Tissue_Lst, function(x){
  x[!(x$ID %in% Sample_Replicates), ]
})

# ______________________________________________________________________________________________________________________
# Organize sample counts by sex and tissue type; drop replicates
# ______________________________________________________________________________________________________________________
# For each individual, make a data frame of gene counts from samples that come 
# from the same person and store in list.
Gene_Cts <- data.frame(Gene_Cts, stringsAsFactors = F) # was both data.table and data frame
colnames(Gene_Cts) <- str_replace_all(colnames(Gene_Cts), pattern = "\\.", replacement = "-")

# Rename columns
names(Gene_Cts)[1:2] <- c("gene_id", "gene_name")

# Remove decimals in gene IDs in gene counts df
Gene_Cts$gene_id <- sub("\\.\\d+$", "", Gene_Cts$gene_id)

# Split tissue lists by sex
Split_Func <- function(x, y){
  x <- x[x[['Sex']] == y,]
  return(x)
}
f.Tissue_Lst <- Map(Split_Func, x=Tissue_Lst, y='Female')
m.Tissue_Lst <- Map(Split_Func, x=Tissue_Lst, y='Male')

# Drop tissues females/males do not have, respectively
f.remove <- c("Prostate", "Testis", "Cells - Leukemia cell line (CML)")
m.remove <- c("Ovary", "Uterus", "Vagina", "Fallopian Tube", "Cervix - Ectocervix",
              "Cervix - Endocervix", "Cells - Leukemia cell line (CML)")

f.Tissue_Lst <- f.Tissue_Lst[!names(f.Tissue_Lst) %in% f.remove]
m.Tissue_Lst <- m.Tissue_Lst[!names(m.Tissue_Lst) %in% m.remove]

# Sort count data into list of dfs by tissue 
Sort_Func <- function(x){
  tmp <- which(colnames(Gene_Cts) %in% x$Sample)
  Gene_Cts[, c(1,2,tmp)]
}
f.Tissue_Counts <- lapply(f.Tissue_Lst, Sort_Func)
m.Tissue_Counts <- lapply(m.Tissue_Lst, Sort_Func)

# Remove sample replicates
Drop_Replicates <- function(x){
  x <- x[, !(names(x) %in% Sample_Replicates)]
  return(x)
}
f.Tissue_Counts <- lapply(f.Tissue_Counts, Drop_Replicates)
m.Tissue_Counts <- lapply(m.Tissue_Counts, Drop_Replicates)

# ______________________________________________________________________________________________________________________
# Get subset of just X chm genes from each data frame in list
# ______________________________________________________________________________________________________________________
# Subset X chr genes; Will not include XIST
XSubset_Func <- function(x){
  res <- filter(x, gene_id %in% X_Genes$gene_id)
  return(res)
}
f.Tissue_XCounts <- lapply(f.Tissue_Counts, XSubset_Func)
m.Tissue_XCounts <- lapply(m.Tissue_Counts, XSubset_Func)

# Remove X chr genes that are unexpressed
# Function to filter list of count dfs by tpm value
Filter_Cts <- function(DF, TPM){
  DF <- DF[apply(DF[,3:ncol(DF)], MARGIN=1, function(x) all(x>TPM)),]
  return(DF)
}
# Keep genes with TPM > 0; counts are seperated by sex
f.Tissue_XCounts <- Map(Filter_Cts, DF=f.Tissue_XCounts, TPM=0)
m.Tissue_XCounts <- Map(Filter_Cts, DF=m.Tissue_XCounts, TPM=0)

# Check for missing values
any(sapply(f.Tissue_XCounts, function(x) sum(is.na(x)))) # FALSE
any(sapply(m.Tissue_XCounts, function(x) sum(is.na(x)))) # FALSE

# Get the median values of expressed X chm genes
Median_Val <- function(x){
  res <- apply(x[,-c(1:2)], 2, median) 
  return(res)
}
f.MedX_Tissue_Counts <- lapply(f.Tissue_XCounts, Median_Val)
m.MedX_Tissue_Counts <- lapply(m.Tissue_XCounts, Median_Val)

# Get XIST values from each data frame
Get_XIST <- function(x){
  res <- as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
  return(res)
}
f.XIST_Tissue_Counts <- lapply(f.Tissue_Counts, Get_XIST)
m.XIST_Tissue_Counts <- lapply(m.Tissue_Counts, Get_XIST)


### This will break the linear model: many male tissues have no samples
### where XIST > 1; thus you will be plotting x=0 for every point.
# Replace XIST counts with TPM < 1 with 0
# Filter_Vec <- function(v){
#   res <- ifelse(v<1, 0.001, v)
#   return(res)
# }
# f.XIST_Tissue_Counts <- lapply(f.XIST_Tissue_Counts, Filter_Vec)
# m.XIST_Tissue_Counts <- lapply(m.XIST_Tissue_Counts, Filter_Vec)

# ______________________________________________________________________________________________________________________
# Sanity check
# ______________________________________________________________________________________________________________________
# Check that all the samples are present in both lists of dfs
length(f.XIST_Tissue_Counts) == length(f.MedX_Tissue_Counts) # TRUE
length(m.XIST_Tissue_Counts) == length(m.MedX_Tissue_Counts) # TRUE

# How many tissue types are there per sex?
length(f.XIST_Tissue_Counts) # 51
length(m.XIST_Tissue_Counts) # 47

# Are all the same individuals listed in both lists?
all(names(f.XIST_Tissue_Counts) == names(f.MedX_Tissue_Counts)) # TRUE
all(names(m.XIST_Tissue_Counts) == names(m.MedX_Tissue_Counts)) # TRUE

# Are they listed in the same order?
identical(names(f.XIST_Tissue_Counts), names(f.MedX_Tissue_Counts)) # TRUE
identical(names(m.XIST_Tissue_Counts), names(m.MedX_Tissue_Counts)) # TRUE

# ______________________________________________________________________________________________________________________
# Table 1; Column 1:2
# Correlate expression from XIST with X chromosome expression within a person across all of their tissues.
# ______________________________________________________________________________________________________________________
# Function to combine vectors into df
Combine_Vectors <- function(a, b){
  data.frame(a, b)
}

# Combine list of median value of silenced genes and XIST expression
Combine_Lsts <- function(x, y, z){
  x <- list()
  for (i in 1:length(z)){
    x[[i]] <- Combine_Vectors(y[[i]], z[[i]])
  }
  names(x) <- names(y)
  return(x)
}
f.MedX_Vs_XIST <- Combine_Lsts(x='f.MedX_Vs_XIST', 
                                y=f.MedX_Tissue_Counts, 
                                z=f.XIST_Tissue_Counts)
m.MedX_Vs_XIST <- Combine_Lsts(x='m.MedX_Vs_XIST', 
                                y=m.MedX_Tissue_Counts, 
                                z=m.XIST_Tissue_Counts)

# Rename col names in each df
Rename_Col <- function(x, a, b){
  x <- setNames(x, c(a, b))
  return(x)
}
f.MedX_Vs_XIST <- Map(Rename_Col,
                       x=f.MedX_Vs_XIST, 
                       a='MedX', 
                       b='XIST')
m.MedX_Vs_XIST <- Map(Rename_Col, 
                       x=m.MedX_Vs_XIST, 
                       a='MedX', 
                       b='XIST')

# Apply lm to each df in list
Linear_Model.1 <- function(x) {
  z <- lm(MedX ~ XIST, data = x)
  return(z)
}
lm_f.MedX_XIST <- lapply(f.MedX_Vs_XIST, Linear_Model.1)
lm_m.MedX_XIST <- lapply(m.MedX_Vs_XIST, Linear_Model.1)

# Function to extract r squared values
Regression_Res <- function(lm){ # expecting object of class 'lm'
  sum. <- summary(lm)
  r.2 <- sum.$r.squared
  p. <- summary(lm)$fstatistic
  p.val <- pf(p.[1], p.[2], p.[3], lower.tail=FALSE, log.p=FALSE) #pf: F distribution; arg 2 and 3 are the deg.freedom 
  attributes(p.val) <- NULL
  res <- data.frame(p_val=p.val, r_2=r.2)
  return(res)
}

# Apply function to list of dfs
Res_f.MedX_XIST <- lapply(lm_f.MedX_XIST, Regression_Res)
Res_m.MedX_XIST <- lapply(lm_m.MedX_XIST, Regression_Res)

# Make table summarizing results
f.Regression <- as.data.frame(do.call(rbind, Res_f.MedX_XIST))
m.Regression <- as.data.frame(do.call(rbind, Res_m.MedX_XIST))

# Label columns
colnames(f.Regression) <- c("pval_MedX", "R2_MedX")
colnames(m.Regression) <- c("pval_MedX","R2_MedX")

# ______________________________________________________________________________________________________________________
# Table 1; Columns 3:12
# Correlation of all genes reported as silenced with XIST
# ______________________________________________________________________________________________________________________
# Categories:
# "Silenced_In_Both", "Silenced_In_Tukiainen", "Silenced_In_Balaton", "Silenced_In_At_Least_One", 
# "Immune_Genes_Silenced_In_At_Least_One"

# Subset list of genes silenced in both studies from each df of x counts in list by sex
Filter_Func <- function(x, y){
  x <- filter(x, gene_name %in% Gene_Lst[[y]])
  return(x)
}

f.Both_Silenced <- Map(Filter_Func, 
                       x=f.Tissue_XCounts, 
                       y='Silenced_In_Both')
f.Tuk_Silenced <- Map(Filter_Func, 
                      x=f.Tissue_XCounts, 
                      y='Silenced_In_Tukiainen')
f.Bal_Silenced <- Map(Filter_Func, 
                      x=f.Tissue_XCounts, 
                      y='Silenced_In_Balaton')
f.One_Silenced <- Map(Filter_Func, 
                      x=f.Tissue_XCounts,
                      y='Silenced_In_At_Least_One')
f.Immune_Silenced <- Map(Filter_Func, 
                         x=f.Tissue_XCounts, 
                         y='Immune_Genes_Silenced_In_At_Least_One')

m.Both_Silenced <- Map(Filter_Func, 
                       x=m.Tissue_XCounts, 
                       y='Silenced_In_Both')
m.Tuk_Silenced <- Map(Filter_Func, 
                      x=m.Tissue_XCounts, 
                      y='Silenced_In_Tukiainan')
m.Bal_Silenced <- Map(Filter_Func, 
                      x=m.Tissue_XCounts, 
                      y='Silenced_In_Balaton')
m.One_Silenced <- Map(Filter_Func, 
                      x=m.Tissue_XCounts, 
                      y='Silenced_In_At_Least_One')
m.Immune_Silenced <- Map(Filter_Func, 
                         x=m.Tissue_XCounts, 
                         y='Immune_Genes_Silenced_In_At_Least_One')

# Get the median values of putatively silenced X chm genes
f.median_Silenced <- lapply(f.Both_Silenced, Median_Val)
f.Tuk_median_Silenced <- lapply(f.Both_Silenced, Median_Val)
f.Bal_median_Silenced <- lapply(f.Bal_Silenced, Median_Val)
f.One_median_Silenced <- lapply(f.One_Silenced, Median_Val)
f.Immune_median_Silenced <- lapply(f.Immune_Silenced, Median_Val)

m.median_Silenced <- lapply(m.Both_Silenced, Median_Val)
m.Tuk_median_Silenced <- lapply(m.Both_Silenced, Median_Val)
m.Bal_median_Silenced <- lapply(m.Bal_Silenced, Median_Val)
m.One_median_Silenced <- lapply(m.One_Silenced, Median_Val)
m.Immune_median_Silenced <- lapply(m.Immune_Silenced, Median_Val)

# Combine list of median value of silenced genes and XIST expression
f.Silenced_median_Vs_XIST <- Combine_Lsts(x='f.Silenced_median_Vs_XIST', 
                                        y=f.median_Silenced, 
                                        z=f.XIST_Tissue_Counts)
f.Tuk_Silenced_median_Vs_XIST <- Combine_Lsts(x='f.Tuk_Silenced_median_Vs_XIST', 
                                            y=f.Tuk_median_Silenced, 
                                            z=f.XIST_Tissue_Counts)
f.Bal_Silenced_median_Vs_XIST <- Combine_Lsts(x='f.Bal_Silenced_median_Vs_XIST', 
                                            y=f.Bal_median_Silenced, 
                                            z=f.XIST_Tissue_Counts)
f.One_Silenced_median_Vs_XIST <- Combine_Lsts(x='f.One_Silenced_median_Vs_XIST', 
                                            y=f.One_median_Silenced, 
                                            z=f.XIST_Tissue_Counts)
f.Immune_Silenced_median_Vs_XIST <- Combine_Lsts(x='f.Immune_Silenced_median_Vs_XIST', 
                                               y=f.Immune_median_Silenced, 
                                               z=f.XIST_Tissue_Counts)

m.Silenced_median_Vs_XIST <- Combine_Lsts(x='m.Silenced_median_Vs_XIST', 
                                        y=m.median_Silenced, 
                                        z=m.XIST_Tissue_Counts)
m.Tuk_Silenced_median_Vs_XIST <- Combine_Lsts(x='m.Tuk_Silenced_median_Vs_XIST', 
                                            y=m.Tuk_median_Silenced, 
                                            z=m.XIST_Tissue_Counts)
m.Bal_Silenced_median_Vs_XIST <- Combine_Lsts(x='m.Bal_Silenced_median_Vs_XIST', 
                                            y=m.Bal_median_Silenced, 
                                            z=m.XIST_Tissue_Counts)
m.One_Silenced_median_Vs_XIST <- Combine_Lsts(x='m.One_Silenced_median_Vs_XIST', 
                                            y=m.One_median_Silenced, 
                                            z=m.XIST_Tissue_Counts)
m.Immune_Silenced_median_Vs_XIST <- Combine_Lsts(x='m.Immune_Silenced_median_Vs_XIST', 
                                               y=m.Immune_median_Silenced, 
                                               z=m.XIST_Tissue_Counts)

# Rename columns in each df
f.Silenced_median_Vs_XIST <- Map(Rename_Col, 
                               x=f.Silenced_median_Vs_XIST, 
                               a='median_Silenced', 
                               b='XIST')
f.Tuk_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                   x=f.Tuk_Silenced_median_Vs_XIST, 
                                   a='median_Silenced', 
                                   b='XIST')
f.Bal_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                   x=f.Bal_Silenced_median_Vs_XIST, 
                                   a='median_Silenced', 
                                   b='XIST')
f.One_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                   x=f.One_Silenced_median_Vs_XIST, 
                                   a='median_Silenced', 
                                   b='XIST')
f.Immune_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                      x=f.Immune_Silenced_median_Vs_XIST, 
                                      a='median_Silenced', 
                                      b='XIST')

m.Silenced_median_Vs_XIST <- Map(Rename_Col, 
                               x=m.Silenced_median_Vs_XIST, 
                               a='median_Silenced', 
                               b='XIST')
m.Tuk_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                   x=m.Tuk_Silenced_median_Vs_XIST, 
                                   a='median_Silenced', 
                                   b='XIST')
m.Bal_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                   x=m.Bal_Silenced_median_Vs_XIST, 
                                   a='median_Silenced', 
                                   b='XIST')
m.One_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                   x=m.One_Silenced_median_Vs_XIST, 
                                   a='median_Silenced', 
                                   b='XIST')
m.Immune_Silenced_median_Vs_XIST <- Map(Rename_Col, 
                                      x=m.Immune_Silenced_median_Vs_XIST, 
                                      a='median_Silenced', 
                                      b='XIST')

# Apply lm to each df in list
Linear_Model.2 <- function(x) {
  z <- lm(median_Silenced ~ XIST, data = x)
  return(z)
}

f.lm_Silenced_XIST <- lapply(f.Silenced_median_Vs_XIST, Linear_Model.2)
f.lm_Tuk_Silenced_XIST <- lapply(f.Tuk_Silenced_median_Vs_XIST, Linear_Model.2)
f.lm_Bal_Silenced_XIST <- lapply(f.Bal_Silenced_median_Vs_XIST, Linear_Model.2)
f.lm_One_Silenced_XIST <- lapply(f.One_Silenced_median_Vs_XIST, Linear_Model.2)
f.lm_Immune_Silenced_XIST <- lapply(f.Immune_Silenced_median_Vs_XIST, Linear_Model.2)

m.lm_Silenced_XIST <- lapply(m.Silenced_median_Vs_XIST, Linear_Model.2)
m.lm_Tuk_Silenced_XIST <- lapply(m.Tuk_Silenced_median_Vs_XIST, Linear_Model.2)
m.lm_Bal_Silenced_XIST <- lapply(m.Bal_Silenced_median_Vs_XIST, Linear_Model.2)
m.lm_One_Silenced_XIST <- lapply(m.One_Silenced_median_Vs_XIST, Linear_Model.2)
m.lm_Immune_Silenced_XIST <- lapply(m.Immune_Silenced_median_Vs_XIST, Linear_Model.2)

# Apply function to list of dfs
f.Res_Silenced_XIST <- lapply(f.lm_Silenced_XIST, Regression_Res)
f.Tuk_Res_Silenced_XIST <- lapply(f.lm_Tuk_Silenced_XIST, Regression_Res)
f.Bal_Res_Silenced_XIST <- lapply(f.lm_Bal_Silenced_XIST, Regression_Res)
f.One_Res_Silenced_XIST <- lapply(f.lm_One_Silenced_XIST, Regression_Res)
f.Immune_Res_Silenced_XIST <- lapply(f.lm_Immune_Silenced_XIST, Regression_Res)

m.Res_Silenced_XIST <- lapply(m.lm_Silenced_XIST, Regression_Res)
m.Tuk_Res_Silenced_XIST <- lapply(m.lm_Tuk_Silenced_XIST, Regression_Res)
m.Bal_Res_Silenced_XIST <- lapply(m.lm_Bal_Silenced_XIST, Regression_Res)
m.One_Res_Silenced_XIST <- lapply(m.lm_One_Silenced_XIST, Regression_Res)
m.Immune_Res_Silenced_XIST <- lapply(m.lm_Immune_Silenced_XIST, Regression_Res)

# Add vector as column to df
Return_R2 <- function(x){
  res <- unlist(x$r_2)
  return(res)
}

Return_pval <- function(x){
  res <- unlist(x$p_val)
  return(res)
}

f.Regression$pval_Silenced_median <- lapply(f.Res_Silenced_XIST, Return_pval)
f.Regression$R2_Silenced_median <- lapply(f.Res_Silenced_XIST, Return_R2)

f.Regression$pval_Tuk_Silenced_median <- lapply(f.Tuk_Res_Silenced_XIST, Return_pval)
f.Regression$R2_Tuk_Silenced_median <- lapply(f.Tuk_Res_Silenced_XIST, Return_R2)

f.Regression$pval_Bal_Silenced_median <- lapply(f.Bal_Res_Silenced_XIST, Return_pval)
f.Regression$R2_Bal_Silenced_median <- lapply(f.Bal_Res_Silenced_XIST, Return_R2)

f.Regression$pval_One_Silenced_median <- lapply(f.One_Res_Silenced_XIST, Return_pval)
f.Regression$R2_One_Silenced_median <- lapply(f.One_Res_Silenced_XIST, Return_R2)

f.Regression$pval_Immune_Silenced_median <- lapply(f.Immune_Res_Silenced_XIST, Return_pval)
f.Regression$R2_Immune_Silenced_median <- lapply(f.Immune_Res_Silenced_XIST, Return_R2)

m.Regression$pval_Silenced_median <- lapply(m.Res_Silenced_XIST, Return_pval)
m.Regression$R2_Silenced_median <- lapply(m.Res_Silenced_XIST, Return_R2)

m.Regression$pval_Tuk_Silenced_median <- lapply(m.Tuk_Res_Silenced_XIST, Return_pval)
m.Regression$R2_Tuk_Silenced_median <- lapply(m.Tuk_Res_Silenced_XIST, Return_R2)

m.Regression$pval_Bal_Silenced_median <- lapply(m.Bal_Res_Silenced_XIST, Return_pval)
m.Regression$R2_Bal_Silenced_median <- lapply(m.Bal_Res_Silenced_XIST, Return_R2)

m.Regression$pval_One_Silenced_median <- lapply(m.One_Res_Silenced_XIST, Return_pval)
m.Regression$R2_One_Silenced_median <- lapply(m.One_Res_Silenced_XIST, Return_R2)

m.Regression$pval_Immune_Silenced_median <- lapply(m.Immune_Res_Silenced_XIST, Return_pval)
m.Regression$R2_Immune_Silenced_median <- lapply(m.Immune_Res_Silenced_XIST, Return_R2)

# ______________________________________________________________________________________________________________________
# Table 1; Columns 13:20
# Correlation of all genes reported as variably silenced with XIST
# ______________________________________________________________________________________________________________________
# Categories:
# "Variable_In_At_Least_One", "Variable_In_Tukiainen", "Variable_In_Balaton", "Immune_Gene_Variable_In_At_Least_One"

# Subset list of genes variably silenced from each df of x counts in list 
f.One_Variable <- Map(Filter_Func, 
                      x=f.Tissue_XCounts, 
                      y='Variable_In_At_Least_One')
f.Tuk_Variable <- Map(Filter_Func, 
                      x=f.Tissue_XCounts, 
                      y='Variable_In_Tukiainen')
f.Bal_Variable <- Map(Filter_Func, 
                      x=f.Tissue_XCounts, 
                      y='Variable_In_Balaton')
f.Immune_Variable <- Map(Filter_Func, 
                         x=f.Tissue_XCounts, 
                         y='Immune_Gene_Variable_In_At_Least_One')

m.One_Variable <- Map(Filter_Func, 
                      x=m.Tissue_XCounts, 
                      y='Variable_In_At_Least_One')
m.Tuk_Variable <- Map(Filter_Func, 
                      x=m.Tissue_XCounts, 
                      y='Variable_In_Tukiainen')
m.Bal_Variable <- Map(Filter_Func, 
                      x=m.Tissue_XCounts, 
                      y='Variable_In_Balaton')
m.Immune_Variable <- Map(Filter_Func, 
                         x=m.Tissue_XCounts, 
                         y='Immune_Gene_Variable_In_At_Least_One')

# Get the median values of putatively variably silenced X chm genes
f.One_median_Variable <- lapply(f.One_Variable, Median_Val)
f.Tuk_median_Variable <- lapply(f.Tuk_Variable, Median_Val)
f.Bal_median_Variable <- lapply(f.Bal_Variable, Median_Val)
f.Immune_median_Variable <- lapply(f.Immune_Variable, Median_Val)

m.One_median_Variable <- lapply(m.One_Variable, Median_Val)
m.Tuk_median_Variable <- lapply(m.Tuk_Variable, Median_Val)
m.Bal_median_Variable <- lapply(m.Bal_Variable, Median_Val)
m.Immune_median_Variable <- lapply(m.Immune_Variable, Median_Val)

# Combine list of median value of variably silenced genes and XIST expression
f.One_median_Variable_Vs_XIST <- Combine_Lsts(x='f.One_median_Variable_Vs_XIST', 
                                            y=f.One_median_Variable, 
                                            z=f.XIST_Tissue_Counts)
f.Tuk_median_Variable_Vs_XIST <- Combine_Lsts(x='f.Tuk_median_Variable_Vs_XIST', 
                                            y=f.Tuk_median_Variable, 
                                            z=f.XIST_Tissue_Counts)
f.Bal_median_Variable_Vs_XIST <- Combine_Lsts(x='f.Bal_median_Variable_Vs_XIST', 
                                            y=f.Bal_median_Variable, 
                                            z=f.XIST_Tissue_Counts)
f.Immune_median_Variable_Vs_XIST <- Combine_Lsts(x='f.Immune_median_Variable_Vs_XIST', 
                                               y=f.Immune_median_Variable, 
                                               z=f.XIST_Tissue_Counts)

m.One_median_Variable_Vs_XIST <- Combine_Lsts(x='m.One_median_Variable_Vs_XIST', 
                                            y=m.One_median_Variable, 
                                            z=m.XIST_Tissue_Counts)
m.Tuk_median_Variable_Vs_XIST <- Combine_Lsts(x='m.Tuk_median_Variable_Vs_XIST', 
                                            y=m.Tuk_median_Variable, 
                                            z=m.XIST_Tissue_Counts)
m.Bal_median_Variable_Vs_XIST <- Combine_Lsts(x='m.Bal_median_Variable_Vs_XIST', 
                                            y=m.Bal_median_Variable, 
                                            z=m.XIST_Tissue_Counts)
m.Immune_median_Variable_Vs_XIST <- Combine_Lsts(x='m.Immune_median_Variable_Vs_XIST',
                                               y=m.Immune_median_Variable, 
                                               z=m.XIST_Tissue_Counts)

# Rename columns in each df
f.One_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                   x=f.One_median_Variable_Vs_XIST,
                                   a='median_Variable', 
                                   b='XIST')
f.Tuk_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                   x=f.Tuk_median_Variable_Vs_XIST, 
                                   a='median_Variable', 
                                   b='XIST')
f.Bal_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                   x=f.Bal_median_Variable_Vs_XIST, 
                                   a='median_Variable', 
                                   b='XIST')
f.Immune_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                      x=f.Immune_median_Variable_Vs_XIST, 
                                      a='median_Variable', 
                                      b='XIST') 

m.One_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                   x=m.One_median_Variable_Vs_XIST, 
                                   a='median_Variable', 
                                   b='XIST')
m.Tuk_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                   x=m.Tuk_median_Variable_Vs_XIST, 
                                   a='median_Variable', 
                                   b='XIST')
m.Bal_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                   x=m.Bal_median_Variable_Vs_XIST, 
                                   a='median_Variable', 
                                   b='XIST')
m.Immune_median_Variable_Vs_XIST <- Map(Rename_Col, 
                                      x=m.Immune_median_Variable_Vs_XIST, 
                                      a='median_Variable', 
                                      b='XIST')

# Apply lm to each df in list
Linear_Model.3 <- function(x) {
  z <- lm(median_Variable ~ XIST, data = x)
  return(z)
}

f.lm_One_Variable_XIST <- lapply(f.One_median_Variable_Vs_XIST, Linear_Model.3)
f.lm_Tuk_Variable_XIST <- lapply(f.Tuk_median_Variable_Vs_XIST, Linear_Model.3)  
f.lm_Bal_Variable_XIST <- lapply(f.Bal_median_Variable_Vs_XIST, Linear_Model.3)  
f.lm_Immune_Variable_XIST <- lapply(f.Immune_median_Variable_Vs_XIST, Linear_Model.3)  

m.lm_One_Variable_XIST <- lapply(m.One_median_Variable_Vs_XIST, Linear_Model.3)
m.lm_Tuk_Variable_XIST <- lapply(m.Tuk_median_Variable_Vs_XIST, Linear_Model.3)  
m.lm_Bal_Variable_XIST <- lapply(m.Bal_median_Variable_Vs_XIST, Linear_Model.3)  
m.lm_Immune_Variable_XIST <- lapply(m.Immune_median_Variable_Vs_XIST, Linear_Model.3) 

# Apply function to list of dfs
f.Res_One_Variable_XIST <- lapply(f.lm_One_Variable_XIST, Regression_Res)
f.Res_Tuk_Variable_XIST <- lapply(f.lm_Tuk_Variable_XIST, Regression_Res)
f.Res_Bal_Variable_XIST <- lapply(f.lm_Bal_Variable_XIST, Regression_Res)
f.Res_Immune_Variable_XIST <- lapply(f.lm_Immune_Variable_XIST, Regression_Res)

m.Res_One_Variable_XIST <- lapply(m.lm_One_Variable_XIST, Regression_Res)
m.Res_Tuk_Variable_XIST <- lapply(m.lm_Tuk_Variable_XIST, Regression_Res)
m.Res_Bal_Variable_XIST <- lapply(m.lm_Bal_Variable_XIST, Regression_Res)
m.Res_Immune_Variable_XIST <- lapply(m.lm_Immune_Variable_XIST, Regression_Res)

# Add vector as column to df
f.Regression$pval_One_Variable_median <- lapply(f.Res_Silenced_XIST, Return_pval)
f.Regression$R2_One_Variable_median <- lapply(f.Res_Silenced_XIST, Return_R2)

f.Regression$pval_Tuk_Variable_median <- lapply(f.Res_Tuk_Variable_XIST, Return_pval)
f.Regression$R2_Tuk_Variable_median <- lapply(f.Res_Tuk_Variable_XIST, Return_R2)

f.Regression$pval_Bal_Variable_median <- lapply(f.Res_Bal_Variable_XIST, Return_pval)
f.Regression$R2_Bal_Variable_median <- lapply(f.Res_Bal_Variable_XIST, Return_R2)

f.Regression$pval_Immune_Variable_median <- lapply(f.Res_Immune_Variable_XIST, Return_pval)
f.Regression$R2_Immune_Variable_median <- lapply(f.Res_Immune_Variable_XIST, Return_R2)

m.Regression$pval_One_Variable_median <- lapply(m.Res_Silenced_XIST, Return_pval)
m.Regression$R2_One_Variable_median <- lapply(m.Res_Silenced_XIST, Return_R2)

m.Regression$pval_Tuk_Variable_median <- lapply(m.Res_Tuk_Variable_XIST, Return_pval)
m.Regression$R2_Tuk_Variable_median <- lapply(m.Res_Tuk_Variable_XIST, Return_R2)

m.Regression$pval_Bal_Variable_median <- lapply(m.Res_Bal_Variable_XIST, Return_pval)
m.Regression$R2_Bal_Variable_median <- lapply(m.Res_Bal_Variable_XIST, Return_R2)

m.Regression$pval_Immune_Variable_median <- lapply(m.Res_Immune_Variable_XIST, Return_pval)
m.Regression$R2_Immune_Variable_median <- lapply(m.Res_Immune_Variable_XIST, Return_R2)
# ______________________________________________________________________________________________________________________
# Table 1; Columns 21:28
# Correlation of all genes reported as incompletely silenced with XIST
# ______________________________________________________________________________________________________________________
# Categories:
# "Incomplete_In_At_Least_One", Incomplete_In_Tukiainen", "Incomplete_In_Balaton", 
# "Immune_Genes_Incomplete_In_At_Least_One"

# Subset list of genes incompletly silenced from each df of x counts in list 
f.One_Incomplete <- Map(Filter_Func, 
                        x=f.Tissue_XCounts, 
                        y='Incomplete_In_At_Least_One')
f.Tuk_Incomplete <- Map(Filter_Func, 
                        x=f.Tissue_XCounts, 
                        y='Incomplete_In_Tukiainen')
f.Bal_Incomplete <- Map(Filter_Func, 
                        x=f.Tissue_XCounts, 
                        y='Incomplete_In_Balaton')
f.Immune_Incomplete <- Map(Filter_Func, 
                           x=f.Tissue_XCounts, 
                           y='Immune_Genes_Incomplete_In_At_Least_One')

m.One_Incomplete <- Map(Filter_Func, 
                        x=m.Tissue_XCounts, 
                        y='Incomplete_In_At_Least_One')
m.Tuk_Incomplete <- Map(Filter_Func, 
                        x=m.Tissue_XCounts, 
                        y='Incomplete_In_Tukiainen')
m.Bal_Incomplete <- Map(Filter_Func, 
                        x=m.Tissue_XCounts, 
                        y='Incomplete_In_Balaton')
m.Immune_Incomplete <- Map(Filter_Func, 
                           x=m.Tissue_XCounts, 
                           y='Immune_Genes_Incomplete_In_At_Least_One')

# Get the median values of putatively incompletely silenced X chm genes
f.One_median_Incomplete <- lapply(f.One_Incomplete, Median_Val)
f.Tuk_median_Incomplete <- lapply(f.Tuk_Incomplete, Median_Val)
f.Bal_median_Incomplete <- lapply(f.Bal_Incomplete, Median_Val)
f.Immune_median_Incomplete <- lapply(f.Immune_Incomplete, Median_Val)

m.One_median_Incomplete <- lapply(m.One_Incomplete, Median_Val)
m.Tuk_median_Incomplete <- lapply(m.Tuk_Incomplete, Median_Val)
m.Bal_median_Incomplete <- lapply(m.Bal_Incomplete, Median_Val)
m.Immune_median_Incomplete <- lapply(m.Immune_Incomplete, Median_Val)

# Combine list of median value of incompletely silenced genes and XIST expression
f.One_median_Incomplete_Vs_XIST <- Combine_Lsts(x='f.One_median_Incomplete_Vs_XIST', 
                                              y=f.One_median_Incomplete, 
                                              z=f.XIST_Tissue_Counts)
f.Tuk_median_Incomplete_Vs_XIST <- Combine_Lsts(x='f.Tuk_median_Incomplete_Vs_XIST', 
                                              y=f.Tuk_median_Incomplete, 
                                              z=f.XIST_Tissue_Counts)
f.Bal_median_Incomplete_Vs_XIST <- Combine_Lsts(x='f.Bal_median_Incomplete_Vs_XIST', 
                                              y=f.Bal_median_Incomplete, 
                                              z=f.XIST_Tissue_Counts)
f.Immune_median_Incomplete_Vs_XIST <- Combine_Lsts(x='f.Immune_median_Incomplete_Vs_XIST', 
                                                 y=f.Immune_median_Incomplete, 
                                                 z=f.XIST_Tissue_Counts)

m.One_median_Incomplete_Vs_XIST <- Combine_Lsts(x='m.One_median_Incomplete_Vs_XIST', 
                                              y=m.One_median_Incomplete, 
                                              z=m.XIST_Tissue_Counts)
m.Tuk_median_Incomplete_Vs_XIST <- Combine_Lsts(x='m.Tuk_median_Incomplete_Vs_XIST', 
                                              y=m.Tuk_median_Incomplete, 
                                              z=m.XIST_Tissue_Counts)
m.Bal_median_Incomplete_Vs_XIST <- Combine_Lsts(x='m.Bal_median_Incomplete_Vs_XIST', 
                                              y=m.Bal_median_Incomplete, 
                                              z=m.XIST_Tissue_Counts)
m.Immune_median_Incomplete_Vs_XIST <- Combine_Lsts(x='m.Immune_median_Incomplete_Vs_XIST', 
                                                 y=m.Immune_median_Incomplete, 
                                                 z=m.XIST_Tissue_Counts)

# Rename columns in each df
f.One_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                     x=f.One_median_Incomplete_Vs_XIST, 
                                     a='median_Incomplete', 
                                     b='XIST') 
f.Tuk_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                     x=f.Tuk_median_Incomplete_Vs_XIST, 
                                     a='median_Incomplete', 
                                     b='XIST') 
f.Bal_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                     x=f.Bal_median_Incomplete_Vs_XIST, 
                                     a='median_Incomplete', 
                                     b='XIST') 
f.Immune_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                        x=f.Immune_median_Incomplete_Vs_XIST, 
                                        a='median_Incomplete', 
                                        b='XIST')  

m.One_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                     x=m.One_median_Incomplete_Vs_XIST, 
                                     a='median_Incomplete',
                                     b='XIST') 
m.Tuk_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                     x=m.Tuk_median_Incomplete_Vs_XIST, 
                                     a='median_Incomplete', 
                                     b='XIST') 
m.Bal_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                     x=m.Bal_median_Incomplete_Vs_XIST, 
                                     a='median_Incomplete', 
                                     b='XIST') 
m.Immune_median_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                        x=m.Immune_median_Incomplete_Vs_XIST, 
                                        a='median_Incomplete', 
                                        b='XIST') 

# Apply lm to each df in list
Linear_Model.4 <- function(x) {
  z <- lm(median_Incomplete ~ XIST, data = x)
  return(z)
}

f.lm_One_Incomplete_XIST <- lapply(f.One_median_Incomplete_Vs_XIST, Linear_Model.4)
f.lm_Tuk_Incomplete_XIST <- lapply(f.Tuk_median_Incomplete_Vs_XIST, Linear_Model.4)
f.lm_Bal_Incomplete_XIST <- lapply(f.Bal_median_Incomplete_Vs_XIST, Linear_Model.4)
f.lm_Immune_Incomplete_XIST <- lapply(f.Immune_median_Incomplete_Vs_XIST, Linear_Model.4)

m.lm_One_Incomplete_XIST <- lapply(m.One_median_Incomplete_Vs_XIST, Linear_Model.4)
m.lm_Tuk_Incomplete_XIST <- lapply(m.Tuk_median_Incomplete_Vs_XIST, Linear_Model.4)
m.lm_Bal_Incomplete_XIST <- lapply(m.Bal_median_Incomplete_Vs_XIST, Linear_Model.4)
m.lm_Immune_Incomplete_XIST <- lapply(m.Immune_median_Incomplete_Vs_XIST, Linear_Model.4)

# Apply function to list of dfs
f.Res_One_Incomplete_XIST <- lapply(f.lm_One_Incomplete_XIST, Regression_Res)
f.Res_Tuk_Incomplete_XIST <- lapply(f.lm_Tuk_Incomplete_XIST, Regression_Res)
f.Res_Bal_Incomplete_XIST <- lapply(f.lm_Bal_Incomplete_XIST, Regression_Res)
f.Res_Immune_Incomplete_XIST <- lapply(f.lm_Immune_Incomplete_XIST, Regression_Res)

m.Res_Tuk_Incomplete_XIST <- lapply(m.lm_Tuk_Incomplete_XIST, Regression_Res)
m.Res_One_Incomplete_XIST <- lapply(m.lm_One_Incomplete_XIST, Regression_Res)
m.Res_Bal_Incomplete_XIST <- lapply(m.lm_Bal_Incomplete_XIST, Regression_Res)
m.Res_Immune_Incomplete_XIST <- lapply(m.lm_Immune_Incomplete_XIST, Regression_Res)

# Add vector as column to df
f.Regression$pval_One_Incomplete_median <- lapply(f.Res_One_Incomplete_XIST, Return_pval) 
f.Regression$R2_One_Incomplete_median <- lapply(f.Res_One_Incomplete_XIST, Return_R2) 

f.Regression$pval_Tuk_Incomplete_median <- lapply(f.Res_Tuk_Incomplete_XIST, Return_pval)
f.Regression$R2_Tuk_Incomplete_median <- lapply(f.Res_Tuk_Incomplete_XIST, Return_R2)

f.Regression$pval_Bal_Incomplete_median <- lapply(f.Res_Bal_Incomplete_XIST, Return_pval)
f.Regression$R2_Bal_Incomplete_median <- lapply(f.Res_Bal_Incomplete_XIST, Return_R2)

f.Regression$pval_Immune_Incomplete_median <- lapply(f.Res_Immune_Incomplete_XIST, Return_pval)
f.Regression$R2_Immune_Incomplete_median <- lapply(f.Res_Immune_Incomplete_XIST, Return_R2)

m.Regression$pval_One_Incomplete_median <- lapply(m.Res_One_Incomplete_XIST, Return_pval) 
m.Regression$R2_One_Incomplete_median <- lapply(m.Res_One_Incomplete_XIST, Return_R2) 

m.Regression$pval_Tuk_Incomplete_median <- lapply(m.Res_Tuk_Incomplete_XIST, Return_pval)
m.Regression$R2_Tuk_Incomplete_median <- lapply(m.Res_Tuk_Incomplete_XIST, Return_R2)

m.Regression$pval_Bal_Incomplete_median <- lapply(m.Res_Bal_Incomplete_XIST, Return_pval)
m.Regression$R2_Bal_Incomplete_median <- lapply(m.Res_Bal_Incomplete_XIST, Return_R2)

m.Regression$pval_Immune_Incomplete_median <- lapply(m.Res_Immune_Incomplete_XIST, Return_pval)
m.Regression$R2_Immune_Incomplete_median <- lapply(m.Res_Immune_Incomplete_XIST, Return_R2)

# ______________________________________________________________________________________________________________________
# Table 1; Columns 29:36
# Correlation of all genes/ all genes not evaluated / PAR genes with XIST
# ______________________________________________________________________________________________________________________
# Categories:
# "All_Evaluated_Balaton_Tukiainen", "Not_Evaluated_In_Either", "Immune_Genes_Not_Evaluated", "PAR_In_Balaton"   

# Subset list of genes from each df of x counts in list 
f.All_Eval <- Map(Filter_Func, 
                  x=f.Tissue_XCounts, 
                  y='All_Evaluated_Balaton_Tukiainen')
f.Not_Eval <- Map(Filter_Func, 
                  x=f.Tissue_XCounts, 
                  y='Not_Evaluated_In_Either')
f.Immune_Not_Eval <- Map(Filter_Func, 
                         x=f.Tissue_XCounts, 
                         y='Immune_Genes_Not_Evaluated')
f.PAR_Bal <-  Map(Filter_Func, 
                  x=f.Tissue_XCounts, 
                  y='PAR_In_Balaton')

m.All_Eval <- Map(Filter_Func, 
                  x=m.Tissue_XCounts, 
                  y='All_Evaluated_Balaton_Tukiainen')
m.Not_Eval <- Map(Filter_Func, 
                  x=m.Tissue_XCounts, 
                  y='Not_Evaluated_In_Either')
m.Immune_Not_Eval <- Map(Filter_Func, 
                         x=m.Tissue_XCounts, 
                         y='Immune_Genes_Not_Evaluated')
m.PAR_Bal <-  Map(Filter_Func, 
                  x=m.Tissue_XCounts,
                  y='PAR_In_Balaton')

# Get the median values X chm genes
f.All_Eval_median <- lapply(f.All_Eval, Median_Val)
f.Not_Eval_median <- lapply(f.Not_Eval, Median_Val)
f.Immune_Not_Eval_median <- lapply(f.Immune_Not_Eval, Median_Val)
f.PAR_median <- lapply(f.PAR_Bal, Median_Val)

m.All_Eval_median <- lapply(m.All_Eval, Median_Val)
m.Not_Eval_median <- lapply(m.Not_Eval, Median_Val)
m.Immune_Not_Eval_median <- lapply(m.Immune_Not_Eval, Median_Val)
m.PAR_median <- lapply(m.PAR_Bal, Median_Val)

# Combine list of median value of genes and XIST expression
f.All_Eval_median_Vs_XIST <- Combine_Lsts(x='f.All_Eval_median_Vs_XIST', 
                                        y=f.All_Eval_median, 
                                        z=f.XIST_Tissue_Counts)
f.Not_Eval_median_Vs_XIST <- Combine_Lsts(x='f.Not_Eval_median_Vs_XIST', 
                                        y=f.Not_Eval_median, 
                                        z=f.XIST_Tissue_Counts)
f.Immune_Not_Eval_median_Vs_XIST <- Combine_Lsts(x='f.Immune_Not_Eval_median_Vs_XIST', 
                                               y=f.Immune_Not_Eval_median, 
                                               z=f.XIST_Tissue_Counts)
f.PAR_median_Vs_XIST <- Combine_Lsts(x='f.PAR_median_Vs_XIST', y=f.PAR_median, z=f.XIST_Tissue_Counts)

m.All_Eval_median_Vs_XIST <- Combine_Lsts(x='m.All_Eval_median_Vs_XIST', 
                                        y=m.All_Eval_median, 
                                        z=m.XIST_Tissue_Counts)
m.Not_Eval_median_Vs_XIST <- Combine_Lsts(x='m.Not_Eval_median_Vs_XIST', 
                                        y=m.Not_Eval_median, 
                                        z=m.XIST_Tissue_Counts)
m.Immune_Not_Eval_median_Vs_XIST <- Combine_Lsts(x='m.Immune_Not_Eval_median_Vs_XIST', 
                                               y=m.Immune_Not_Eval_median, 
                                               z=m.XIST_Tissue_Counts)
m.PAR_median_Vs_XIST <- Combine_Lsts(x='m.PAR_median_Vs_XIST', 
                                   y=m.PAR_median, 
                                   z=m.XIST_Tissue_Counts)

# Rename columns in each df
f.All_Eval_median_Vs_XIST <- Map(Rename_Col, 
                               x=f.All_Eval_median_Vs_XIST, 
                               a='Misc', 
                               b='XIST') 
f.Not_Eval_median_Vs_XIST <- Map(Rename_Col, 
                               x=f.Not_Eval_median_Vs_XIST, 
                               a='Misc', 
                               b='XIST') 
f.Immune_Not_Eval_median_Vs_XIST <- Map(Rename_Col, 
                                      x=f.Immune_Not_Eval_median_Vs_XIST, 
                                      a='Misc', 
                                      b='XIST') 
f.PAR_median_Vs_XIST <- Map(Rename_Col, 
                          x=f.PAR_median_Vs_XIST, 
                          a='Misc',
                          b='XIST') 

m.All_Eval_median_Vs_XIST <- Map(Rename_Col, 
                               x=m.All_Eval_median_Vs_XIST, 
                               a='Misc', 
                               b='XIST') 
m.Not_Eval_median_Vs_XIST <- Map(Rename_Col, 
                               x=m.Not_Eval_median_Vs_XIST,
                               a='Misc', 
                               b='XIST') 
m.Immune_Not_Eval_median_Vs_XIST <- Map(Rename_Col, 
                                      x=m.Immune_Not_Eval_median_Vs_XIST, 
                                      a='Misc', 
                                      b='XIST') 
m.PAR_median_Vs_XIST <- Map(Rename_Col, 
                          x=m.PAR_median_Vs_XIST, 
                          a='Misc', 
                          b='XIST') 

# Apply lm to each df in list
Linear_Model.5 <- function(x) {
  z <- lm(Misc ~ XIST, data = x)
  return(z)
}

f.lm_All_Eval <- lapply(f.All_Eval_median_Vs_XIST, Linear_Model.5)
f.lm_Not_Eval <- lapply(f.Not_Eval_median_Vs_XIST, Linear_Model.5)
f.lm_Immune_Not_Eval <- lapply(f.Immune_Not_Eval_median_Vs_XIST, Linear_Model.5)
f.lm_PAR <- lapply(f.PAR_median_Vs_XIST, Linear_Model.5)

m.lm_All_Eval <- lapply(m.All_Eval_median_Vs_XIST, Linear_Model.5)
m.lm_Not_Eval <- lapply(m.Not_Eval_median_Vs_XIST, Linear_Model.5)
m.lm_Immune_Not_Eval <- lapply(m.Immune_Not_Eval_median_Vs_XIST, Linear_Model.5)
m.lm_PAR <- lapply(m.PAR_median_Vs_XIST, Linear_Model.5)

# Apply function to list of dfs
f.Res_All_Eval <- lapply(f.lm_All_Eval, Regression_Res)
f.Res_Not_Eval <- lapply(f.lm_Not_Eval, Regression_Res)
f.Res_Immune_Not_Eval <- lapply(f.lm_Immune_Not_Eval, Regression_Res)
f.Res_PAR <- lapply(f.lm_PAR, Regression_Res)

m.Res_All_Eval <- lapply(m.lm_All_Eval, Regression_Res)
m.Res_Not_Eval <- lapply(m.lm_Not_Eval, Regression_Res)
m.Res_Immune_Not_Eval <- lapply(m.lm_Immune_Not_Eval, Regression_Res)
m.Res_PAR <- lapply(m.lm_PAR, Regression_Res)

# Add vector as column to df
f.Regression$pval_All_Eval <- lapply(f.Res_All_Eval, Return_pval)
f.Regression$R2_All_Eval <- lapply(f.Res_All_Eval, Return_R2)

f.Regression$pval_Not_Eval <- lapply(f.Res_Not_Eval, Return_pval)
f.Regression$R2_Not_Eval <- lapply(f.Res_Not_Eval, Return_R2)

f.Regression$pval_Immune_Not_Eval <- lapply(f.Res_Immune_Not_Eval, Return_pval)
f.Regression$R2_Immune_Not_Eval <- lapply(f.Res_Immune_Not_Eval, Return_R2)

f.Regression$pval_PAR <- lapply(f.Res_PAR, Return_pval)
f.Regression$R2_PAR <- lapply(f.Res_PAR, Return_R2)

m.Regression$pval_All_Eval <- lapply(m.Res_All_Eval, Return_pval)
m.Regression$R2_All_Eval <- lapply(m.Res_All_Eval, Return_R2)

m.Regression$pval_Not_Eval <- lapply(m.Res_Not_Eval, Return_pval)
m.Regression$R2_Not_Eval <- lapply(m.Res_Not_Eval, Return_R2)

m.Regression$pval_Immune_Not_Eval <- lapply(m.Res_Immune_Not_Eval, Return_pval)
m.Regression$R2_Immune_Not_Eval <- lapply(m.Res_Immune_Not_Eval, Return_R2)

m.Regression$pval_PAR <- lapply(m.Res_PAR, Return_pval)
m.Regression$R2_PAR <- lapply(m.Res_PAR, Return_R2)

# ______________________________________________________________________________________________________________________
#  Write table 1 and summary tables; Columns 37:39
# ______________________________________________________________________________________________________________________
# Add col with median(XIST) and sd(XIST)
f.Regression$median_XIST <- lapply(f.XIST_Tissue_Counts, median)
m.Regression$median_XIST <- lapply(m.XIST_Tissue_Counts, median)

f.Regression$sd_XIST <- lapply(f.XIST_Tissue_Counts, sd)
m.Regression$sd_XIST <- lapply(m.XIST_Tissue_Counts, sd)

# Add column with number of tissues
Count_Rows <- function(x){
  res <- nrow(x)
  return(res)
}

Num_Fem <- as.character(lapply(f.Tissue_Lst, Count_Rows))
Num_Male <- as.character(lapply(m.Tissue_Lst, Count_Rows))

# Add column with number of samples per tissue
f.Regression <- cbind(Num_Tissues=Num_Fem, f.Regression) 
m.Regression <- cbind(Num_Tissues=Num_Male, m.Regression)

# Convert Num_Tissues from factor to numeric
f.Regression$Num_Tissues <- as.numeric(as.character(f.Regression$Num_Tissues))
m.Regression$Num_Tissues <- as.numeric(as.character(m.Regression$Num_Tissues))

# Convert rest of cols from list to numeric
f.Regression[,4:ncol(f.Regression)] <- lapply(f.Regression[,4:ncol(f.Regression)], function(x) unlist(x))
m.Regression[,4:ncol(m.Regression)] <- lapply(m.Regression[,4:ncol(m.Regression)], function(x) unlist(x))

# Add column indicating sample tissue type as first column
f.Regression <- cbind(Tissue=names(f.Tissue_Lst), f.Regression) 
m.Regression <- cbind(Tissue=names(m.Tissue_Lst), m.Regression) 

# Sort by number of tissues
f.Regression <- f.Regression[order(f.Regression$Num_Tissues, decreasing=TRUE),]
m.Regression <- m.Regression[order(m.Regression$Num_Tissues, decreasing=TRUE),]

# Write to file
write.csv(f.Regression, LM_FEM, row.names=FALSE)
write.csv(m.Regression, LM_MALE, row.names=FALSE)

# ______________________________________________________________________________________________________________________
#  Session Data
# ______________________________________________________________________________________________________________________
save.image(file=DATA)

