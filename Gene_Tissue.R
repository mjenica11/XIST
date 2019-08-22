# Perform linear regression on the gene count data (version 7).
# Find correlation between mean X chm expression and XIST per tissue across individuals. 

# Constants
COUNTS <- "~/XIST_Vs_TSIX/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct" # TPM normalized
# COUNTS <- "~/XIST_Vs_TSIX/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct" # These are unnormalized counts
METRICS <- "~/XIST_Vs_TSIX/Files/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"
PHENOTYPES <- "~/XIST_Vs_TSIX/Files/GTEX_v7_Annotations_SubjectPhenotypesDS.txt"
GENCODE <- "gencode.v19.genes.v7.patched_contigs.gff3"
GENE_LST <- "~/XIST_Vs_TSIX/Files/X_Genes_Status.json"

# Load libraries
library(readr) 
library(refGenome) # Parse .gff file
library(dplyr) 
library(data.table) 
library(stringr)
library(broom)
library(rjson)

# read.gtf will only accept relative path
setwd("~/XIST_Vs_TSIX/Files/")

# Read in files
Metrics <- read_tsv(METRICS) # Contains tissue sample info
Phenotypes <- read_tsv(PHENOTYPES) # Contains sex info
Gene_Cts <- fread(COUNTS) # df of counts
Gene_Lst <- fromJSON(file=GENE_LST)

# Create ensemblGenome object for storing Ensembl genomic annotation data
ENS <- ensemblGenome()

# Read in .gff annotation file as ensemblGenome object
read.gtf(ENS, GENCODE)

# ________________________________________________________________________________________________________
# Get list of genes on X chromosome
# ________________________________________________________________________________________________________

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
# Don't want to include XIST in mean X chromsome count values
X_Genes <- X_Genes[!X_Genes$gene_name == 'XIST',]

# ________________________________________________________________________________________________________
# Organize samples into list of dfs by tissue type
# ________________________________________________________________________________________________________

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
sapply(Tissue_Lst, function(x) sum(is.na(x)))

# Get list of female IDs
Fem.IDs <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Get list of sample sex for each tissue
Sex_Lst <- lapply(Tissue_Lst, function(x) {
  x <- with(x['ID'], ifelse(ID %in% Fem.IDs, "Female", "Male"))
})

# Add column containing sample sex to each df in list 
Tissue_Lst <- mapply(cbind, Tissue_Lst, "Sex"=Sex_Lst, SIMPLIFY=F)

# For each individual, make a data frame of samples that comes from the same person and store in list.
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

# ________________________________________________________________________________________________________
# Get X chromosome counts for each sample organized by tissue type
# ________________________________________________________________________________________________________

# For each individual, make a data frame of gene counts from samples that come from the same person and store in list.
Gene_Cts <- data.frame(Gene_Cts, stringsAsFactors = F) # was both data.table and data frame
colnames(Gene_Cts) <- str_replace_all(colnames(Gene_Cts), pattern = "\\.", replacement = "-")

# Rename columns
names(Gene_Cts)[1:2] <- c("gene_id", "gene_name")

# Remove decimals in gene IDs in gene counts df
Gene_Cts$gene_id <- sub("\\.\\d+$", "", Gene_Cts$gene_id)

# Split tissue lists by sex
f.Tissue_Lst <- lapply(Tissue_Lst, function(x){
  x[x$Sex == 'Female',]
})

m.Tissue_Lst <- lapply(Tissue_Lst, function(x){
  x[x$Sex == 'Male',]
})

# Sort count data into list of dfs by tissue 
f.Tissue_Counts <- lapply(f.Tissue_Lst, function(x){
  tmp <- which(colnames(Gene_Cts) %in% x$Sample)
  Gene_Cts[, c(1,2,tmp)]
})

m.Tissue_Counts <- lapply(m.Tissue_Lst, function(x){
  tmp <- which(colnames(Gene_Cts) %in% x$Sample)
  Gene_Cts[, c(1,2,tmp)]
})

# Check for any empty data frames
# e.g. individuals who have GTEx IDs, but have no tissue samples (no cols in df)
length(f.Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, f.Tissue_Counts)) # FALSE
length(m.Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, m.Tissue_Counts)) # FALSE

# Remove individuals who do not have count data
# i.e. only have 'gene_id' and 'gene_name' cols in df
f.Tissue_Counts <- f.Tissue_Counts[sapply(f.Tissue_Counts, function(x) dim(x)[2]) > 2] 
m.Tissue_Counts <- m.Tissue_Counts[sapply(m.Tissue_Counts, function(x) dim(x)[2]) > 2] 

# Check for any empty data frames
length(f.Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, f.Tissue_Counts)) # TRUE
length(m.Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, m.Tissue_Counts)) # TRUE

# Remove sample replicates
f.Tissue_Counts <- lapply(f.Tissue_Counts, function(x){
  x[, !(names(x) %in% Sample_Replicates)]
})

m.Tissue_Counts <- lapply(m.Tissue_Counts, function(x){
  x[, !(names(x) %in% Sample_Replicates)]
})

# Get subset of just X chm genes from each data frame in list
# Will not include XIST
f.Tissue_XCounts <- lapply(f.Tissue_Counts, function(x) {
  filter(x, gene_id %in% X_Genes$gene_id)
})

m.Tissue_XCounts <- lapply(m.Tissue_Counts, function(x) {
  filter(x, gene_id %in% X_Genes$gene_id)
})

# Check for missing values
sapply(f.Tissue_XCounts, function(x) sum(is.na(x)))
sapply(m.Tissue_XCounts, function(x) sum(is.na(x)))

# Get the mean values of the X chm genes
f.MeanX_Tissue_Counts <- lapply(f.Tissue_XCounts, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

m.MeanX_Tissue_Counts <- lapply(m.Tissue_XCounts, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Get XIST values from each data frame
f.XIST_Tissue_Counts <- lapply(f.Tissue_Counts, function(x) {
  as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
})

m.XIST_Tissue_Counts <- lapply(m.Tissue_Counts, function(x) {
  as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
})

# Sanity check
# Check that all the samples are present in both lists of dfs
length(f.XIST_Tissue_Counts) == length(f.MeanX_Tissue_Counts) # TRUE
length(m.XIST_Tissue_Counts) == length(m.MeanX_Tissue_Counts) # TRUE

# How many tissue types are there per sex?
length(f.XIST_Tissue_Counts) # 51
length(m.XIST_Tissue_Counts) # 47

# Are all the same individuals listed in both lists?
all(names(f.XIST_Tissue_Counts) == names(f.MeanX_Tissue_Counts)) # TRUE
all(names(m.XIST_Tissue_Counts) == names(m.MeanX_Tissue_Counts)) # TRUE

# Are they listed in the same order?
identical(names(f.XIST_Tissue_Counts), names(f.MeanX_Tissue_Counts)) # TRUE
identical(names(m.XIST_Tissue_Counts), names(m.MeanX_Tissue_Counts)) # TRUE

# ________________________________________________________________________________________________________
# Table 1; Column 1:
# Correlate expression from XIST with X chromosome expression within a person across all of their tissues.
# ________________________________________________________________________________________________________
# Function to combine vectors into df
Combine_Vectors <- function(a, b){
  data.frame(a, b)
}

# Combine lists of named character vectors to a list of dataframes
f.MeanX_Vs_XIST <- list()
for (i in 1:length(f.XIST_Tissue_Counts)){
  f.MeanX_Vs_XIST[[i]] <- Combine_Vectors(f.MeanX_Tissue_Counts[[i]], f.XIST_Tissue_Counts[[i]])
}
names(f.MeanX_Vs_XIST) <- names(f.XIST_Tissue_Counts)

m.MeanX_Vs_XIST <- list()
for (i in 1:length(m.XIST_Tissue_Counts)){
  m.MeanX_Vs_XIST[[i]] <- Combine_Vectors(m.MeanX_Tissue_Counts[[i]], m.XIST_Tissue_Counts[[i]])
}
names(m.MeanX_Vs_XIST) <- names(m.XIST_Tissue_Counts)

# Rename col names in each df
f.MeanX_Vs_XIST <- lapply(f.MeanX_Vs_XIST, setNames, c("MeanX", "XIST")) 
m.MeanX_Vs_XIST <- lapply(m.MeanX_Vs_XIST, setNames, c("MeanX", "XIST")) 

# Check for any missing values
sapply(f.MeanX_Vs_XIST, function(x) sum(is.na(x)))
sapply(m.MeanX_Vs_XIST, function(x) sum(is.na(x)))

# Apply lm to each df in list
lm_f.MeanX_XIST <- lapply(f.MeanX_Vs_XIST, function(x) {
  z <- lm(MeanX ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_m.MeanX_XIST <- lapply(m.MeanX_Vs_XIST, function(x) {
  z <- lm(MeanX ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Function to extract r squared values
Regression_Res <- function(lm){ # expecting object of class 'lm'
  sum. <- summary(lm)
  r.2 <- sum.$r.squared
  p. <- summary(lm)$fstatistic
  p.val <- pf(p.[1], p.[2], p.[3], lower.tail=FALSE, log.p=FALSE) #pf: F distribution ; arg 2 and 3 are the degrees of freedom 
  attributes(p.val) <- NULL
  res <- data.frame(p_val=p.val, r_2=r.2)
  return(res)
}

# Apply function to list of dfs
Res_f.MeanX_XIST <- lapply(lm_f.MeanX_XIST, Regression_Res)
Res_m.MeanX_XIST <- lapply(lm_m.MeanX_XIST, Regression_Res)

# Make table summarizing results
f.MeanX_XIST.df <- as.data.frame(do.call(rbind, Res_f.MeanX_XIST))
m.MeanX_XIST.df <- as.data.frame(do.call(rbind, Res_m.MeanX_XIST))

# Add column with sex
# Sample_Sex <- sapply(rownames(MeanX_XIST.df), function(x) {
#   if (x %in% Fem.IDs > 0) {"Female"}
#   else {"Male"}
# }, simplify = TRUE)
# 
# MeanX_XIST.df <- cbind(Sample_Sex, MeanX_XIST.df)


# Label columns
colnames(f.MeanX_XIST.df) <- c("p_Value", "R2_MeanX")
colnames(m.MeanX_XIST.df) <- c("p_Value","R2_MeanX")

# ________________________________________________________________________________________________________
# Table 1; Columns 2-6
# Correlation of all genes reported as silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Silenced_In_Both", "Silenced_In_Tukiainen", "Silenced_In_Balaton", "Silenced_In_At_Least_One", 
# "Immune_Genes_Silenced_In_At_Least_One"

# Subset list of genes silenced in both studies from each df in list by sex
Filter <- function(x, y){
  x <- filter(x, gene_name %in% Gene_Lst[[y]])
  return(x)
}

f.Both_Silenced <- Map(Filter, x=f.Tissue_Counts, y='Silenced_In_Both')
f.Tuk_Silenced <- Map(Filter, x=f.Tissue_Counts, y='Silenced_In_Tukiainen')
f.Bal_Silenced <- Map(Filter, x=f.Tissue_Counts, y='Silenced_In_Balaton')
f.One_Silenced <- Map(Filter, x=f.Tissue_Counts, y='Silenced_In_At_Least_One')
f.Immune_Silenced <- Map(Filter, x=f.Tissue_Counts, y='Immune_Genes_Silenced_In_At_Least_One')

m.Both_Silenced <- Map(Filter, x=m.Tissue_Counts, y='Silenced_In_Both')
m.Tuk_Silenced <- Map(Filter, x=m.Tissue_Counts, y='Silenced_In_Tukiainan')
m.Bal_Silenced <- Map(Filter, x=m.Tissue_Counts, y='Silenced_In_Balaton')
m.One_Silenced <- Map(Filter, x=m.Tissue_Counts, y='Silenced_In_At_Least_One')
m.Immune_Silenced <- Map(Filter, x=m.Tissue_Counts, y='Immune_Genes_Silenced_In_At_Least_One')

# Get the mean values of putatively silenced X chm genes
Mean_Val <- function(x){
  res <- colMeans(x[sapply(x, is.numeric)]) 
  return(res)
}

f.Mean_Silenced <- lapply(f.Both_Silenced, Mean_Val)
f.Tuk_Mean_Silenced <- lapply(f.Both_Silenced, Mean_Val)
f.Bal_Mean_Silenced <- lapply(f.Bal_Silenced, Mean_Val)
f.One_Mean_Silenced <- lapply(f.One_Silenced, Mean_Val)
f.Immune_Mean_Silenced <- lapply(f.Immune_Silenced, Mean_Val)

m.Mean_Silenced <- lapply(m.Both_Silenced, Mean_Val)
m.Tuk_Mean_Silenced <- lapply(m.Both_Silenced, Mean_Val)
m.Bal_Mean_Silenced <- lapply(m.Bal_Silenced, Mean_Val)
m.One_Mean_Silenced <- lapply(m.One_Silenced, Mean_Val)
m.Immune_Mean_Silenced <- lapply(m.Immune_Silenced, Mean_Val)

# Combine list of mean value of silenced genes and XIST expression
Combine_Lsts <- function(x, y, z){
  y <- list()
  for (i in 1:length(z)){
    y[[i]] <- Combine_Vectors(x[[i]], z[[i]])
  }
  names(y) <- names(z)
  return(y)
}

f.Silenced_Mean_Vs_XIST <- Combine_Lsts(x=f.Mean_Silenced, y='f.Silenced_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Tuk_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=f.Tuk_Mean_Silenced, y='f.Tuk_Silenced_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Bal_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=f.Bal_Mean_Silenced, y='f.Bal_Silenced_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)
f.One_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=f.One_Mean_Silenced, y='f.One_Silenced_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Immune_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=f.Immune_Mean_Silenced, y='f.Immune_Silenced_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)

m.Silenced_Mean_Vs_XIST <- Combine_Lsts(x=m.Mean_Silenced, y='m.Silenced_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Tuk_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=m.Tuk_Mean_Silenced, y='m.Tuk_Silenced_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Bal_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=m.Bal_Mean_Silenced, y='m.Bal_Silenced_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)
m.One_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=m.One_Mean_Silenced, y='m.One_Silenced_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Immune_Silenced_Mean_Vs_XIST <- Combine_Lsts(x=m.Immune_Mean_Silenced, y='m.Immune_Silenced_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)

# Rename columns in each df
Rename_Col <- function(x){
  x <- setNames(x, c("Mean_Silenced", "XIST"))
}

f.Silenced_Mean_Vs_XIST <- lapply(f.Silenced_Mean_Vs_XIST, Rename_Col)
f.Tuk_Silenced_Mean_Vs_XIST <- lapply(f.Tuk_Silenced_Mean_Vs_XIST, Rename_Col)
f.Bal_Silenced_Mean_Vs_XIST <- lapply(f.Bal_Silenced_Mean_Vs_XIST, Rename_Col)
f.One_Silenced_Mean_Vs_XIST <- lapply(f.One_Silenced_Mean_Vs_XIST, Rename_Col)
f.Immune_Silenced_Mean_Vs_XIST <- lapply(f.Immune_Silenced_Mean_Vs_XIST, Rename_Col)

m.Silenced_Mean_Vs_XIST <- lapply(m.Silenced_Mean_Vs_XIST, Rename_Col)
m.Tuk_Silenced_Mean_Vs_XIST <- lapply(m.Tuk_Silenced_Mean_Vs_XIST, Rename_Col)
m.Bal_Silenced_Mean_Vs_XIST <- lapply(m.Bal_Silenced_Mean_Vs_XIST, Rename_Col)
m.One_Silenced_Mean_Vs_XIST <- lapply(m.One_Silenced_Mean_Vs_XIST, Rename_Col)
m.Immune_Silenced_Mean_Vs_XIST <- lapply(m.Immune_Silenced_Mean_Vs_XIST, Rename_Col)

# Apply lm to each df in list
Linear_Model <- function(x) {
  z <- lm(Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
}

f.lm_Silenced_XIST <- lapply(f.Silenced_Mean_Vs_XIST, Linear_Model)
f.lm_Tuk_Silenced_XIST <- lapply(f.Tuk_Silenced_Mean_Vs_XIST, Linear_Model)
f.lm_Bal_Silenced_XIST <- lapply(f.Bal_Silenced_Mean_Vs_XIST, Linear_Model)
f.lm_One_Silenced_XIST <- lapply(f.One_Silenced_Mean_Vs_XIST, Linear_Model)
f.lm_Immune_Silenced_XIST <- lapply(f.Immune_Silenced_Mean_Vs_XIST, Linear_Model)

m.lm_Silenced_XIST <- lapply(m.Silenced_Mean_Vs_XIST, Linear_Model)
m.lm_Tuk_Silenced_XIST <- lapply(m.Tuk_Silenced_Mean_Vs_XIST, Linear_Model)
m.lm_Bal_Silenced_XIST <- lapply(m.Bal_Silenced_Mean_Vs_XIST, Linear_Model)
m.lm_One_Silenced_XIST <- lapply(m.One_Silenced_Mean_Vs_XIST, Linear_Model)
m.lm_Immune_Silenced_XIST <- lapply(m.Immune_Silenced_Mean_Vs_XIST, Linear_Model)

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
Unlist_Col <- function(x){
  res <- unlist(x$r_2)
  return(res)
}

f.MeanX_XIST.df$R2_Silenced_Mean <- lapply(f.Res_Silenced_XIST, Unlist_Col)
f.MeanX_XIST.df$R2_Tuk_Silenced_Mean <- lapply(f.Tuk_Res_Silenced_XIST, Unlist_Col)
f.MeanX_XIST.df$R2_Bal_Silenced_Mean <- lapply(f.Bal_Res_Silenced_XIST, Unlist_Col)
f.MeanX_XIST.df$R2_One_Silenced_Mean <- lapply(f.One_Res_Silenced_XIST, Unlist_Col)
f.MeanX_XIST.df$R2_Immune_Silenced_Mean <- lapply(f.Immune_Res_Silenced_XIST, Unlist_Col)

m.MeanX_XIST.df$R2_Silenced_Mean <- lapply(m.Res_Silenced_XIST, Unlist_Col)
m.MeanX_XIST.df$R2_Tuk_Silenced_Mean <- lapply(m.Tuk_Res_Silenced_XIST, Unlist_Col)
m.MeanX_XIST.df$R2_Bal_Silenced_Mean <- lapply(m.Bal_Res_Silenced_XIST, Unlist_Col)
m.MeanX_XIST.df$R2_One_Silenced_Mean <- lapply(m.One_Res_Silenced_XIST, Unlist_Col)
m.MeanX_XIST.df$R2_Immune_Silenced_Mean <- lapply(m.Immune_Res_Silenced_XIST, Unlist_Col)

# ________________________________________________________________________________________________________
# Table 1; Columns 7-10
# Correlation of all genes reported as variably silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Variable_In_At_Least_One", "Variable_In_Tukiainen", "Variable_In_Balaton", "Immune_Gene_Variable_In_At_Least_One"

# Subset list of genes variably silenced from each df in list 
f.One_Variable <- Map(Filter, x=f.Tissue_Counts, y='Variable_In_At_Least_One')
f.Tuk_Variable <- Map(Filter, x=f.Tissue_Counts, y='Variable_In_Tukiainen')
f.Bal_Variable <- Map(Filter, x=f.Tissue_Counts, y='Variable_In_Balaton')
f.Immune_Variable <- Map(Filter, x=f.Tissue_Counts, y='Immune_Gene_Variable_In_At_Least_One')

m.One_Variable <- Map(Filter, x=m.Tissue_Counts, y='Variable_In_At_Least_One')
m.Tuk_Variable <- Map(Filter, x=m.Tissue_Counts, y='Variable_In_Tukiainen')
m.Bal_Variable <- Map(Filter, x=m.Tissue_Counts, y='Variable_In_Balaton')
m.Immune_Variable <- Map(Filter, x=m.Tissue_Counts, y='Immune_Gene_Variable_In_At_Least_One')

# Get the mean values of putatively variably silenced X chm genes
f.One_Mean_Variable <- lapply(f.One_Variable, Mean_Val)
f.Tuk_Mean_Variable <- lapply(f.Tuk_Variable, Mean_Val)
f.Bal_Mean_Variable <- lapply(f.Bal_Variable, Mean_Val)
f.Immune_Mean_Variable <- lapply(f.Immune_Variable, Mean_Val)

m.One_Mean_Variable <- lapply(m.One_Variable, Mean_Val)
m.Tuk_Mean_Variable <- lapply(m.Tuk_Variable, Mean_Val)
m.Bal_Mean_Variable <- lapply(m.Bal_Variable, Mean_Val)
m.Immune_Mean_Variable <- lapply(m.Immune_Variable, Mean_Val)

# Combine list of mean value of variably silenced genes and XIST expression
f.One_Mean_Variable_Vs_XIST <- Combine_Lsts(x=f.One_Mean_Variable, y='f.One_Mean_Variable_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Tuk_Mean_Variable_Vs_XIST <- Combine_Lsts(x=f.Tuk_Mean_Variable, y='f.Tuk_Mean_Variable_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Bal_Mean_Variable_Vs_XIST <- Combine_Lsts(x=f.Bal_Mean_Variable, y='f.Bal_Mean_Variable_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Immune_Mean_Variable_Vs_XIST <- Combine_Lsts(x=f.Immune_Mean_Variable, y='f.Immune_Mean_Variable_Vs_XIST', z=f.XIST_Tissue_Counts)

m.One_Mean_Variable_Vs_XIST <- Combine_Lsts(x=m.One_Mean_Variable, y='m.One_Mean_Variable_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Tuk_Mean_Variable_Vs_XIST <- Combine_Lsts(x=m.Tuk_Mean_Variable, y='m.Tuk_Mean_Variable_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Bal_Mean_Variable_Vs_XIST <- Combine_Lsts(x=m.Bal_Mean_Variable, y='m.Bal_Mean_Variable_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Immune_Mean_Variable_Vs_XIST <- Combine_Lsts(x=m.Immune_Mean_Variable, y='m.Immune_Mean_Variable_Vs_XIST', z=m.XIST_Tissue_Counts)

# Rename columns in each df
f.One_Mean_Variable_Vs_XIST <- lapply(f.One_Mean_Variable_Vs_XIST, Rename_Col)
f.Tuk_Mean_Variable_Vs_XIST <- lapply(f.Tuk_Mean_Variable_Vs_XIST, Rename_Col) 
f.Bal_Mean_Variable_Vs_XIST <- lapply(f.Bal_Mean_Variable_Vs_XIST, Rename_Col)
f.Immune_Mean_Variable_Vs_XIST <- lapply(f.Immune_Mean_Variable_Vs_XIST, Rename_Col) 

m.One_Mean_Variable_Vs_XIST <- lapply(m.One_Mean_Variable_Vs_XIST, Rename_Col)
m.Tuk_Mean_Variable_Vs_XIST <- lapply(m.Tuk_Mean_Variable_Vs_XIST, Rename_Col) 
m.Bal_Mean_Variable_Vs_XIST <- lapply(m.Bal_Mean_Variable_Vs_XIST, Rename_Col)
m.Immune_Mean_Variable_Vs_XIST <- lapply(m.Immune_Mean_Variable_Vs_XIST, Rename_Col) 

# Apply lm to each df in list
f.lm_One_Variable_XIST <- lapply(f.One_Mean_Variable_Vs_XIST, Linear_Model)
f.lm_Tuk_Variable_XIST <- lapply(f.Tuk_Mean_Variable_Vs_XIST, Linear_Model)  
f.lm_Bal_Variable_XIST <- lapply(f.Bal_Mean_Variable_Vs_XIST, Linear_Model)  
f.lm_Immune_Variable_XIST <- lapply(f.Immune_Mean_Variable_Vs_XIST, Linear_Model)  

m.lm_One_Variable_XIST <- lapply(m.One_Mean_Variable_Vs_XIST, Linear_Model)
m.lm_Tuk_Variable_XIST <- lapply(m.Tuk_Mean_Variable_Vs_XIST, Linear_Model)  
m.lm_Bal_Variable_XIST <- lapply(m.Bal_Mean_Variable_Vs_XIST, Linear_Model)  
m.lm_Immune_Variable_XIST <- lapply(m.Immune_Mean_Variable_Vs_XIST, Linear_Model) 

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
f.MeanX_XIST.df$One_Variable_Mean <- lapply(f.Res_Silenced_XIST, Unlist_Col)
f.MeanX_XIST.df$Tuk_Variable_Mean <- lapply(f.Res_Tuk_Variable_XIST, Unlist_Col)
f.MeanX_XIST.df$Bal_Variable_Mean <- lapply(f.Res_Bal_Variable_XIST, Unlist_Col)
f.MeanX_XIST.df$Immune_Variable_Mean <- lapply(f.Res_Immune_Variable_XIST, Unlist_Col)

m.MeanX_XIST.df$One_Variable_Mean <- lapply(m.Res_Silenced_XIST, Unlist_Col)
m.MeanX_XIST.df$Tuk_Variable_Mean <- lapply(m.Res_Tuk_Variable_XIST, Unlist_Col)
m.MeanX_XIST.df$Bal_Variable_Mean <- lapply(m.Res_Bal_Variable_XIST, Unlist_Col)
m.MeanX_XIST.df$Immune_Variable_Mean <- lapply(m.Res_Immune_Variable_XIST, Unlist_Col)

# ________________________________________________________________________________________________________
# Table 1; Columns 11-14
# Correlation of all genes reported as incompletely silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Incomplete_In_At_Least_One", Incomplete_In_Tukiainen", "Incomplete_In_Balaton", "Immune_Genes_Incomplete_In_At_Least_One"

# Subset list of genes incompletly silenced from each df in list 
f.One_Incomplete <- Map(Filter, f.Tissue_Counts, y='Incomplete_In_At_Least_One')
f.Tuk_Incomplete <- Map(Filter, f.Tissue_Counts, y='Incomplete_In_Tukiainen')
f.Bal_Incomplete <- Map(Filter, f.Tissue_Counts, y='Incomplete_In_Balaton')
f.Immune_Incomplete <- Map(Filter, f.Tissue_Counts, y='Immune_Genes_Incomplete_In_At_Least_One')

m.One_Incomplete <- Map(Filter, m.Tissue_Counts, y='Incomplete_In_At_Least_One')
m.Tuk_Incomplete <- Map(Filter, m.Tissue_Counts, y='Incomplete_In_Tukiainen')
m.Bal_Incomplete <- Map(Filter, m.Tissue_Counts, y='Incomplete_In_Balaton')
m.Immune_Incomplete <- Map(Filter, m.Tissue_Counts, y='Immune_Genes_Incomplete_In_At_Least_One')

# Get the mean values of putatively incompletely silenced X chm genes
f.One_Mean_Incomplete <- lapply(f.One_Incomplete, Mean_Val)
f.Tuk_Mean_Incomplete <- lapply(f.Tuk_Incomplete, Mean_Val)
f.Bal_Mean_Incomplete <- lapply(f.Bal_Incomplete, Mean_Val)
f.Immune_Mean_Incomplete <- lapply(f.Immune_Incomplete, Mean_Val)

m.One_Mean_Incomplete <- lapply(m.One_Incomplete, Mean_Val)
m.Tuk_Mean_Incomplete <- lapply(m.Tuk_Incomplete, Mean_Val)
m.Bal_Mean_Incomplete <- lapply(m.Bal_Incomplete, Mean_Val)
m.Immune_Mean_Incomplete <- lapply(m.Immune_Incomplete, Mean_Val)

# Combine list of mean value of incompletely silenced genes and XIST expression
f.One_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=f.One_Mean_Incomplete, y='f.One_Mean_Incomplete_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Tuk_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=f.Tuk_Mean_Incomplete, y='f.Tuk_Mean_Incomplete_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Bal_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=f.Bal_Mean_Incomplete, y='f.Bal_Mean_Incomplete_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Immune_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=f.Immune_Mean_Incomplete, y='f.Immune_Mean_Incomplete_Vs_XIST', z=f.XIST_Tissue_Counts)

m.One_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=m.One_Mean_Incomplete, y='m.One_Mean_Incomplete_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Tuk_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=m.Tuk_Mean_Incomplete, y='m.Tuk_Mean_Incomplete_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Bal_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=m.Bal_Mean_Incomplete, y='m.Bal_Mean_Incomplete_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Immune_Mean_Incomplete_Vs_XIST <- Combine_Lsts(x=m.Immune_Mean_Incomplete, y='m.Immune_Mean_Incomplete_Vs_XIST', z=m.XIST_Tissue_Counts)

# Rename columns in each df
f.One_Mean_Incomplete_Vs_XIST <- lapply(f.One_Mean_Incomplete_Vs_XIST, Rename_Col) 
f.Tuk_Mean_Incomplete_Vs_XIST <- lapply(f.Tuk_Mean_Incomplete_Vs_XIST, Rename_Col)
f.Bal_Mean_Incomplete_Vs_XIST <- lapply(f.Bal_Mean_Incomplete_Vs_XIST, Rename_Col)
f.Immune_Mean_Incomplete_Vs_XIST <- lapply(f.Immune_Mean_Incomplete_Vs_XIST, Rename_Col) 

m.One_Mean_Incomplete_Vs_XIST <- lapply(m.One_Mean_Incomplete_Vs_XIST, Rename_Col) 
m.Tuk_Mean_Incomplete_Vs_XIST <- lapply(m.Tuk_Mean_Incomplete_Vs_XIST, Rename_Col)
m.Bal_Mean_Incomplete_Vs_XIST <- lapply(m.Bal_Mean_Incomplete_Vs_XIST, Rename_Col)
m.Immune_Mean_Incomplete_Vs_XIST <- lapply(m.Immune_Mean_Incomplete_Vs_XIST, Rename_Col) 

# Apply lm to each df in list
f.lm_One_Incomplete_XIST <- lapply(f.One_Mean_Incomplete_Vs_XIST, Linear_Model)
f.lm_Tuk_Incomplete_XIST <- lapply(f.Tuk_Mean_Incomplete_Vs_XIST, Linear_Model)
f.lm_Bal_Incomplete_XIST <- lapply(f.Bal_Mean_Incomplete_Vs_XIST, Linear_Model)
f.lm_Immune_Incomplete_XIST <- lapply(f.Immune_Mean_Incomplete_Vs_XIST, Linear_Model)

m.lm_One_Incomplete_XIST <- lapply(m.One_Mean_Incomplete_Vs_XIST, Linear_Model)
m.lm_Tuk_Incomplete_XIST <- lapply(m.Tuk_Mean_Incomplete_Vs_XIST, Linear_Model)
m.lm_Bal_Incomplete_XIST <- lapply(m.Bal_Mean_Incomplete_Vs_XIST, Linear_Model)
m.lm_Immune_Incomplete_XIST <- lapply(m.Immune_Mean_Incomplete_Vs_XIST, Linear_Model)

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
f.MeanX_XIST.df$One_Incomplete_Mean <- lapply(f.Res_One_Incomplete_XIST, Unlist_Col) 
f.MeanX_XIST.df$Tuk_Incomplete_Mean <- lapply(f.Res_Tuk_Incomplete_XIST, Unlist_Col)
f.MeanX_XIST.df$Bal_Incomplete_Mean <- lapply(f.Res_Bal_Incomplete_XIST, Unlist_Col)
f.MeanX_XIST.df$Immune_Incomplete_Mean <- lapply(f.Res_Immune_Incomplete_XIST, Unlist_Col)

m.MeanX_XIST.df$One_Incomplete_Mean <- lapply(m.Res_One_Incomplete_XIST, Unlist_Col) 
m.MeanX_XIST.df$Tuk_Incomplete_Mean <- lapply(m.Res_Tuk_Incomplete_XIST, Unlist_Col)
m.MeanX_XIST.df$Bal_Incomplete_Mean <- lapply(m.Res_Bal_Incomplete_XIST, Unlist_Col)
m.MeanX_XIST.df$Immune_Incomplete_Mean <- lapply(m.Res_Immune_Incomplete_XIST, Unlist_Col)

# ________________________________________________________________________________________________________
# Table 1; Columns 15-18
# Correlation of all genes/ all genes not evaluated / PAR genes with XIST
# ________________________________________________________________________________________________________

# Categories:
# "All_Evaluated_Balaton_Tukiainen", "Not_Evaluated_In_Either", "Immune_Genes_Not_Evaluated", "PAR_In_Balaton"   

# Subset list of genes from each df in list 
f.All_Eval <- Map(Filter, f.Tissue_Counts, y='All_Evaluated_Balaton_Tukiainen')
f.Not_Eval <- Map(Filter, f.Tissue_Counts, y='Not_Evaluated_In_Either')
f.Immune_Not_Eval <- Map(Filter, f.Tissue_Counts, y='Immune_Genes_Not_Evaluated')
f.PAR_Bal <-  Map(Filter, f.Tissue_Counts, y='PAR_In_Balaton')

m.All_Eval <- Map(Filter, m.Tissue_Counts, y='All_Evaluated_Balaton_Tukiainen')
m.Not_Eval <- Map(Filter, m.Tissue_Counts, y='Not_Evaluated_In_Either')
m.Immune_Not_Eval <- Map(Filter, m.Tissue_Counts, y='Immune_Genes_Not_Evaluated')
m.PAR_Bal <-  Map(Filter, m.Tissue_Counts, y='PAR_In_Balaton')

# Get the mean values X chm genes
f.All_Eval_Mean <- lapply(f.All_Eval, Mean_Val)
f.Not_Eval_Mean <- lapply(f.Not_Eval, Mean_Val)
f.Immune_Not_Eval_Mean <- lapply(f.Immune_Not_Eval, Mean_Val)
f.PAR_Mean <- lapply(f.PAR_Bal, Mean_Val)

m.All_Eval_Mean <- lapply(m.All_Eval, Mean_Val)
m.Not_Eval_Mean <- lapply(m.Not_Eval, Mean_Val)
m.Immune_Not_Eval_Mean <- lapply(m.Immune_Not_Eval, Mean_Val)
m.PAR_Mean <- lapply(m.PAR_Bal, Mean_Val)

# Combine list of mean value of genes and XIST expression
f.All_Eval_Mean_Vs_XIST <- Combine_Lsts(x=f.All_Eval_Mean, y='f.All_Eval_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Not_Eval_Mean_Vs_XIST <- Combine_Lsts(x=f.Not_Eval_Mean, y='f.Not_Eval_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)
f.Immune_Not_Eval_Mean_Vs_XIST <- Combine_Lsts(x=f.Immune_Not_Eval_Mean, y='f.Immune_Not_Eval_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)
f.PAR_Mean_Vs_XIST <- Combine_Lsts(x=f.PAR_Mean, y='f.PAR_Mean_Vs_XIST', z=f.XIST_Tissue_Counts)

m.All_Eval_Mean_Vs_XIST <- Combine_Lsts(x=m.All_Eval_Mean, y='m.All_Eval_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Not_Eval_Mean_Vs_XIST <- Combine_Lsts(x=m.Not_Eval_Mean, y='m.Not_Eval_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)
m.Immune_Not_Eval_Mean_Vs_XIST <- Combine_Lsts(x=m.Immune_Not_Eval_Mean, y='m.Immune_Not_Eval_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)
m.PAR_Mean_Vs_XIST <- Combine_Lsts(x=m.PAR_Mean, y='m.PAR_Mean_Vs_XIST', z=m.XIST_Tissue_Counts)

# Rename columns in each df
f.All_Eval_Mean_Vs_XIST <- lapply(f.All_Eval_Mean_Vs_XIST, Rename_Col) 
f.Not_Eval_Mean_Vs_XIST <- lapply(f.Not_Eval_Mean_Vs_XIST, Rename_Col) 
f.Immune_Not_Eval_Mean_Vs_XIST <- lapply(f.Immune_Not_Eval_Mean_Vs_XIST, Rename_Col) 
f.PAR_Mean_Vs_XIST <- lapply(f.PAR_Mean_Vs_XIST, Rename_Col) 

m.All_Eval_Mean_Vs_XIST <- lapply(m.All_Eval_Mean_Vs_XIST, Rename_Col) 
m.Not_Eval_Mean_Vs_XIST <- lapply(m.Not_Eval_Mean_Vs_XIST, Rename_Col) 
m.Immune_Not_Eval_Mean_Vs_XIST <- lapply(m.Immune_Not_Eval_Mean_Vs_XIST, Rename_Col) 
m.PAR_Mean_Vs_XIST <- lapply(m.PAR_Mean_Vs_XIST, Rename_Col) 

# Apply lm to each df in list
f.lm_All_Eval <- lapply(f.All_Eval_Mean_Vs_XIST, Linear_Model)
f.lm_Not_Eval <- lapply(f.Not_Eval_Mean_Vs_XIST, Linear_Model)
f.lm_Immune_Not_Eval <- lapply(f.Immune_Not_Eval_Mean_Vs_XIST, Linear_Model)
f.lm_PAR <- lapply(f.PAR_Mean_Vs_XIST, Linear_Model)

m.lm_All_Eval <- lapply(m.All_Eval_Mean_Vs_XIST, Linear_Model)
m.lm_Not_Eval <- lapply(m.Not_Eval_Mean_Vs_XIST, Linear_Model)
m.lm_Immune_Not_Eval <- lapply(m.Immune_Not_Eval_Mean_Vs_XIST, Linear_Model)
m.lm_PAR <- lapply(m.PAR_Mean_Vs_XIST, Linear_Model)

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
f.MeanX_XIST.df$All_Eval <- lapply(f.Res_All_Eval, Unlist_Col)
f.MeanX_XIST.df$Not_Eval <- lapply(f.Res_Not_Eval, Unlist_Col)
f.MeanX_XIST.df$Immune_Not_Eval <- lapply(f.Res_Immune_Not_Eval, Unlist_Col)
f.MeanX_XIST.df$PAR <- lapply(f.Res_PAR, Unlist_Col)

m.MeanX_XIST.df$All_Eval <- lapply(m.Res_All_Eval, Unlist_Col)
m.MeanX_XIST.df$Not_Eval <- lapply(m.Res_Not_Eval, Unlist_Col)
m.MeanX_XIST.df$Immune_Not_Eval <- lapply(m.Res_Immune_Not_Eval, Unlist_Col)
m.MeanX_XIST.df$PAR <- lapply(m.Res_PAR, Unlist_Col)

# ________________________________________________________________________________________________________
#  Write table 1 and summary 
# ________________________________________________________________________________________________________
# Average R^2 of silenced genes reported in both studies for females and males
f.df1 <- f.MeanX_XIST.df %>% summarise(Silenced = mean(Silenced_Mean))
m.df1 <- m.MeanX_XIST.df %>% summarise(Silenced = mean(Silenced_Mean))

f.df2 <- f.MeanX_XIST.df %>% summarise(Tuk_Silenced = mean(Tuk_Silenced_Mean))
m.df2 <- m.MeanX_XIST.df %>% summarise(Tuk_Silenced = mean(Tuk_Silenced_Mean))

f.df3 <- f.MeanX_XIST.df %>% summarise(Bal_Silenced = mean(Bal_Silenced_Mean))
m.df3 <- m.MeanX_XIST.df %>% summarise(Bal_Silenced = mean(Bal_Silenced_Mean))

f.df4 <- f.MeanX_XIST.df %>% summarise(One_Silenced = mean(One_Silenced_Mean))
m.df4 <- m.MeanX_XIST.df %>% summarise(One_Silenced = mean(One_Silenced_Mean))

f.df5 <- f.MeanX_XIST.df %>% summarise(Immune_Silenced = mean(Immune_Silenced_Mean))
m.df5 <- m.MeanX_XIST.df %>% summarise(Immune_Silenced = mean(Immune_Silenced_Mean))

# Average R^2 of silenced genes reported in both studies for females and males
f.df6 <- f.MeanX_XIST.df %>% summarise(One_Variable = mean(One_Variable_Mean))
m.df6 <- m.MeanX_XIST.df %>% summarise(One_Variable = mean(One_Variable_Mean))

f.df7 <- f.MeanX_XIST.df %>% summarise(Tuk_Variable = mean(Tuk_Variable_Mean))
m.df7 <- m.MeanX_XIST.df %>% summarise(Tuk_Variable = mean(Tuk_Variable_Mean))

f.df8 <- f.MeanX_XIST.df %>% summarise(Bal_Variable = mean(Bal_Variable_Mean))
m.df8 <- m.MeanX_XIST.df %>% summarise(Bal_Variable = mean(Bal_Variable_Mean))

f.df9 <- f.MeanX_XIST.df %>% summarise(Immune_Variable = mean(Immune_Variable_Mean))
m.df9 <- m.MeanX_XIST.df %>% summarise(Immune_Variable = mean(Immune_Variable_Mean))

# Average R^2 of silenced genes reported in both studies for females and males
f.df10 <- f.MeanX_XIST.df %>% summarise(One_Incomplete = mean(One_Incomplete_Mean))
m.df10 <- m.MeanX_XIST.df %>% summarise(One_Incomplete = mean(One_Incomplete_Mean))

f.df11 <- f.MeanX_XIST.df %>% summarise(Tuk_Incomplete = mean(Tuk_Incomplete_Mean))
m.df11 <- m.MeanX_XIST.df %>% summarise(Tuk_Incomplete = mean(Tuk_Incomplete_Mean))

f.df12 <- f.MeanX_XIST.df %>% summarise(Bal_Incomplete = mean(Bal_Incomplete_Mean))
m.df12 <- m.MeanX_XIST.df %>% summarise(Bal_Incomplete = mean(Bal_Incomplete_Mean))

f.df13 <- f.MeanX_XIST.df %>% summarise(Immune_Incomplete = mean(Immune_Incomplete_Mean))
m.df13 <- m.MeanX_XIST.df %>% summarise(Immune_Incomplete = mean(Immune_Incomplete_Mean))

# Average R^2 of genes for females and males
f.df14 <- f.MeanX_XIST.df %>% summarise(All_Eval = mean(All_Eval))
m.df14 <- m.MeanX_XIST.df %>% summarise(All_Eval = mean(All_Eval))

f.df15 <- f.MeanX_XIST.df %>% summarise(Not_Eval = mean(Not_Eval))
m.df15 <- m.MeanX_XIST.df %>% summarise(Not_Eval = mean(Not_Eval))

f.df16 <- f.MeanX_XIST.df %>% summarise(Immune_Not_Eval = mean(Immune_Not_Eval))
m.df16 <- m.MeanX_XIST.df %>% summarise(Immune_Not_Eval = mean(Immune_Not_Eval))

f.df17 <- f.MeanX_XIST.df %>% summarise(Bal_PAR = mean(PAR))
m.df17 <- m.MeanX_XIST.df %>% summarise(Bal_PAR = mean(PAR))

# Add column with number of tissues
Num_Tissues <- function(x){
  res <- length(x$Sample)
  return(res)
}

Num_Fem <- lapply(f.Tissue_Lst, Num_Tissues)
Num_Male <- lapply(m.Tissue_Lst, Num_Tissues)

# Drop tissues females/males do not have
f.remove <- c("Prostate", "Testis", "Cells - Leukemia cell line (CML)")
m.remove <- c("Ovary", "Uterus", "Vagina", "Fallopian Tube", "Cervix - Ectocervix", "Cervix - Endocervix", "Cells - Leukemia cell line (CML)")

Num_Fem <- Num_Fem[!names(Num_Fem) %in% f.remove]
Num_Male <- Num_Male[!names(Num_Male) %in% m.remove]

f.MeanX_XIST.df$Num_Tissues <- Num_Fem
m.MeanX_XIST.df$Num_Tissues <- Num_Male

# Add col with mean(XIST) and sd(XIST)
f.MeanX_XIST.df$Mean_XIST <- lapply(f.XIST_Tissue_Counts, mean)
m.MeanX_XIST.df$Mean_XIST <- lapply(m.XIST_Tissue_Counts, mean)

f.MeanX_XIST.df$sd_XIST <- lapply(f.XIST_Tissue_Counts, sd)
m.MeanX_XIST.df$sd_XIST <- lapply(m.XIST_Tissue_Counts, sd)

#f.MeanX_XIST.df <- apply(f.MeanX_XIST.df, 2, as.character)
#m.MeanX_XIST.df <- apply(m.MeanX_XIST.df, 2, as.character)

# Combine dfs of averages to make quick summary 
f.df_lst <- list(f.df1, f.df2, f.df3, f.df4, f.df5, f.df6, f.df7, f.df8, f.df9, f.df10, f.df11, f.df12, f.df13, f.df14, f.df15, f.df16, f.df17)
m.df_lst <- list(m.df1, m.df2, m.df3, m.df4, m.df5, m.df6, m.df7, m.df8, m.df9, m.df10, m.df11, m.df12, m.df13, m.df14, m.df15, m.df16, m.df17)

f.Combined <- Reduce(merge, lapply(f.df_lst, function(x) data.frame(x, rn = row.names(x))))
m.Combined <- Reduce(merge, lapply(m.df_lst, function(x) data.frame(x, rn = row.names(x))))

# Write to file
# write.csv(f.Combined, "Fem.Tissue_Linear_Models_Summary.csv")
# write.csv(m.Combined, "m.Tissue_Linear_Models_Summary.csv")
# write.csv(f.MeanX_XIST.df, "Fem.Tissue_Correlations.csv")
# write.csv(m.MeanX_XIST.df, "m.Tissue_Correlations.csv")

# Table of slopes
# Get list of female and male tissue type samples
f.Tissues <- rownames(f.MeanX_XIST.df)
m.Tissues <- rownames(m.MeanX_XIST.df)

# Extract list of slopes
f.Slopes <- lapply(lm_f.MeanX_XIST, function(x){
  res <- x[[1]][[2]]
  return(res)
})

m.Slopes <- lapply(lm_m.MeanX_XIST, function(x){
  res <- x[[1]][[2]]
  return(res)
})

# Make and write df
l <- list(f.Slopes, m.Slopes)
Slopes.df <- rbindlist(l, use.names=TRUE, fill=TRUE, idcol="Sex")
Slopes.df$Sex <- c("Female", "Male")

#write.csv(Slopes.df, "Tissue_Slopes_Table.csv")

# ________________________________________________________________________________________________________
#  TO-DO: adjust Regression_Res function to also extract p-values
# ________________________________________________________________________________________________________
# Function to extract r squared  and p values values
### QUESTION: Do I want to set lower.tail to F?
Regression_Res <- function(lm){ # expecting object of class 'lm'
  sum. <- summary(lm)
  r.2 <- sum.$r.squared
  p. <- summary(lm)$fstatistic
  p.val <- pf(p.[1], p.[2], p.[3], lower.tail=FALSE, log.p=FALSE) #pf: F distribution ; arg 2 and 3 are the degrees of freedom 
  attributes(p.val) <- NULL
  res <- data.frame(p_val=p.val, r_2=r.2)
  return(res)
}

# Apply function to list of dfs
Res_f.MeanX_XIST <- lapply(lm_f.MeanX_XIST, Regression_Res)
Res_m.MeanX_XIST <- lapply(lm_m.MeanX_XIST, Regression_Res)

# ________________________________________________________________________________________________________
#  Scatter Plots
# ________________________________________________________________________________________________________
# Find the max x/y values to set as x/y lim
Max_Func <- function(x, y){
  lim <- max(x[[y]])
  return(lim)
}

# x limits
f.xmax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_XIST, y='XIST'))))
m.xmax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_XIST, y='XIST')))
# y limits
f.ymax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_XIST, y='MeanX'))))
m.ymax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_XIST, y='MeanX')))

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, XMAX, YMAX, RESULTS){
  plot(LM$model$XIST, LM$model$MeanX, main=TITLE, xlab='XIST', ylab='Mean X chromosome',
       xlim=c(0, XMAX), ylim=c(0, YMAX))
  legend("bottomright", bty="n", legend=paste("R^2: ", format(RESULTS$r_2, digits=3), "; p_Val: ", format(RESULTS$p_val, digits=3)))
  abline(LM)
}

# Print plots
#pdf('Fem.Tissue_Scatter_Plots.pdf')
f.Scatter <- Map(Scatter_Func, LM=lm_f.MeanX_XIST, TITLE=names(lm_f.MeanX_XIST), XMAX=f.xmax, YMAX=f.ymax, RESULTS=Res_f.MeanX_XIST)
#dev.off()

#pdf('m.Tissue_Scatter_Plots.pdf')
# set x lim to male max(XIST) but keep y lim as female max(MeanX)
m.Scatter <- Map(Scatter_Func, LM=lm_m.MeanX_XIST, TITLE=names(lm_m.MeanX_XIST), XMAX=m.xmax, YMAX=f.ymax, RESULTS=Res_m.MeanX_XIST)
#dev.off()

# ________________________________________________________________________________________________________
#  Violin Plots
# ________________________________________________________________________________________________________
library(ggplot2)
library(plyr)

# Create a Violin plot
ggplot(f.MeanX_Vs_XIST$`Adipose - Subcutaneous`, aes(x =XIST, y = MeanX)) + 
  geom_violin()


# Collapse list of dfs into one df
f.MeanX_Vs_XIST <- ldply(f.MeanX_Vs_XIST, data.frame)

# Violin plot function
Violin <- function(x){
  ggplot(x, aes(x =XIST, y = MeanX)) + 
    geom_violin()
}
test <- lapply(f.MeanX_Vs_XIST, Violin)


