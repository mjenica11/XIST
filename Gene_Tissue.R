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
Female_IDs <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Get list of sample sex for each tissue
Sex_Lst <- lapply(Tissue_Lst, function(x) {
  x <- with(x['ID'], ifelse(ID %in% Female_IDs, "Female", "Male"))
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
Fem_Tissue_Lst <- lapply(Tissue_Lst, function(x){
  x[x$Sex == 'Female',]
})

Male_Tissue_Lst <- lapply(Tissue_Lst, function(x){
  x[x$Sex == 'Male',]
})

# Sort count data into list of dfs by tissue 
Fem_Tissue_Counts <- lapply(Fem_Tissue_Lst, function(x){
  tmp <- which(colnames(Gene_Cts) %in% x$Sample)
  Gene_Cts[, c(1,2,tmp)]
})

Male_Tissue_Counts <- lapply(Male_Tissue_Lst, function(x){
  tmp <- which(colnames(Gene_Cts) %in% x$Sample)
  Gene_Cts[, c(1,2,tmp)]
})

# Check for any empty data frames
# e.g. individuals who have GTEx IDs, but have no tissue samples (no cols in df)
length(Fem_Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, Fem_Tissue_Counts)) # FALSE
length(Male_Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, Male_Tissue_Counts)) # FALSE

# Remove individuals who do not have count data
# i.e. only have 'gene_id' and 'gene_name' cols in df
Fem_Tissue_Counts <- Fem_Tissue_Counts[sapply(Fem_Tissue_Counts, function(x) dim(x)[2]) > 2] 
Male_Tissue_Counts <- Male_Tissue_Counts[sapply(Male_Tissue_Counts, function(x) dim(x)[2]) > 2] 

# Check for any empty data frames
length(Fem_Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, Fem_Tissue_Counts)) # TRUE
length(Male_Tissue_Counts) == length(Filter(function(x) dim(x)[2] > 2, Male_Tissue_Counts)) # TRUE

# Remove sample replicates
Fem_Tissue_Counts <- lapply(Fem_Tissue_Counts, function(x){
  x[, !(names(x) %in% Sample_Replicates)]
})

Male_Tissue_Counts <- lapply(Male_Tissue_Counts, function(x){
  x[, !(names(x) %in% Sample_Replicates)]
})

# Get subset of just X chm genes from each data frame in list
# Will not include XIST
Fem_Tissue_XCounts <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_id %in% X_Genes$gene_id)
})

Male_Tissue_XCounts <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_id %in% X_Genes$gene_id)
})

# Check for missing values
sapply(Fem_Tissue_XCounts, function(x) sum(is.na(x)))
sapply(Male_Tissue_XCounts, function(x) sum(is.na(x)))

# Get the mean values of the X chm genes
Fem_MeanX_Tissue_Counts <- lapply(Fem_Tissue_XCounts, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Male_MeanX_Tissue_Counts <- lapply(Male_Tissue_XCounts, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Get XIST values from each data frame
Fem_XIST_Tissue_Counts <- lapply(Fem_Tissue_Counts, function(x) {
  as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
})

Male_XIST_Tissue_Counts <- lapply(Male_Tissue_Counts, function(x) {
  as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
})

# ________________________________________________________________________________________________________
# Table 1; Column 1:
# Correlate expression from XIST with X chromosome expression within a person across all of their tissues.
# ________________________________________________________________________________________________________

# Check that all the samples are present in both lists of dfs
length(Fem_XIST_Tissue_Counts) == length(Fem_MeanX_Tissue_Counts) # TRUE
length(Male_XIST_Tissue_Counts) == length(Male_MeanX_Tissue_Counts) # TRUE

# How many tissue types are there per sex?
length(Fem_XIST_Tissue_Counts) # 51
length(Male_XIST_Tissue_Counts) # 47

# Are all the same individuals listed in both lists?
all(names(Fem_XIST_Tissue_Counts) == names(Fem_MeanX_Tissue_Counts)) # TRUE
all(names(Male_XIST_Tissue_Counts) == names(Male_MeanX_Tissue_Counts)) # TRUE

# Are they listed in the same order?
identical(names(Fem_XIST_Tissue_Counts), names(Fem_MeanX_Tissue_Counts)) # TRUE
identical(names(Male_XIST_Tissue_Counts), names(Male_MeanX_Tissue_Counts)) # TRUE

# Function to combine vectors into df
Combine_Vectors <- function(a, b){
  data.frame(a, b)
}

# Combine lists of named character vectors to a list of dataframes
Fem_MeanX_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_MeanX_Vs_XIST[[i]] <- Combine_Vectors(Fem_MeanX_Tissue_Counts[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_MeanX_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_MeanX_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_MeanX_Vs_XIST[[i]] <- Combine_Vectors(Male_MeanX_Tissue_Counts[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_MeanX_Vs_XIST) <- names(Male_XIST_Tissue_Counts)

# Rename col names in each df
Fem_MeanX_Vs_XIST <- lapply(Fem_MeanX_Vs_XIST, setNames, c("MeanX", "XIST")) 

Male_MeanX_Vs_XIST <- lapply(Male_MeanX_Vs_XIST, setNames, c("MeanX", "XIST")) 

# Check for any missing values
sapply(Fem_MeanX_Vs_XIST, function(x) sum(is.na(x)))

sapply(Male_MeanX_Vs_XIST, function(x) sum(is.na(x)))

# Apply lm to each df in list
lm_Fem_MeanX_XIST <- lapply(Fem_MeanX_Vs_XIST, function(x) {
  z <- lm(MeanX ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Male_MeanX_XIST <- lapply(Male_MeanX_Vs_XIST, function(x) {
  z <- lm(MeanX ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Function to extract r squared values
R_Squared <- function(lm){ # expecting object of class 'lm'
  sum. <- summary(lm)
  r.2 <- sum.$r.squared
  r.2
}

# Apply function to list of dfs
Res_Fem_MeanX_XIST <- lapply(lm_Fem_MeanX_XIST, R_Squared)
Res_Male_MeanX_XIST <- lapply(lm_Male_MeanX_XIST, R_Squared)

# Make table summarizing results
Fem_MeanX_XIST.df <- as.data.frame(do.call(rbind, Res_Fem_MeanX_XIST))
Male_MeanX_XIST.df <- as.data.frame(do.call(rbind, Res_Male_MeanX_XIST))

# Add column with sex
# Sample_Sex <- sapply(rownames(MeanX_XIST.df), function(x) {
#   if (x %in% Female_IDs > 0) {"Female"}
#   else {"Male"}
# }, simplify = TRUE)
# 
# MeanX_XIST.df <- cbind(Sample_Sex, MeanX_XIST.df)


# Label columns
colnames(Fem_MeanX_XIST.df) <- c("R2_MeanX")
colnames(Male_MeanX_XIST.df) <- c("R2_MeanX")

# Average R^2 for females vs males (probs not that informative)
Fem_MeanX_XIST.df %>% summarise(mean = mean(R2_MeanX)) 
Male_MeanX_XIST.df %>% summarise(mean = mean(R2_MeanX))


# ________________________________________________________________________________________________________
# Table 1; Columns 2-6
# Correlation of all genes reported as silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Silenced_In_Both", "Silenced_In_Tukiainen", "Silenced_In_Balaton", "Silenced_In_At_Least_One", 
# "Immune_Genes_Silenced_In_At_Least_One"

# Subset list of genes silenced in both studies from each df in list 
# Silenced in both 
Fem_Both_Silenced <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Both)
})
Male_Both_Silenced <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Both)
})

#  Silenced in Tukainen
Fem_Tuk_Silenced <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Tukiainen)
})
Male_Tuk_Silenced <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Tukiainen)
})

# Silenced in Balaton
Fem_Bal_Silenced <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Balaton)
})
Male_Bal_Silenced <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Balaton)
})

# Silenced in at least one study
Fem_One_Silenced <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_At_Least_One)
})
Male_One_Silenced <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_At_Least_One)
})

# Silenced immune genes
Fem_Immune_Silenced <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Silenced_In_At_Least_One)
})
Male_Immune_Silenced <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Silenced_In_At_Least_One)
})

# Get the mean values of putatively silenced X chm genes
Fem_Mean_Silenced <- lapply(Fem_Both_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Mean_Silenced <- lapply(Male_Both_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Tuk_Mean_Silenced <- lapply(Fem_Tuk_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Tuk_Mean_Silenced <- lapply(Male_Tuk_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Bal_Mean_Silenced <- lapply(Fem_Bal_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Bal_Mean_Silenced <- lapply(Male_Bal_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_One_Mean_Silenced <- lapply(Fem_One_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_One_Mean_Silenced <- lapply(Male_One_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Immune_Mean_Silenced <- lapply(Fem_Immune_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Immune_Mean_Silenced <- lapply(Male_Immune_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of silenced genes and XIST expression
Fem_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_Mean_Silenced[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Silenced_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_Mean_Silenced[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Silenced_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Tuk_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Tuk_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_Tuk_Mean_Silenced[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Tuk_Silenced_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Tuk_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Tuk_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_Tuk_Mean_Silenced[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Tuk_Silenced_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Bal_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Bal_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_Bal_Mean_Silenced[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Bal_Silenced_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Bal_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Bal_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_Bal_Mean_Silenced[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Bal_Silenced_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_One_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_One_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_One_Mean_Silenced[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_One_Silenced_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_One_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_One_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_One_Mean_Silenced[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_One_Silenced_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Immune_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Immune_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_Immune_Mean_Silenced[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Immune_Silenced_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Immune_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Immune_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_Immune_Mean_Silenced[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Immune_Silenced_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)

# Rename columns in each df
Fem_Silenced_Mean_Vs_XIST <- lapply(Fem_Silenced_Mean_Vs_XIST, setNames, c("Mean_Silenced", "XIST")) 
Male_Silenced_Mean_Vs_XIST <- lapply(Male_Silenced_Mean_Vs_XIST, setNames, c("Mean_Silenced", "XIST")) 

Fem_Tuk_Silenced_Mean_Vs_XIST <- lapply(Fem_Tuk_Silenced_Mean_Vs_XIST, setNames, c("Tuk_Mean_Silenced", "XIST")) 
Male_Tuk_Silenced_Mean_Vs_XIST <- lapply(Male_Tuk_Silenced_Mean_Vs_XIST, setNames, c("Tuk_Mean_Silenced", "XIST")) 

Fem_Bal_Silenced_Mean_Vs_XIST <- lapply(Fem_Bal_Silenced_Mean_Vs_XIST, setNames, c("Bal_Mean_Silenced", "XIST")) 
Male_Bal_Silenced_Mean_Vs_XIST <- lapply(Male_Bal_Silenced_Mean_Vs_XIST, setNames, c("Bal_Mean_Silenced", "XIST")) 

Fem_One_Silenced_Mean_Vs_XIST <- lapply(Fem_One_Silenced_Mean_Vs_XIST, setNames, c("One_Mean_Silenced", "XIST")) 
Male_One_Silenced_Mean_Vs_XIST <- lapply(Male_One_Silenced_Mean_Vs_XIST, setNames, c("One_Mean_Silenced", "XIST")) 

Fem_Immune_Silenced_Mean_Vs_XIST <- lapply(Fem_Immune_Silenced_Mean_Vs_XIST, setNames, c("Immune_Mean_Silenced", "XIST"))
Male_Immune_Silenced_Mean_Vs_XIST <- lapply(Male_Immune_Silenced_Mean_Vs_XIST, setNames, c("Immune_Mean_Silenced", "XIST"))

# Apply lm to each df in list
Fem_lm_Silenced_XIST <- lapply(Fem_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Silenced_XIST <- lapply(Male_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})


Fem_lm_Tuk_Silenced_XIST <- lapply(Fem_Tuk_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Tuk_Silenced_XIST <- lapply(Male_Tuk_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})


Fem_lm_Bal_Silenced_XIST <- lapply(Fem_Bal_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Bal_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Bal_Silenced_XIST <- lapply(Male_Bal_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Bal_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})


Fem_lm_One_Silenced_XIST <- lapply(Fem_One_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(One_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_One_Silenced_XIST <- lapply(Male_One_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(One_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})


Fem_lm_Immune_Silenced_XIST <- lapply(Fem_Immune_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Immune_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Immune_Silenced_XIST <- lapply(Male_Immune_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Immune_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Fem_Res_Silenced_XIST <- lapply(Fem_lm_Silenced_XIST, R_Squared)
Male_Res_Silenced_XIST <- lapply(Male_lm_Silenced_XIST, R_Squared)

Fem_Tuk_Res_Silenced_XIST <- lapply(Fem_lm_Tuk_Silenced_XIST, R_Squared)
Male_Tuk_Res_Silenced_XIST <- lapply(Male_lm_Tuk_Silenced_XIST, R_Squared)

Fem_Bal_Res_Silenced_XIST <- lapply(Fem_lm_Bal_Silenced_XIST, R_Squared)
Male_Bal_Res_Silenced_XIST <- lapply(Male_lm_Bal_Silenced_XIST, R_Squared)

Fem_One_Res_Silenced_XIST <- lapply(Fem_lm_One_Silenced_XIST, R_Squared)
Male_One_Res_Silenced_XIST <- lapply(Male_lm_One_Silenced_XIST, R_Squared)

Fem_Immune_Res_Silenced_XIST <- lapply(Fem_lm_Immune_Silenced_XIST, R_Squared)
Male_Immune_Res_Silenced_XIST <- lapply(Male_lm_Immune_Silenced_XIST, R_Squared)

# Add vector as column to df
Fem_MeanX_XIST.df$Silenced_Mean <- unlist(Fem_Res_Silenced_XIST)
Male_MeanX_XIST.df$Silenced_Mean <- unlist(Male_Res_Silenced_XIST)

Fem_MeanX_XIST.df$Tuk_Silenced_Mean <- unlist(Fem_Tuk_Res_Silenced_XIST)
Male_MeanX_XIST.df$Tuk_Silenced_Mean <- unlist(Male_Tuk_Res_Silenced_XIST)

Fem_MeanX_XIST.df$Bal_Silenced_Mean <- unlist(Fem_Bal_Res_Silenced_XIST)
Male_MeanX_XIST.df$Bal_Silenced_Mean <- unlist(Male_Bal_Res_Silenced_XIST)

Fem_MeanX_XIST.df$One_Silenced_Mean <- unlist(Fem_One_Res_Silenced_XIST)
Male_MeanX_XIST.df$One_Silenced_Mean <- unlist(Male_One_Res_Silenced_XIST)

Fem_MeanX_XIST.df$Immune_Silenced_Mean <- unlist(Fem_Immune_Res_Silenced_XIST)
Male_MeanX_XIST.df$Immune_Silenced_Mean <- unlist(Male_Immune_Res_Silenced_XIST)

# ________________________________________________________________________________________________________
# Table 1; Columns 7-10
# Correlation of all genes reported as variably silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Variable_In_At_Least_One", "Variable_In_Tukiainen", "Variable_In_Balaton", "Immune_Gene_Variable_In_At_Least_One"

# Subset list of genes variably silenced from each df in list 
Fem_One_Variable <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_At_Least_One)
})
Male_One_Variable <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_At_Least_One)
})

Fem_Tuk_Variable <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_Tukiainen)
})
Male_Tuk_Variable <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_Tukiainen)
})

Fem_Bal_Variable <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_Balaton)
})
Male_Bal_Variable <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_Balaton)
})

Fem_Immune_Variable <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Gene_Variable_In_At_Least_One)
})
Male_Immune_Variable <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Gene_Variable_In_At_Least_One)
})

# Get the mean values of putatively variably silenced X chm genes
Fem_One_Mean_Variable <- lapply(Fem_One_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_One_Mean_Variable <- lapply(Male_One_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Tuk_Mean_Variable <- lapply(Fem_Tuk_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Tuk_Mean_Variable <- lapply(Male_Tuk_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Bal_Mean_Variable <- lapply(Fem_Bal_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Bal_Mean_Variable <- lapply(Male_Bal_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Immune_Mean_Variable <- lapply(Fem_Immune_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Immune_Mean_Variable <- lapply(Male_Immune_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of variably silenced genes and XIST expression
Fem_One_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_One_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Fem_One_Mean_Variable[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_One_Mean_Variable_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_One_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_One_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Male_One_Mean_Variable[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_One_Mean_Variable_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Tuk_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Tuk_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Fem_Tuk_Mean_Variable[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Tuk_Mean_Variable_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Tuk_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Tuk_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Male_Tuk_Mean_Variable[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Tuk_Mean_Variable_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Bal_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Bal_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Fem_Bal_Mean_Variable[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Bal_Mean_Variable_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Bal_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Bal_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Male_Bal_Mean_Variable[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Bal_Mean_Variable_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Immune_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Immune_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Fem_Immune_Mean_Variable[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Immune_Mean_Variable_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Immune_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Immune_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Male_Immune_Mean_Variable[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Immune_Mean_Variable_Vs_XIST) <- names(Male_XIST_Tissue_Counts)

# Rename columns in each df
Fem_One_Mean_Variable_Vs_XIST <- lapply(Fem_One_Mean_Variable_Vs_XIST, setNames, c("One_Variable_Mean", "XIST")) 
Male_One_Mean_Variable_Vs_XIST <- lapply(Male_One_Mean_Variable_Vs_XIST, setNames, c("One_Variable_Mean", "XIST"))

Fem_Tuk_Mean_Variable_Vs_XIST <- lapply(Fem_Tuk_Mean_Variable_Vs_XIST, setNames, c("Tuk_Variable_Mean", "XIST")) 
Male_Tuk_Mean_Variable_Vs_XIST <- lapply(Male_Tuk_Mean_Variable_Vs_XIST, setNames, c("Tuk_Variable_Mean", "XIST"))

Fem_Bal_Mean_Variable_Vs_XIST <- lapply(Fem_Bal_Mean_Variable_Vs_XIST, setNames, c("Bal_Variable_Mean", "XIST"))
Male_Bal_Mean_Variable_Vs_XIST <- lapply(Male_Bal_Mean_Variable_Vs_XIST, setNames, c("Bal_Variable_Mean", "XIST"))

Fem_Immune_Mean_Variable_Vs_XIST <- lapply(Fem_Immune_Mean_Variable_Vs_XIST, setNames, c("Immune_Variable_Mean", "XIST")) 
Male_Immune_Mean_Variable_Vs_XIST <- lapply(Male_Immune_Mean_Variable_Vs_XIST, setNames, c("Immune_Variable_Mean", "XIST"))

# Apply lm to each df in list
Fem_lm_One_Variable_XIST <- lapply(Fem_One_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(One_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_One_Variable_XIST <- lapply(Male_One_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(One_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Tuk_Variable_XIST <- lapply(Fem_Tuk_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Tuk_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Tuk_Variable_XIST <- lapply(Male_Tuk_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Tuk_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Bal_Variable_XIST <- lapply(Fem_Bal_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Bal_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Bal_Variable_XIST <- lapply(Male_Bal_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Bal_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Immune_Variable_XIST <- lapply(Fem_Immune_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Immune_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Immune_Variable_XIST <- lapply(Male_Immune_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Immune_Variable_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Fem_Res_One_Variable_XIST <- lapply(Fem_lm_One_Variable_XIST, R_Squared)
Male_Res_One_Variable_XIST <- lapply(Male_lm_One_Variable_XIST, R_Squared)

Fem_Res_Tuk_Variable_XIST <- lapply(Fem_lm_Tuk_Variable_XIST, R_Squared)
Male_Res_Tuk_Variable_XIST <- lapply(Male_lm_Tuk_Variable_XIST, R_Squared)

Fem_Res_Bal_Variable_XIST <- lapply(Fem_lm_Bal_Variable_XIST, R_Squared)
Male_Res_Bal_Variable_XIST <- lapply(Male_lm_Bal_Variable_XIST, R_Squared)

Fem_Res_Immune_Variable_XIST <- lapply(Fem_lm_Immune_Variable_XIST, R_Squared)
Male_Res_Immune_Variable_XIST <- lapply(Male_lm_Immune_Variable_XIST, R_Squared)

# Add vector as column to df
Fem_MeanX_XIST.df$One_Variable_Mean <- unlist(Fem_Res_One_Variable_XIST)
Male_MeanX_XIST.df$One_Variable_Mean <- unlist(Male_Res_One_Variable_XIST)

Fem_MeanX_XIST.df$Tuk_Variable_Mean <- unlist(Fem_Res_Tuk_Variable_XIST)
Male_MeanX_XIST.df$Tuk_Variable_Mean <- unlist(Male_Res_Tuk_Variable_XIST)

Fem_MeanX_XIST.df$Bal_Variable_Mean <- unlist(Fem_Res_Bal_Variable_XIST)
Male_MeanX_XIST.df$Bal_Variable_Mean <- unlist(Male_Res_Bal_Variable_XIST)

Fem_MeanX_XIST.df$Immune_Variable_Mean <- unlist(Fem_Res_Immune_Variable_XIST)
Male_MeanX_XIST.df$Immune_Variable_Mean <- unlist(Male_Res_Immune_Variable_XIST)

# ________________________________________________________________________________________________________
# Table 1; Columns 11-14
# Correlation of all genes reported as incompletely silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Incomplete_In_At_Least_One", Incomplete_In_Tukiainen", "Incomplete_In_Balaton", "Immune_Genes_Incomplete_In_At_Least_One"

# Subset list of genes incompletly silenced from each df in list 
Fem_One_Incomplete <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_At_Least_One)
})
Male_One_Incomplete <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_At_Least_One)
})

Fem_Tuk_Incomplete <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_Tukiainen)
})
Male_Tuk_Incomplete <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_Tukiainen)
})

Fem_Bal_Incomplete <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_Balaton)
})
Male_Bal_Incomplete <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_Balaton)
})

Fem_Immune_Incomplete <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Incomplete_In_At_Least_One)
})
Male_Immune_Incomplete <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Incomplete_In_At_Least_One)
})

# Get the mean values of putatively incompletely silenced X chm genes
Fem_One_Mean_Incomplete <- lapply(Fem_One_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_One_Mean_Incomplete <- lapply(Male_One_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Tuk_Mean_Incomplete <- lapply(Fem_Tuk_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Tuk_Mean_Incomplete <- lapply(Male_Tuk_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Bal_Mean_Incomplete <- lapply(Fem_Bal_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Bal_Mean_Incomplete <- lapply(Male_Bal_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Immune_Mean_Incomplete <- lapply(Fem_Immune_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Immune_Mean_Incomplete <- lapply(Male_Immune_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of incompletely silenced genes and XIST expression
Fem_One_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_One_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Fem_One_Mean_Incomplete[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_One_Mean_Incomplete_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_One_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_One_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Male_One_Mean_Incomplete[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_One_Mean_Incomplete_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Tuk_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Tuk_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Fem_Tuk_Mean_Incomplete[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Tuk_Mean_Incomplete_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Tuk_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Tuk_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Male_Tuk_Mean_Incomplete[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Tuk_Mean_Incomplete_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Bal_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Bal_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Fem_Bal_Mean_Incomplete[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Bal_Mean_Incomplete_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Bal_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Bal_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Male_Bal_Mean_Incomplete[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Bal_Mean_Incomplete_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Immune_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Immune_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Fem_Immune_Mean_Incomplete[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Immune_Mean_Incomplete_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Immune_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Immune_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Male_Immune_Mean_Incomplete[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Immune_Mean_Incomplete_Vs_XIST) <- names(Male_XIST_Tissue_Counts)

# Rename columns in each df
Fem_One_Mean_Incomplete_Vs_XIST <- lapply(Fem_One_Mean_Incomplete_Vs_XIST, setNames, c("One_Mean_Incomplete", "XIST")) 
Male_One_Mean_Incomplete_Vs_XIST <- lapply(Male_One_Mean_Incomplete_Vs_XIST, setNames, c("One_Mean_Incomplete", "XIST"))

Fem_Tuk_Mean_Incomplete_Vs_XIST <- lapply(Fem_Tuk_Mean_Incomplete_Vs_XIST, setNames, c("Tuk_Mean_Incomplete", "XIST"))
Male_Tuk_Mean_Incomplete_Vs_XIST <- lapply(Male_Tuk_Mean_Incomplete_Vs_XIST, setNames, c("Tuk_Mean_Incomplete", "XIST"))

Fem_Bal_Mean_Incomplete_Vs_XIST <- lapply(Fem_Bal_Mean_Incomplete_Vs_XIST, setNames, c("Bal_Mean_Incomplete", "XIST"))
Male_Bal_Mean_Incomplete_Vs_XIST <- lapply(Male_Bal_Mean_Incomplete_Vs_XIST, setNames, c("Bal_Mean_Incomplete", "XIST")) 

Fem_Immune_Mean_Incomplete_Vs_XIST <- lapply(Fem_Immune_Mean_Incomplete_Vs_XIST, setNames, c("Immune_Mean_Incomplete", "XIST")) 
Male_Immune_Mean_Incomplete_Vs_XIST <- lapply(Male_Immune_Mean_Incomplete_Vs_XIST, setNames, c("Immune_Mean_Incomplete", "XIST")) 

# Apply lm to each df in list
Fem_lm_One_Incomplete_XIST <- lapply(Fem_One_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(One_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_One_Incomplete_XIST <- lapply(Male_One_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(One_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Tuk_Incomplete_XIST <- lapply(Fem_Tuk_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Tuk_Incomplete_XIST <- lapply(Male_Tuk_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Bal_Incomplete_XIST <- lapply(Fem_Bal_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Bal_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Bal_Incomplete_XIST <- lapply(Male_Bal_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Bal_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Immune_Incomplete_XIST <- lapply(Fem_Immune_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Immune_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Immune_Incomplete_XIST <- lapply(Male_Immune_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Immune_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Fem_Res_One_Incomplete_XIST <- lapply(Fem_lm_One_Incomplete_XIST, R_Squared)
Male_Res_One_Incomplete_XIST <- lapply(Male_lm_One_Incomplete_XIST, R_Squared)

Fem_Res_Tuk_Incomplete_XIST <- lapply(Fem_lm_Tuk_Incomplete_XIST, R_Squared)
Male_Res_Tuk_Incomplete_XIST <- lapply(Male_lm_Tuk_Incomplete_XIST, R_Squared)

Fem_Res_Bal_Incomplete_XIST <- lapply(Fem_lm_Bal_Incomplete_XIST, R_Squared)
Male_Res_Bal_Incomplete_XIST <- lapply(Male_lm_Bal_Incomplete_XIST, R_Squared)

Fem_Res_Immune_Incomplete_XIST <- lapply(Fem_lm_Immune_Incomplete_XIST, R_Squared)
Male_Res_Immune_Incomplete_XIST <- lapply(Male_lm_Immune_Incomplete_XIST, R_Squared)

# Add vector as column to df
Fem_MeanX_XIST.df$One_Incomplete_Mean <- unlist(Fem_Res_One_Incomplete_XIST)
Male_MeanX_XIST.df$One_Incomplete_Mean <- unlist(Male_Res_One_Incomplete_XIST)

Fem_MeanX_XIST.df$Tuk_Incomplete_Mean <- unlist(Fem_Res_Tuk_Incomplete_XIST)
Male_MeanX_XIST.df$Tuk_Incomplete_Mean <- unlist(Male_Res_Tuk_Incomplete_XIST)

Fem_MeanX_XIST.df$Bal_Incomplete_Mean <- unlist(Fem_Res_Bal_Incomplete_XIST)
Male_MeanX_XIST.df$Bal_Incomplete_Mean <- unlist(Male_Res_Bal_Incomplete_XIST)

Fem_MeanX_XIST.df$Immune_Incomplete_Mean <- unlist(Fem_Res_Immune_Incomplete_XIST)
Male_MeanX_XIST.df$Immune_Incomplete_Mean <- unlist(Male_Res_Immune_Incomplete_XIST)

# ________________________________________________________________________________________________________
# Table 1; Columns 15-18
# Correlation of all genes/ all genes not evaluated / PAR genes with XIST
# ________________________________________________________________________________________________________

# Categories:
# "All_Evaluated_Balaton_Tukiainen", "Not_Evaluated_In_Either", "Immune_Genes_Not_Evaluated", "PAR_In_Balaton"   

# Subset list of genes from each df in list 
Fem_All_Eval <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$All_Evaluated_Balaton_Tukiainen)
})
Male_All_Eval <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$All_Evaluated_Balaton_Tukiainen)
})

Fem_Not_Eval <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Not_Evaluated_In_Either)
})
Male_Not_Eval <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Not_Evaluated_In_Either)
})

Fem_Immune_Not_Eval <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Not_Evaluated)
})
Male_Immune_Not_Eval <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Not_Evaluated)
})

Fem_PAR_Bal <- lapply(Fem_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$PAR_In_Balaton)
})
Male_PAR_Bal <- lapply(Male_Tissue_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$PAR_In_Balaton)
})

# Get the mean values X chm genes
Fem_All_Eval_Mean <- lapply(Fem_All_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_All_Eval_Mean <- lapply(Male_All_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Not_Eval_Mean <- lapply(Fem_Not_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Not_Eval_Mean <- lapply(Male_Not_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_Immune_Not_Eval_Mean <- lapply(Fem_Immune_Not_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_Immune_Not_Eval_Mean <- lapply(Male_Immune_Not_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Fem_PAR_Mean <- lapply(Fem_PAR_Bal, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})
Male_PAR_Mean <- lapply(Male_PAR_Bal, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of genes and XIST expression
Fem_All_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_All_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_All_Eval_Mean[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_All_Eval_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_All_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_All_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_All_Eval_Mean[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_All_Eval_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Not_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Not_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_Not_Eval_Mean[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Not_Eval_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Not_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Not_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_Not_Eval_Mean[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Not_Eval_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_Immune_Not_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_Immune_Not_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_Immune_Not_Eval_Mean[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_Immune_Not_Eval_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_Immune_Not_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_Immune_Not_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_Immune_Not_Eval_Mean[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_Immune_Not_Eval_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)


Fem_PAR_Mean_Vs_XIST <- list()
for (i in 1:length(Fem_XIST_Tissue_Counts)){
  Fem_PAR_Mean_Vs_XIST[[i]] <- Combine_Vectors(Fem_PAR_Mean[[i]], Fem_XIST_Tissue_Counts[[i]])
}
names(Fem_PAR_Mean_Vs_XIST) <- names(Fem_XIST_Tissue_Counts)

Male_PAR_Mean_Vs_XIST <- list()
for (i in 1:length(Male_XIST_Tissue_Counts)){
  Male_PAR_Mean_Vs_XIST[[i]] <- Combine_Vectors(Male_PAR_Mean[[i]], Male_XIST_Tissue_Counts[[i]])
}
names(Male_PAR_Mean_Vs_XIST) <- names(Male_XIST_Tissue_Counts)

# Rename columns in each df
Fem_All_Eval_Mean_Vs_XIST <- lapply(Fem_All_Eval_Mean_Vs_XIST, setNames, c("All_Eval_Mean", "XIST")) 
Male_All_Eval_Mean_Vs_XIST <- lapply(Male_All_Eval_Mean_Vs_XIST, setNames, c("All_Eval_Mean", "XIST")) 

Fem_Not_Eval_Mean_Vs_XIST <- lapply(Fem_Not_Eval_Mean_Vs_XIST, setNames, c("Not_Eval_Mean", "XIST")) 
Male_Not_Eval_Mean_Vs_XIST <- lapply(Male_Not_Eval_Mean_Vs_XIST, setNames, c("Not_Eval_Mean", "XIST"))

Fem_Immune_Not_Eval_Mean_Vs_XIST <- lapply(Fem_Immune_Not_Eval_Mean_Vs_XIST, setNames, c("Immune_Not_Eval_Mean", "XIST")) 
Male_Immune_Not_Eval_Mean_Vs_XIST <- lapply(Male_Immune_Not_Eval_Mean_Vs_XIST, setNames, c("Immune_Not_Eval_Mean", "XIST")) 

Fem_PAR_Mean_Vs_XIST <- lapply(Fem_PAR_Mean_Vs_XIST, setNames, c("PAR_Mean", "XIST")) 
Male_PAR_Mean_Vs_XIST <- lapply(Male_PAR_Mean_Vs_XIST, setNames, c("PAR_Mean", "XIST")) 

# Apply lm to each df in list
Fem_lm_All_Eval <- lapply(Fem_All_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(All_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_All_Eval <- lapply(Male_All_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(All_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Not_Eval <- lapply(Fem_Not_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(Not_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Not_Eval <- lapply(Male_Not_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(Not_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_Immune_Not_Eval <- lapply(Fem_Immune_Not_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(Immune_Not_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_Immune_Not_Eval <- lapply(Male_Immune_Not_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(Immune_Not_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

Fem_lm_PAR <- lapply(Fem_PAR_Mean_Vs_XIST, function(x) {
  z <- lm(PAR_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})
Male_lm_PAR <- lapply(Male_PAR_Mean_Vs_XIST, function(x) {
  z <- lm(PAR_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Fem_Res_All_Eval <- lapply(Fem_lm_All_Eval, R_Squared)
Male_Res_All_Eval <- lapply(Male_lm_All_Eval, R_Squared)

Fem_Res_Not_Eval <- lapply(Fem_lm_Not_Eval, R_Squared)
Male_Res_Not_Eval <- lapply(Male_lm_Not_Eval, R_Squared)

Fem_Res_Immune_Not_Eval <- lapply(Fem_lm_Immune_Not_Eval, R_Squared)
Male_Res_Immune_Not_Eval <- lapply(Male_lm_Immune_Not_Eval, R_Squared)

Fem_Res_PAR <- lapply(Fem_lm_PAR, R_Squared)
Male_Res_PAR <- lapply(Male_lm_PAR, R_Squared)

# Add vector as column to df
Fem_MeanX_XIST.df$All_Eval <- unlist(Fem_Res_All_Eval)
Male_MeanX_XIST.df$All_Eval <- unlist(Male_Res_All_Eval)

Fem_MeanX_XIST.df$Not_Eval <- unlist(Fem_Res_Not_Eval)
Male_MeanX_XIST.df$Not_Eval <- unlist(Male_Res_Not_Eval)

Fem_MeanX_XIST.df$Immune_Not_Eval <- unlist(Fem_Res_Immune_Not_Eval)
Male_MeanX_XIST.df$Immune_Not_Eval <- unlist(Male_Res_Immune_Not_Eval)

Fem_MeanX_XIST.df$PAR <- unlist(Fem_Res_PAR)
Male_MeanX_XIST.df$PAR <- unlist(Male_Res_PAR)

# ________________________________________________________________________________________________________
#  Write table 1 and summary 
# ________________________________________________________________________________________________________

# Average R^2 of silenced genes reported in both studies for females and males
f.df1 <- Fem_MeanX_XIST.df %>% summarise(Silenced = mean(Silenced_Mean))
m.df1 <- Male_MeanX_XIST.df %>% summarise(Silenced = mean(Silenced_Mean))

f.df2 <- Fem_MeanX_XIST.df %>% summarise(Tuk_Silenced = mean(Tuk_Silenced_Mean))
m.df2 <- Male_MeanX_XIST.df %>% summarise(Tuk_Silenced = mean(Tuk_Silenced_Mean))

f.df3 <- Fem_MeanX_XIST.df %>% summarise(Bal_Silenced = mean(Bal_Silenced_Mean))
m.df3 <- Male_MeanX_XIST.df %>% summarise(Bal_Silenced = mean(Bal_Silenced_Mean))

f.df4 <- Fem_MeanX_XIST.df %>% summarise(One_Silenced = mean(One_Silenced_Mean))
m.df4 <- Male_MeanX_XIST.df %>% summarise(One_Silenced = mean(One_Silenced_Mean))

f.df5 <- Fem_MeanX_XIST.df %>% summarise(Immune_Silenced = mean(Immune_Silenced_Mean))
m.df5 <- Male_MeanX_XIST.df %>% summarise(Immune_Silenced = mean(Immune_Silenced_Mean))

# Average R^2 of silenced genes reported in both studies for females and males
f.df6 <- Fem_MeanX_XIST.df %>% summarise(One_Variable = mean(One_Variable_Mean))
m.df6 <- Male_MeanX_XIST.df %>% summarise(One_Variable = mean(One_Variable_Mean))

f.df7 <- Fem_MeanX_XIST.df %>% summarise(Tuk_Variable = mean(Tuk_Variable_Mean))
m.df7 <- Male_MeanX_XIST.df %>% summarise(Tuk_Variable = mean(Tuk_Variable_Mean))

f.df8 <- Fem_MeanX_XIST.df %>% summarise(Bal_Variable = mean(Bal_Variable_Mean))
m.df8 <- Male_MeanX_XIST.df %>% summarise(Bal_Variable = mean(Bal_Variable_Mean))

f.df9 <- Fem_MeanX_XIST.df %>% summarise(Immune_Variable = mean(Immune_Variable_Mean))
m.df9 <- Male_MeanX_XIST.df %>% summarise(Immune_Variable = mean(Immune_Variable_Mean))

# Average R^2 of silenced genes reported in both studies for females and males
f.df10 <- Fem_MeanX_XIST.df %>% summarise(One_Incomplete = mean(One_Incomplete_Mean))
m.df10 <- Male_MeanX_XIST.df %>% summarise(One_Incomplete = mean(One_Incomplete_Mean))

f.df11 <- Fem_MeanX_XIST.df %>% summarise(Tuk_Incomplete = mean(Tuk_Incomplete_Mean))
m.df11 <- Male_MeanX_XIST.df %>% summarise(Tuk_Incomplete = mean(Tuk_Incomplete_Mean))

f.df12 <- Fem_MeanX_XIST.df %>% summarise(Bal_Incomplete = mean(Bal_Incomplete_Mean))
m.df12 <- Male_MeanX_XIST.df %>% summarise(Bal_Incomplete = mean(Bal_Incomplete_Mean))

f.df13 <- Fem_MeanX_XIST.df %>% summarise(Immune_Incomplete = mean(Immune_Incomplete_Mean))
m.df13 <- Male_MeanX_XIST.df %>% summarise(Immune_Incomplete = mean(Immune_Incomplete_Mean))

# Average R^2 of genes for females and males
f.df14 <- Fem_MeanX_XIST.df %>% summarise(All_Eval = mean(All_Eval))
m.df14 <- Male_MeanX_XIST.df %>% summarise(All_Eval = mean(All_Eval))

f.df15 <- Fem_MeanX_XIST.df %>% summarise(Not_Eval = mean(Not_Eval))
m.df15 <- Male_MeanX_XIST.df %>% summarise(Not_Eval = mean(Not_Eval))

f.df16 <- Fem_MeanX_XIST.df %>% summarise(Immune_Not_Eval = mean(Immune_Not_Eval))
m.df16 <- Male_MeanX_XIST.df %>% summarise(Immune_Not_Eval = mean(Immune_Not_Eval))

f.df17 <- Fem_MeanX_XIST.df %>% summarise(Bal_PAR = mean(PAR))
m.df17 <- Male_MeanX_XIST.df %>% summarise(Bal_PAR = mean(PAR))

# Add column with number of tissues
Num_Tissues <- function(x){
  res <- length(x$Sample)
  return(res)
}

Num_Fem <- lapply(Fem_Tissue_Lst, Num_Tissues)
Num_Male <- lapply(Male_Tissue_Lst, Num_Tissues)

# Drop tissues females/males do not have
f.remove <- c("Prostate", "Testis", "Cells - Leukemia cell line (CML)")
m.remove <- c("Ovary", "Uterus", "Vagina", "Fallopian Tube", "Cervix - Ectocervix", "Cervix - Endocervix", "Cells - Leukemia cell line (CML)")

Num_Fem <- Num_Fem[!names(Num_Fem) %in% f.remove]
Num_Male <- Num_Male[!names(Num_Male) %in% m.remove]

Fem_MeanX_XIST.df$Num_Tissues <- Num_Fem
Male_MeanX_XIST.df$Num_Tissues <- Num_Male

# Add col with mean(XIST) and sd(XIST)
Fem_MeanX_XIST.df$Mean_XIST <- lapply(Fem_XIST_Tissue_Counts, mean)
Male_MeanX_XIST.df$Mean_XIST <- lapply(Male_XIST_Tissue_Counts, mean)

Fem_MeanX_XIST.df$sd_XIST <- lapply(Fem_XIST_Tissue_Counts, sd)
Male_MeanX_XIST.df$sd_XIST <- lapply(Male_XIST_Tissue_Counts, sd)

Fem_MeanX_XIST.df <- apply(Fem_MeanX_XIST.df, 2, as.character)
Male_MeanX_XIST.df <- apply(Male_MeanX_XIST.df, 2, as.character)

# Combine dfs of averages to make quick summary 
f.df_lst <- list(f.df1, f.df2, f.df3, f.df4, f.df5, f.df6, f.df7, f.df8, f.df9, f.df10, f.df11, f.df12, f.df13, f.df14, f.df15, f.df16, f.df17)
m.df_lst <- list(m.df1, m.df2, m.df3, m.df4, m.df5, m.df6, m.df7, m.df8, m.df9, m.df10, m.df11, m.df12, m.df13, m.df14, m.df15, m.df16, m.df17)

f.Combined <- Reduce(merge, lapply(f.df_lst, function(x) data.frame(x, rn = row.names(x))))
m.Combined <- Reduce(merge, lapply(m.df_lst, function(x) data.frame(x, rn = row.names(x))))

# Write to file
# write.csv(f.Combined, "Female_Tissue_Linear_Models_Summary.csv")
# write.csv(m.Combined, "Male_Tissue_Linear_Models_Summary.csv")
# write.csv(Fem_MeanX_XIST.df, "Female_Tissue_Correlations.csv")
# write.csv(Male_MeanX_XIST.df, "Male_Tissue_Correlations.csv")

# Table of slopes
# Get list of female and male tissue type samples
Fem_Tissues <- rownames(Fem_MeanX_XIST.df)
Male_Tissues <- rownames(Male_MeanX_XIST.df)

# Extract list of slopes
Fem_Slopes <- lapply(lm_Fem_MeanX_XIST, function(x){
  res <- x[[1]][[2]]
  return(res)
})

Male_Slopes <- lapply(lm_Male_MeanX_XIST, function(x){
  res <- x[[1]][[2]]
  return(res)
})

# Make and write df
l <- list(Fem_Slopes, Male_Slopes)
Slopes.df <- rbindlist(l, use.names=TRUE, fill=TRUE, idcol="Sex")
Slopes.df$Sex <- c("Female", "Male")

#write.csv(Slopes.df, "Tissue_Slopes_Table.csv")

# ________________________________________________________________________________________________________
#  TO-DO: adjust R_Squared function to also extract p-values
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
Res_Fem_MeanX_XIST <- lapply(lm_Fem_MeanX_XIST, Regression_Res)
Res_Male_MeanX_XIST <- lapply(lm_Male_MeanX_XIST, Regression_Res)

# ________________________________________________________________________________________________________
#  Scatter Plots
# ________________________________________________________________________________________________________
# Find the max x/y values to set as x/y lim
Max_Func <- function(x, y){
  lim <- max(x[[y]])
  return(lim)
}

# x limits
f.xmax <- round(max(unlist(Map(Max_Func, x=Fem_MeanX_Vs_XIST, y='XIST'))))
m.xmax <- max(unlist(Map(Max_Func, x=Male_MeanX_Vs_XIST, y='XIST')))
# y limits
f.ymax <- round(max(unlist(Map(Max_Func, x=Fem_MeanX_Vs_XIST, y='MeanX'))))
m.ymax <- max(unlist(Map(Max_Func, x=Male_MeanX_Vs_XIST, y='MeanX')))

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, XMAX, YMAX, RESULTS){
  plot(LM$model$XIST, LM$model$MeanX, main=TITLE, xlab='XIST', ylab='Mean X chromosome',
       xlim=c(0, XMAX), ylim=c(0, YMAX))
  legend("bottomright", bty="n", legend=paste("R^2: ", format(RESULTS$r_2, digits=3), "; p_Val: ", format(RESULTS$p_val, digits=3)))
  abline(LM)
}

# Print plots
#pdf('Female_Tissue_Scatter_Plots.pdf')
f.Scatter <- Map(Scatter_Func, LM=lm_Fem_MeanX_XIST, TITLE=names(lm_Fem_MeanX_XIST), XMAX=f.xmax, YMAX=f.ymax, RESULTS=Res_Fem_MeanX_XIST)
#dev.off()

#pdf('Male_Tissue_Scatter_Plots.pdf')
# set x lim to male max(XIST) but keep y lim as female max(MeanX)
m.Scatter <- Map(Scatter_Func, LM=lm_Male_MeanX_XIST, TITLE=names(lm_Male_MeanX_XIST), XMAX=m.xmax, YMAX=f.ymax, RESULTS=Res_Male_MeanX_XIST)
#dev.off()

# ________________________________________________________________________________________________________
#  Violin Plots
# ________________________________________________________________________________________________________
library(ggplot2)
library(Hmisc)

# test violin plot
ggplot(Fem_MeanX_Vs_XIST$`Adipose - Subcutaneous`, aes(x=XIST, y=MeanX)) + 
  geom_violin() + stat_summary(fun.y=mean, geom="point", shape=23, size=2) + # add mean point
  stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 ) +
  stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")

# Violin plot function
### idk if I need this step
Drop_Rownames <- function(x){
  rownames(x) <- c()
  return(x)
}
f.Violin_MeanX_XIST <- lapply(Fem_MeanX_Vs_XIST, Drop_Rownames)

# Have to merge list of dfs to one df
f.plotData <- do.call("rbind", f.Violin_MeanX_XIST)

Violin_Func <- function(x){
  ggplot(x, aes(x=XIST, y=MeanX)) + 
    geom_violin()  
}
Map(Violin_Func, x=f.Violin_MeanX_XIST)




