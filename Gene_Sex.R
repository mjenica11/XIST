# Perform linear regression on the gene count data (version 7).
#  Find correlation between mean X chm expression and XIST across all tissues per individual

# Constants
COUNTS <- "~/XIST_Vs_TSIX/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct"
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
# Tissue summary
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

# ________________________________________________________________________________________________________
# X chromosome counts
# ________________________________________________________________________________________________________

# For each individual, make a data frame of gene counts from samples that come from the same person and store in list.
Gene_Cts <- data.frame(Gene_Cts, stringsAsFactors = F) # was both data.table and data frame
colnames(Gene_Cts) <- str_replace_all(colnames(Gene_Cts), pattern = "\\.", replacement = "-")

# Rename columns
names(Gene_Cts)[1:2] <- c("gene_id", "gene_name")

# Remove decimals in gene IDs in gene counts df
Gene_Cts$gene_id <- sub("\\.\\d+$", "", Gene_Cts$gene_id)

# Sort count data into list of dfs by individual IDs
Ind_Counts <- lapply(Individual_IDs, function(x) {
  # Get col index that contains individuals GTEx ID
  ind <- str_which(colnames(Gene_Cts), pattern = paste("^", x, sep = ""))
  Gene_Cts[, c(1,2,ind)] # keep gene IDs and gene name columns
})
names(Ind_Counts) <- Individual_IDs

# Check for any empty data frames
# e.g. individuals who have GTEx IDs, but have no tissue samples (no cols in df)
length(Ind_Counts) == length(Filter(function(x) dim(x)[2] > 2, Ind_Counts)) # FALSE

# Remove individuals who do not have count data
# i.e. only have 'gene_id' and 'gene_name' cols in df
Ind_Counts <- Ind_Counts[sapply(Ind_Counts, function(x) dim(x)[2]) > 2] 

# Check for any empty data frames
length(Ind_Counts) == length(Filter(function(x) dim(x)[2] > 2, Ind_Counts)) # TRUE

# Get list of sample replicates
Sample_Replicates <- lapply(Dup_Ind_Tissues, function(x) {
  if (as.character(x['Sample']) != "character(0)") {
    as.character(x['Sample'])
  } else {}
})

Sample_Replicates <- unlist(Sample_Replicates)

# Remove sample replicates
Ind_Counts <- lapply(Ind_Counts, function(x){
  x[, !(names(x) %in% Sample_Replicates)]
})

# Get subset of just X chm genes from each data frame in list
# Will not include XIST
X_Ind_Counts <- lapply(Ind_Counts, function(x) {
  filter(x, gene_id %in% X_Genes$gene_id)
})

# Check for missing values
sapply(X_Ind_Counts, function(x) sum(is.na(x)))

# Get the mean values of the X chm genes
MeanX_Ind_Counts <- lapply(X_Ind_Counts, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Get XIST values from each data frame
XIST_Ind_Counts <- lapply(Ind_Counts, function(x) {
  as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
})

# ________________________________________________________________________________________________________
# Table 1; Column 1:
# Correlate expression from XIST with X chromosome expression within a person across all of their tissues.
# ________________________________________________________________________________________________________

# Check that all the samples are present in both lists of dfs
length(XIST_Ind_Counts) == length(MeanX_Ind_Counts) # TRUE

# How many individuals are there?
length(XIST_Ind_Counts) # 714

# Are all the same individuals listed in both lists?
all(names(XIST_Ind_Counts) == names(MeanX_Ind_Counts)) # TRUE

# Are they listed in the same order?
identical(names(XIST_Ind_Counts), names(MeanX_Ind_Counts)) # TRUE

# Function to combine vectors into df
Combine_Vectors <- function(a, b){
  data.frame(a, b)
}

# Combine lists of named character vectors to a list of dataframes
MeanX_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  MeanX_Vs_XIST[[i]] <- Combine_Vectors(MeanX_Ind_Counts[[i]], XIST_Ind_Counts[[i]])
}
names(MeanX_Vs_XIST) <- names(XIST_Ind_Counts)

# Rename col names in each df
MeanX_Vs_XIST <- lapply(MeanX_Vs_XIST, setNames, c("MeanX", "XIST")) 

# Check for any missing values
sapply(MeanX_Vs_XIST, function(x) sum(is.na(x)))

# Apply lm to each df in list
lm_MeanX_XIST <- lapply(MeanX_Vs_XIST, function(x) {
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
Res_MeanX_XIST <- lapply(lm_MeanX_XIST, R_Squared)

# Get list of female IDs
Female_IDs <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Make table summarizing results
MeanX_XIST.df <- as.data.frame(do.call(rbind, Res_MeanX_XIST))

# Add column with sex
Sample_Sex <- sapply(rownames(MeanX_XIST.df), function(x) {
  if (x %in% Female_IDs > 0) {"Female"}
  else {"Male"}
}, simplify = TRUE)

MeanX_XIST.df <- cbind(Sample_Sex, MeanX_XIST.df)

# Label columns
colnames(MeanX_XIST.df) <- c("Sex", "R2_MeanX")

# Average R^2 for females vs males
MeanX_XIST.df %>% group_by(Sex) %>% summarise(mean = mean(R2_MeanX))

# ________________________________________________________________________________________________________
# Table 1; Columns 2-6
# Correlation of all genes reported as silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Silenced_In_Both", "Silenced_In_Tukiainen", "Silenced_In_Balaton", "Silenced_In_At_Least_One", 
# "Immune_Genes_Silenced_In_At_Least_One"

# Subset list of genes silenced in both studies from each df in list 
Both_Silenced <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Both)
})

Tuk_Silenced <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Tukiainen)
})

Bal_Silenced <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Balaton)
})

One_Silenced <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_At_Least_One)
})

Immune_Silenced <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Silenced_In_At_Least_One)
})

# Get the mean values of putatively silenced X chm genes
Mean_Silenced <- lapply(Both_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Tuk_Mean_Silenced <- lapply(Tuk_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Bal_Mean_Silenced <- lapply(Bal_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

One_Mean_Silenced <- lapply(One_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Immune_Mean_Silenced <- lapply(Immune_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of silenced genes and XIST expression
Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Mean_Silenced[[i]], XIST_Ind_Counts[[i]])
}
names(Silenced_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

Tuk_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Tuk_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Tuk_Mean_Silenced[[i]], XIST_Ind_Counts[[i]])
}
names(Tuk_Silenced_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

Bal_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Bal_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Bal_Mean_Silenced[[i]], XIST_Ind_Counts[[i]])
}
names(Bal_Silenced_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

One_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  One_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(One_Mean_Silenced[[i]], XIST_Ind_Counts[[i]])
}
names(One_Silenced_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

Immune_Silenced_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Immune_Silenced_Mean_Vs_XIST[[i]] <- Combine_Vectors(Immune_Mean_Silenced[[i]], XIST_Ind_Counts[[i]])
}
names(Immune_Silenced_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

# Rename columns in each df
Silenced_Mean_Vs_XIST <- lapply(Silenced_Mean_Vs_XIST, setNames, c("Mean_Silenced", "XIST")) 

Tuk_Silenced_Mean_Vs_XIST <- lapply(Tuk_Silenced_Mean_Vs_XIST, setNames, c("Tuk_Mean_Silenced", "XIST")) 

Bal_Silenced_Mean_Vs_XIST <- lapply(Bal_Silenced_Mean_Vs_XIST, setNames, c("Bal_Mean_Silenced", "XIST")) 

One_Silenced_Mean_Vs_XIST <- lapply(One_Silenced_Mean_Vs_XIST, setNames, c("One_Mean_Silenced", "XIST")) 

Immune_Silenced_Mean_Vs_XIST <- lapply(Immune_Silenced_Mean_Vs_XIST, setNames, c("Immune_Mean_Silenced", "XIST"))

# Apply lm to each df in list
lm_Silenced_XIST <- lapply(Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Tuk_Silenced_XIST <- lapply(Tuk_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Bal_Silenced_XIST <- lapply(Bal_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Bal_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_One_Silenced_XIST <- lapply(One_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(One_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Immune_Silenced_XIST <- lapply(Immune_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Immune_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Res_Silenced_XIST <- lapply(lm_Silenced_XIST, R_Squared)

Tuk_Res_Silenced_XIST <- lapply(lm_Tuk_Silenced_XIST, R_Squared)

Bal_Res_Silenced_XIST <- lapply(lm_Bal_Silenced_XIST, R_Squared)

One_Res_Silenced_XIST <- lapply(lm_One_Silenced_XIST, R_Squared)

Immune_Res_Silenced_XIST <- lapply(lm_Immune_Silenced_XIST, R_Squared)

# Add vector as column to df
MeanX_XIST.df$Silenced_Mean <- unlist(Res_Silenced_XIST)

MeanX_XIST.df$Tuk_Silenced_Mean <- unlist(Tuk_Res_Silenced_XIST)

MeanX_XIST.df$Bal_Silenced_Mean <- unlist(Bal_Res_Silenced_XIST)

MeanX_XIST.df$One_Silenced_Mean <- unlist(One_Res_Silenced_XIST)

MeanX_XIST.df$Immune_Silenced_Mean <- unlist(Immune_Res_Silenced_XIST)

# ________________________________________________________________________________________________________
# Table 1; Columns 7-10
# Correlation of all genes reported as variably silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Variable_In_At_Least_One", "Variable_In_Tukiainen", "Variable_In_Balaton", "Immune_Gene_Variable_In_At_Least_One"

# Subset list of genes variably silenced from each df in list 
One_Variable <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_At_Least_One)
})

Tuk_Variable <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_Tukiainen)
})

Bal_Variable <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Variable_In_Balaton)
})

Immune_Variable <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Gene_Variable_In_At_Least_One)
})

# Get the mean values of putatively variably silenced X chm genes
One_Mean_Variable <- lapply(One_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Tuk_Mean_Variable <- lapply(Tuk_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Bal_Mean_Variable <- lapply(Bal_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Immune_Mean_Variable <- lapply(Immune_Variable, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of variably silenced genes and XIST expression
One_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  One_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(One_Mean_Variable[[i]], XIST_Ind_Counts[[i]])
}
names(One_Mean_Variable_Vs_XIST) <- names(XIST_Ind_Counts)

Tuk_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Tuk_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Tuk_Mean_Variable[[i]], XIST_Ind_Counts[[i]])
}
names(Tuk_Mean_Variable_Vs_XIST) <- names(XIST_Ind_Counts)

Bal_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Bal_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Bal_Mean_Variable[[i]], XIST_Ind_Counts[[i]])
}
names(Bal_Mean_Variable_Vs_XIST) <- names(XIST_Ind_Counts)

Immune_Mean_Variable_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Immune_Mean_Variable_Vs_XIST[[i]] <- Combine_Vectors(Immune_Mean_Variable[[i]], XIST_Ind_Counts[[i]])
}
names(Immune_Mean_Variable_Vs_XIST) <- names(XIST_Ind_Counts)

# Rename columns in each df
One_Mean_Variable_Vs_XIST <- lapply(One_Mean_Variable_Vs_XIST, setNames, c("One_Mean_Variable", "XIST")) 

Tuk_Mean_Variable_Vs_XIST <- lapply(Tuk_Mean_Variable_Vs_XIST, setNames, c("Tuk_Mean_Variable", "XIST")) 

Bal_Mean_Variable_Vs_XIST <- lapply(Bal_Mean_Variable_Vs_XIST, setNames, c("Bal_Mean_Variable", "XIST")) 

Immune_Mean_Variable_Vs_XIST <- lapply(Immune_Mean_Variable_Vs_XIST, setNames, c("Immune_Mean_Variable", "XIST")) 

# Apply lm to each df in list
lm_One_Variable_XIST <- lapply(One_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(One_Mean_Variable ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Tuk_Variable_XIST <- lapply(Tuk_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Variable ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Bal_Variable_XIST <- lapply(Bal_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Bal_Mean_Variable ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Immune_Variable_XIST <- lapply(Immune_Mean_Variable_Vs_XIST, function(x) {
  z <- lm(Immune_Mean_Variable ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Res_One_Variable_XIST <- lapply(lm_One_Variable_XIST, R_Squared)

Res_Tuk_Variable_XIST <- lapply(lm_Tuk_Variable_XIST, R_Squared)

Res_Bal_Variable_XIST <- lapply(lm_Bal_Variable_XIST, R_Squared)

Res_Immune_Variable_XIST <- lapply(lm_Immune_Variable_XIST, R_Squared)

# Add vector as column to df
MeanX_XIST.df$One_Variable_Mean <- unlist(Res_One_Variable_XIST)

MeanX_XIST.df$Tuk_Variable_Mean <- unlist(Res_Tuk_Variable_XIST)

MeanX_XIST.df$Bal_Variable_Mean <- unlist(Res_Bal_Variable_XIST)

MeanX_XIST.df$Immune_Variable_Mean <- unlist(Res_Immune_Variable_XIST)

# ________________________________________________________________________________________________________
# Table 1; Columns 11-14
# Correlation of all genes reported as incompletely silenced with XIST
# ________________________________________________________________________________________________________

# Categories:
# "Incomplete_In_At_Least_One", Incomplete_In_Tukiainen", "Incomplete_In_Balaton", "Immune_Genes_Incomplete_In_At_Least_One"

# Subset list of genes incompletly silenced from each df in list 
One_Incomplete <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_At_Least_One)
})

Tuk_Incomplete <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_Tukiainen)
})

Bal_Incomplete <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Incomplete_In_Balaton)
})

Immune_Incomplete <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Incomplete_In_At_Least_One)
})

# Get the mean values of putatively incompletely silenced X chm genes
One_Mean_Incomplete <- lapply(One_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Tuk_Mean_Incomplete <- lapply(Tuk_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Bal_Mean_Incomplete <- lapply(Bal_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Immune_Mean_Incomplete <- lapply(Immune_Incomplete, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of incompletely silenced genes and XIST expression
One_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  One_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(One_Mean_Incomplete[[i]], XIST_Ind_Counts[[i]])
}
names(One_Mean_Incomplete_Vs_XIST) <- names(XIST_Ind_Counts)

Tuk_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Tuk_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Tuk_Mean_Incomplete[[i]], XIST_Ind_Counts[[i]])
}
names(Tuk_Mean_Incomplete_Vs_XIST) <- names(XIST_Ind_Counts)

Bal_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Bal_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Bal_Mean_Incomplete[[i]], XIST_Ind_Counts[[i]])
}
names(Bal_Mean_Incomplete_Vs_XIST) <- names(XIST_Ind_Counts)

Immune_Mean_Incomplete_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Immune_Mean_Incomplete_Vs_XIST[[i]] <- Combine_Vectors(Immune_Mean_Incomplete[[i]], XIST_Ind_Counts[[i]])
}
names(Immune_Mean_Incomplete_Vs_XIST) <- names(XIST_Ind_Counts)

# Rename columns in each df
One_Mean_Incomplete_Vs_XIST <- lapply(One_Mean_Incomplete_Vs_XIST, setNames, c("One_Mean_Incomplete", "XIST")) 

Tuk_Mean_Incomplete_Vs_XIST <- lapply(Tuk_Mean_Incomplete_Vs_XIST, setNames, c("Tuk_Mean_Incomplete", "XIST")) 

Bal_Mean_Incomplete_Vs_XIST <- lapply(Bal_Mean_Incomplete_Vs_XIST, setNames, c("Bal_Mean_Incomplete", "XIST")) 

Immune_Mean_Incomplete_Vs_XIST <- lapply(Immune_Mean_Incomplete_Vs_XIST, setNames, c("Immune_Mean_Incomplete", "XIST")) 

# Apply lm to each df in list
lm_One_Incomplete_XIST <- lapply(One_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(One_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Tuk_Incomplete_XIST <- lapply(Tuk_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Bal_Incomplete_XIST <- lapply(Bal_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Bal_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Immune_Incomplete_XIST <- lapply(Immune_Mean_Incomplete_Vs_XIST, function(x) {
  z <- lm(Immune_Mean_Incomplete ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Res_One_Incomplete_XIST <- lapply(lm_One_Incomplete_XIST, R_Squared)

Res_Tuk_Incomplete_XIST <- lapply(lm_Tuk_Incomplete_XIST, R_Squared)

Res_Bal_Incomplete_XIST <- lapply(lm_Bal_Incomplete_XIST, R_Squared)

Res_Immune_Incomplete_XIST <- lapply(lm_Immune_Incomplete_XIST, R_Squared)

# Add vector as column to df
MeanX_XIST.df$One_Incomplete_Mean <- unlist(Res_One_Incomplete_XIST)

MeanX_XIST.df$Tuk_Incomplete_Mean <- unlist(Res_Tuk_Incomplete_XIST)

MeanX_XIST.df$Bal_Incomplete_Mean <- unlist(Res_Bal_Incomplete_XIST)

MeanX_XIST.df$Immune_Incomplete_Mean <- unlist(Res_Immune_Incomplete_XIST)

# ________________________________________________________________________________________________________
# Table 1; Columns 15-18
# Correlation of all genes/ all genes not evaluated / PAR genes with XIST
# ________________________________________________________________________________________________________

# Categories:
# "All_Evaluated_Balaton_Tukiainen", "Not_Evaluated_In_Either", "Immune_Genes_Not_Evaluated", "PAR_In_Balaton"   

# Subset list of genes from each df in list 
All_Eval <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$All_Evaluated_Balaton_Tukiainen)
})

Not_Eval <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Not_Evaluated_In_Either)
})

Immune_Not_Eval <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Immune_Genes_Not_Evaluated)
})

PAR_Bal <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$PAR_In_Balaton)
})

# Get the mean values X chm genes
All_Eval_Mean <- lapply(All_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Not_Eval_Mean <- lapply(Not_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Immune_Not_Eval_Mean <- lapply(Immune_Not_Eval, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

PAR_Mean <- lapply(PAR_Bal, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

# Combine list of mean value of genes and XIST expression
All_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  All_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(All_Eval_Mean[[i]], XIST_Ind_Counts[[i]])
}
names(All_Eval_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

Not_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Not_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Not_Eval_Mean[[i]], XIST_Ind_Counts[[i]])
}
names(Not_Eval_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

Immune_Not_Eval_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  Immune_Not_Eval_Mean_Vs_XIST[[i]] <- Combine_Vectors(Immune_Not_Eval_Mean[[i]], XIST_Ind_Counts[[i]])
}
names(Immune_Not_Eval_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

PAR_Mean_Vs_XIST <- list()
for (i in 1:length(XIST_Ind_Counts)){
  PAR_Mean_Vs_XIST[[i]] <- Combine_Vectors(PAR_Mean[[i]], XIST_Ind_Counts[[i]])
}
names(PAR_Mean_Vs_XIST) <- names(XIST_Ind_Counts)

# Rename columns in each df
All_Eval_Mean_Vs_XIST <- lapply(All_Eval_Mean_Vs_XIST, setNames, c("All_Eval_Mean", "XIST")) 

Not_Eval_Mean_Vs_XIST <- lapply(Not_Eval_Mean_Vs_XIST, setNames, c("Not_Eval_Mean", "XIST")) 

Immune_Not_Eval_Mean_Vs_XIST <- lapply(Immune_Not_Eval_Mean_Vs_XIST, setNames, c("Immune_Not_Eval_Mean", "XIST")) 

PAR_Mean_Vs_XIST <- lapply(PAR_Mean_Vs_XIST, setNames, c("PAR_Mean", "XIST")) 

# Apply lm to each df in list
lm_All_Eval <- lapply(All_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(All_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Not_Eval <- lapply(Not_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(Not_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Immune_Not_Eval <- lapply(Immune_Not_Eval_Mean_Vs_XIST, function(x) {
  z <- lm(Immune_Not_Eval_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_PAR <- lapply(PAR_Mean_Vs_XIST, function(x) {
  z <- lm(PAR_Mean ~ XIST, data = x, na.action = na.omit)
  return(z)
})

# Apply function to list of dfs
Res_All_Eval <- lapply(lm_All_Eval, R_Squared)

Res_Not_Eval <- lapply(lm_Not_Eval, R_Squared)

Res_Immune_Not_Eval <- lapply(lm_Immune_Not_Eval, R_Squared)

Res_PAR <- lapply(lm_PAR, R_Squared)

# Add vector as column to df
MeanX_XIST.df$All_Eval <- unlist(Res_All_Eval)

MeanX_XIST.df$Not_Eval <- unlist(Res_Not_Eval)

MeanX_XIST.df$Immune_Not_Eval <- unlist(Res_Immune_Not_Eval)

MeanX_XIST.df$PAR <- unlist(Res_PAR)

# ________________________________________________________________________________________________________
#  Write table 1 and summary 
# ________________________________________________________________________________________________________

# Average R^2 of silenced genes reported in both studies for females vs males
df1 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Silenced = mean(Silenced_Mean))

df2 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Tuk_Silenced = mean(Tuk_Silenced_Mean))

df3 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Bal_Silenced = mean(Bal_Silenced_Mean))

df4 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(One_Silenced = mean(One_Silenced_Mean))

df5 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Immune_Silenced = mean(Immune_Silenced_Mean))

# Average R^2 of silenced genes reported in both studies for females vs males
df6 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(One_Variable = mean(One_Variable_Mean))

df7 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Tuk_Variable = mean(Tuk_Variable_Mean))

df8 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Bal_Variable = mean(Bal_Variable_Mean))

df9 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Immune_Variable = mean(Immune_Variable_Mean))

# Average R^2 of silenced genes reported in both studies for females vs males
df10 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(One_Incomplete = mean(One_Incomplete_Mean))

df11 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Tuk_Incomplete = mean(Tuk_Incomplete_Mean))

df12 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Bal_Incomplete = mean(Bal_Incomplete_Mean))

df13 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Immune_Incomplete = mean(Immune_Incomplete_Mean))

# Average R^2 of genes for females vs males
df14 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(All_Eval = mean(All_Eval))

df15 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Not_Eval = mean(Not_Eval))

df16 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Immune_Not_Eval = mean(Immune_Not_Eval))

df17 <- MeanX_XIST.df %>% group_by(Sex) %>% summarise(Bal_PAR = mean(PAR))

# Combine dfs of averages
df_lst <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17)

Combined <- Reduce(merge, lapply(df_lst, function(x) data.frame(x, rn = row.names(x))))

# Add column listing the total number of tissues for each individual
# Get number of columns for each person in Ind_Counts excluding the first 2
Num_Tissues <- sapply(Ind_Counts, ncol)
Num_Tissues <- Num_Tissues - 2

MeanX_XIST.df$Number_Tissues <- Num_Tissues

# Write to file
write.csv(Combined, "Linear_Models_Summary.csv")

write.csv(MeanX_XIST.df, "Individual_Correlation.csv")

# ________________________________________________________________________________________________________
#  Plots
# ________________________________________________________________________________________________________

# Get individual with highest average X chromsomes expression R^2 and plot
# ASK: Do I just remove the people with only a handful of samples? What would be a reasonable cut-off?
# For now: Cutoff of min 5 tissue samples
Sub_MeanX_XIST.df <- MeanX_XIST.df[which(MeanX_XIST.df$Number_Tissues > 5),]
Highest_R2 <- Sub_MeanX_XIST.df[which.max(Sub_MeanX_XIST.df$R2_MeanX),]

# Plot 
plot(XIST_Ind_Counts$'GTEX-R55G', MeanX_Ind_Counts$'GTEX-R55G', 
     main = "Scatterplot of individual with highest R^2",
     xlab = 'XIST count',
     ylab = 'Mean X chromosome gene count')


# Make a histogram that shows the female and male R^2 values
Fem_R2 <- subset(MeanX_XIST.df, Sex == "Female")$R2_Mean
Male_R2 <- subset(MeanX_XIST.df, Sex == "Male")$R2_Mean


hist(Fem_R2)
hist(Male_R2)
