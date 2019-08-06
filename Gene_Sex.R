# Perform linear regression on the gene count data (version 7).
#  Find correlation between mean X chm expression and XIST across tissues per individual

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
# Table 1; Column 2-6
# Correlation of all genes reported as silenced in both studies
# ________________________________________________________________________________________________________

# Categories:
# "Immune_Genes_Silenced_In_At_Least_One", "Silenced_In_Tukiainen", "Silenced_In_Balaton", "Silenced_In_Both", "Silenced_In_At_Least_One"  

# Subset list of genes silenced in both studies from each df in list 
Both_Silenced <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Both)
})

Tuk_Silenced <- lapply(Ind_Counts, function(x) {
  filter(x, gene_name %in% Gene_Lst$Silenced_In_Tukiainen)
})

# Get the mean values of putatively silenced X chm genes
Mean_Silenced <- lapply(Both_Silenced, function(x){
  colMeans(x[sapply(x, is.numeric)]) 
})

Tuk_Mean_Silenced <- lapply(Tuk_Silenced, function(x){
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

# Rename columns in each df
Silenced_Mean_Vs_XIST <- lapply(Silenced_Mean_Vs_XIST, setNames, c("Mean_Silenced", "XIST")) 

Tuk_Silenced_Mean_Vs_XIST <- lapply(Tuk_Silenced_Mean_Vs_XIST, setNames, c("Tuk_Mean_Silenced", "XIST")) 

# Apply lm to each df in list
lm_Silenced_XIST <- lapply(Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})

lm_Tuk_Silenced_XIST <- lapply(Tuk_Silenced_Mean_Vs_XIST, function(x) {
  z <- lm(Tuk_Mean_Silenced ~ XIST, data = x, na.action = na.omit)
  return(z)
})


# Apply function to list of dfs
Res_Silenced_XIST <- lapply(lm_Silenced_XIST, R_Squared)

Tuk_Res_Silenced_XIST <- lapply(lm_Tuk_Silenced_XIST, R_Squared)

# Add vector as column to df
MeanX_XIST.df$Silenced_Mean <- unlist(Res_Silenced_XIST)

MeanX_XIST.df$Tuk_Silenced_Mean <- unlist(Tuk_Res_Silenced_XIST)

# Average R^2 of silenced genes reported in both studies for females vs males
MeanX_XIST.df %>% group_by(Sex) %>% summarise(mean = mean(Silenced_Mean))

MeanX_XIST.df %>% group_by(Sex) %>% summarise(mean = mean(Tuk_Silenced_Mean))

