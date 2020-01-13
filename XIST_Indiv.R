# Perform linear regression on the gene count data (version 7).
# Fit linear models by individual: 
# Reponse variable: Mean X chromosome expression and various sets of X-linked genes 
# Predictor variable: XIST
setwd("~/XIST")

# Constants
COUNTS <- "~/XIST/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct" # TPM normalized
METRICS <- "~/XIST/Files/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"
PHENOTYPES <- "~/XIST/Files/GTEX_v7_Annotations_SubjectPhenotypesDS.txt"
GENCODE <- "gencode.v19.genes.v7.patched_contigs.gff3" # has to be in same dir
GENE_LST <- "~/XIST/Files/X_Genes_Status.json"

# Results
TABLE <- "~/XIST/Indiv/Individual_Correlations.csv"

# Data 
DATA <- "XIST_Indiv_121719.RData"

# Load libraries
library(readr)
library(refGenome) # Parse .gff file
library(dplyr)
library(data.table)
library(stringr)
library(broom)
library(rjson)

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
# Don't want to include XIST in mean X chromsome count values
X_Genes <- X_Genes[!X_Genes$gene_name == 'XIST',]

# ______________________________________________________________________________________________________________________
# Tissue summary
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

# ______________________________________________________________________________________________________________________
# X chromosome counts
# ______________________________________________________________________________________________________________________
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

# Check for any empty data frames/have at least 3 samples
# Min to do regression
# i.e. > 5 cols including 'gene_id' and 'gene_name'
length(Ind_Counts) == length(Filter(function(x) dim(x)[2] > 5, Ind_Counts)) # FALSE

# Remove individuals who have at least 3 samples
Ind_Counts <- Ind_Counts[sapply(Ind_Counts, function(x) dim(x)[2]) > 5]

# Check for any empty data frames
length(Ind_Counts) == length(Filter(function(x) dim(x)[2] > 5, Ind_Counts)) # TRUE

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
Ind_XCounts <- lapply(Ind_Counts, function(x) {
  filter(x, gene_id %in% X_Genes$gene_id)
})

# Check for missing values
any(sapply(Ind_XCounts, function(x) sum(is.na(x)))) # FALSE

# Get the mean values of the X chm genes
Mean_Ind_XCounts <- lapply(Ind_XCounts, function(x){
  colMeans(x[sapply(x, is.numeric)])
})

# Get XIST values from each data frame
XIST_Ind_Counts <- lapply(Ind_Counts, function(x) {
  as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
})

# ______________________________________________________________________________________________________________________
# Table 1; Column 1:
# Correlate expression from XIST with X chromosome expression within a person across all of their tissues.
# ______________________________________________________________________________________________________________________
# Check that all the samples are present in both lists of dfs
length(XIST_Ind_Counts) == length(Mean_Ind_XCounts) # TRUE

# How many individuals are there?
length(XIST_Ind_Counts) # 680

# Are all the same individuals listed in both lists?
all(names(XIST_Ind_Counts) == names(Mean_Ind_XCounts)) # TRUE

# Are they listed in the same order?
identical(names(XIST_Ind_Counts), names(Mean_Ind_XCounts)) # TRUE

# Function to combine vectors into df
Combine_Vectors <- function(a, b){
  data.frame(a, b)
}

# Combine list of mean value of silenced genes and XIST expression
Combine_Lsts <- function(x, y, z){
  x <- list()
  for (i in 1:length(z)){
    x[[i]] <- Combine_Vectors(y[[i]], z[[i]])
  }
  names(x) <- names(y)
  return(x)
}
MeanX_Vs_XIST <- Combine_Lsts(x='MeanX_Vs_XIST', y=Mean_Ind_XCounts, z=XIST_Ind_Counts)

# Rename cols
Rename_Col <- function(x, a, b){
  x <- setNames(x, c(a, b))
  return(x)
}
MeanX_Vs_XIST <- Map(Rename_Col, x=MeanX_Vs_XIST, a='MeanX', b='XIST')

# Apply lm to each df in list
Linear_Model.1 <- function(x) {
  z <- lm(MeanX ~ XIST, data = x)
  return(z)
}
lm.MeanX_XIST <- lapply(MeanX_Vs_XIST, Linear_Model.1)

# Function to extract r squared values
Regression_Res <- function(lm){ # expecting object of class 'lm'
  sum. <- summary(lm)
  r.2 <- sum.$r.squared
  p. <- summary(lm)$fstatistic
  p.val <- pf(p.[1], p.[2], p.[3], lower.tail=FALSE, log.p=FALSE) #pf: F distribution ; arg 2 and 3 are the deg. freedom
  attributes(p.val) <- NULL
  res <- data.frame(p_val=p.val, r_2=r.2)
  return(res)
}

# Apply function to list of dfs
Res.MeanX_XIST <- lapply(lm.MeanX_XIST, Regression_Res)

# Make table summarizing results
Regression <- as.data.frame(do.call(rbind, Res.MeanX_XIST))

# Label columns
colnames(Regression) <- c("pval_MeanX", "R2_MeanX")

# Add column with sex
# Get list of female IDs
Female_IDs <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

Sample_Sex <- sapply(rownames(Regression), function(x) {
  if (x %in% Female_IDs > 0) {"Female"}
  else {"Male"}
}, simplify = TRUE)
Regression <- cbind(Sample_Sex, Regression)

# Average R^2 for females vs males
Regression %>% group_by(Sample_Sex) %>% summarise(mean = mean(R2_MeanX)) # female: 0.36, male:0.0662

# ______________________________________________________________________________________________________________________
# Table 1; Columns 2-12
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
Both_Silenced <- Map(Filter_Func, 
                     x=Ind_XCounts, 
                     y='Silenced_In_Both')
Tuk_Silenced <- Map(Filter_Func, 
                    x=Ind_XCounts,
                    y='Silenced_In_Tukianen')
Bal_Silenced <- Map(Filter_Func, 
                    x=Ind_XCounts, 
                    y='Silenced_In_Balaton')
One_Silenced <- Map(Filter_Func, 
                    x=Ind_XCounts, 
                    y='Silenced_In_At_Least_One')
Immune_Silenced <- Map(Filter_Func, 
                       x=Ind_XCounts, 
                       y='Immune_Genes_Silenced_In_At_Least_One')

# Get the mean values of putatively silenced X chm genes
Mean_Val <- function(x){
  res <- colMeans(x[sapply(x, is.numeric)]) 
  return(res)
}
Mean_Silenced <- lapply(Both_Silenced, Mean_Val)
Tuk_Mean_Silenced <- lapply(Both_Silenced, Mean_Val)
Bal_Mean_Silenced <- lapply(Bal_Silenced, Mean_Val)
One_Mean_Silenced <- lapply(One_Silenced, Mean_Val)
Immune_Mean_Silenced <- lapply(Immune_Silenced, Mean_Val)

# Combine list of mean value of silenced genes and XIST expression
Silenced_Mean_XIST <- Combine_Lsts(x='Silenced_Mean_XIST', 
                                   y=Mean_Silenced, 
                                   z=XIST_Ind_Counts)
Tuk_Silenced_Mean_XIST <- Combine_Lsts(x='Tuk_Silenced_Mean_XIST', 
                                       y=Tuk_Mean_Silenced, 
                                       z=XIST_Ind_Counts)
Bal_Silenced_Mean_XIST <- Combine_Lsts(x='Bal_Silenced_Mean_XIST', 
                                       y=Bal_Mean_Silenced, 
                                       z=XIST_Ind_Counts)
One_Silenced_Mean_XIST <- Combine_Lsts(x='One_Silenced_Mean_XIST', 
                                       y=One_Mean_Silenced, 
                                       z=XIST_Ind_Counts)
Immune_Silenced_Mean_XIST <- Combine_Lsts(x='Immune_Silenced_Mean_XIST', 
                                          y=Immune_Mean_Silenced, 
                                          z=XIST_Ind_Counts)

# Rename columns in each df
Silenced_Mean_XIST <- Map(Rename_Col, 
                          x=Silenced_Mean_XIST, 
                          a='Mean_Silenced', 
                          b='XIST')
Tuk_Silenced_Mean_XIST <- Map(Rename_Col, 
                              x=Tuk_Silenced_Mean_XIST, 
                              a='Mean_Silenced', 
                              b='XIST')
Bal_Silenced_Mean_XIST <- Map(Rename_Col, 
                              x=Bal_Silenced_Mean_XIST, 
                              a='Mean_Silenced', 
                              b='XIST')
One_Silenced_Mean_XIST <- Map(Rename_Col, 
                              x=One_Silenced_Mean_XIST, 
                              a='Mean_Silenced', 
                              b='XIST')
Immune_Silenced_Mean_XIST <- Map(Rename_Col, 
                                 x=Immune_Silenced_Mean_XIST, 
                                 a='Mean_Silenced', 
                                 b='XIST')

# Apply lm to each df in list
Linear_Model.2 <- function(x) {
  z <- lm(Mean_Silenced ~ XIST, data = x)
  return(z)
}
lm.Silenced_XIST <- lapply(Silenced_Mean_XIST, Linear_Model.2)
lm.Tuk_Silenced_XIST <- lapply(Tuk_Silenced_Mean_XIST, Linear_Model.2)
lm.Bal_Silenced_XIST <- lapply(Bal_Silenced_Mean_XIST, Linear_Model.2)
lm.One_Silenced_XIST <- lapply(One_Silenced_Mean_XIST, Linear_Model.2)
lm.Immune_Silenced_XIST <- lapply(Immune_Silenced_Mean_XIST, Linear_Model.2)

# Apply function to list of dfs
Res.Silenced_XIST <- lapply(lm.Silenced_XIST, Regression_Res)
Tuk.Res_Silenced_XIST <- lapply(lm.Tuk_Silenced_XIST, Regression_Res)
Bal.Res_Silenced_XIST <- lapply(lm.Bal_Silenced_XIST, Regression_Res)
One.Res_Silenced_XIST <- lapply(lm.One_Silenced_XIST, Regression_Res)
Immune.Res_Silenced_XIST <- lapply(lm.Immune_Silenced_XIST, Regression_Res)

# Add vector as column to df
Return_R2 <- function(x){
  res <- unlist(x$r_2)
  return(res)
}

Return_pval <- function(x){
  res <- unlist(x$p_val)
  return(res)
}

Regression$pval_Silenced_Mean <- lapply(Res.Silenced_XIST, Return_pval)
Regression$R2_Silenced_Mean <- lapply(Res.Silenced_XIST, Return_R2)

Regression$pval_Tuk_Silenced_Mean <- lapply(Tuk.Res_Silenced_XIST, Return_pval)
Regression$R2_Tuk_Silenced_Mean <- lapply(Tuk.Res_Silenced_XIST, Return_R2)

Regression$pval_Bal_Silenced_Mean <- lapply(Bal.Res_Silenced_XIST, Return_pval)
Regression$R2_Bal_Silenced_Mean <- lapply(Bal.Res_Silenced_XIST, Return_R2)

Regression$pval_One_Silenced_Mean <- lapply(One.Res_Silenced_XIST, Return_pval)
Regression$R2_One_Silenced_Mean <- lapply(One.Res_Silenced_XIST, Return_R2)

Regression$pval_Immune_Silenced_Mean <- lapply(Immune.Res_Silenced_XIST, Return_pval)
Regression$R2_Immune_Silenced_Mean <- lapply(Immune.Res_Silenced_XIST, Return_R2)

# ______________________________________________________________________________________________________________________
# Table 1; Columns 13-21
# Correlation of all genes reported as variably silenced with XIST
# ______________________________________________________________________________________________________________________
# Categories:
# "Variable_In_At_Least_One", "Variable_In_Tukiainen", "Variable_In_Balaton", "Immune_Gene_Variable_In_At_Least_One"

# Subset list of genes variably silenced from each df of x counts in list 
One.Variable <- Map(Filter_Func, 
                    x=Ind_XCounts, 
                    y='Variable_In_At_Least_One')
Tuk.Variable <- Map(Filter_Func, 
                    x=Ind_XCounts, 
                    y='Variable_In_Tukiainen')
Bal.Variable <- Map(Filter_Func, 
                    x=Ind_XCounts, 
                    y='Variable_In_Balaton')
Immune.Variable <- Map(Filter_Func, 
                       x=Ind_XCounts, 
                       y='Immune_Gene_Variable_In_At_Least_One')

# Get the mean values of putatively variably silenced X chm genes
One.Mean_Variable <- lapply(One.Variable, Mean_Val)
Tuk.Mean_Variable <- lapply(Tuk.Variable, Mean_Val)
Bal.Mean_Variable <- lapply(Bal.Variable, Mean_Val)
Immune.Mean_Variable <- lapply(Immune.Variable, Mean_Val)

# Combine list of mean value of variably silenced genes and XIST expression
One.Mean_Variable_Vs_XIST <- Combine_Lsts(x='One_Mean_Variable_XIST', 
                                          y=One.Mean_Variable, 
                                          z=XIST_Ind_Counts)
Tuk.Mean_Variable_Vs_XIST <- Combine_Lsts(x='Tuk_Mean_Variable_XIST', 
                                          y=Tuk.Mean_Variable, 
                                          z=XIST_Ind_Counts)
Bal.Mean_Variable_Vs_XIST <- Combine_Lsts(x='Bal_Mean_Variable_XIST', 
                                          y=Bal.Mean_Variable, 
                                          z=XIST_Ind_Counts)
Immune.Mean_Variable_Vs_XIST <- Combine_Lsts(x='Immune_Mean_Variable_XIST', 
                                             y=Immune.Mean_Variable, 
                                             z=XIST_Ind_Counts)

# Rename columns in each df
One.Mean_Variable_Vs_XIST <- Map(Rename_Col, 
                                 x=One.Mean_Variable_Vs_XIST, 
                                 a='Mean_Variable', 
                                 b='XIST')
Tuk.Mean_Variable_Vs_XIST <- Map(Rename_Col, 
                                 x=Tuk.Mean_Variable_Vs_XIST, 
                                 a='Mean_Variable', 
                                 b='XIST')
Bal.Mean_Variable_Vs_XIST <- Map(Rename_Col, 
                                 x=Bal.Mean_Variable_Vs_XIST, 
                                 a='Mean_Variable',
                                 b='XIST')
Immune.Mean_Variable_Vs_XIST <- Map(Rename_Col, 
                                    x=Immune.Mean_Variable_Vs_XIST, 
                                    a='Mean_Variable', 
                                    b='XIST') 

# Apply lm to each df in list
Linear_Model.3 <- function(x) {
  z <- lm(Mean_Variable ~ XIST, data = x)
  return(z)
}

lm.One_Variable_XIST <- lapply(One.Mean_Variable_Vs_XIST, Linear_Model.3)
lm.Tuk_Variable_XIST <- lapply(Tuk.Mean_Variable_Vs_XIST, Linear_Model.3)  
lm.Bal_Variable_XIST <- lapply(Bal.Mean_Variable_Vs_XIST, Linear_Model.3)  
lm.Immune_Variable_XIST <- lapply(Immune.Mean_Variable_Vs_XIST, Linear_Model.3)  

# Apply function to list of dfs
Res.One_Variable_XIST <- lapply(lm.One_Variable_XIST, Regression_Res)
Res.Tuk_Variable_XIST <- lapply(lm.Tuk_Variable_XIST, Regression_Res)
Res.Bal_Variable_XIST <- lapply(lm.Bal_Variable_XIST, Regression_Res)
Res.Immune_Variable_XIST <- lapply(lm.Immune_Variable_XIST, Regression_Res)

# Add vector as column to df
Regression$pval_One_Variable_Mean <- lapply(Res.Silenced_XIST, Return_pval)
Regression$R2_One_Variable_Mean <- lapply(Res.Silenced_XIST, Return_R2)

Regression$pval_Tuk_Variable_Mean <- lapply(Res.Tuk_Variable_XIST, Return_pval)
Regression$R2_Tuk_Variable_Mean <- lapply(Res.Tuk_Variable_XIST, Return_R2)

Regression$pval_Bal_Variable_Mean <- lapply(Res.Bal_Variable_XIST, Return_pval)
Regression$R2_Bal_Variable_Mean <- lapply(Res.Bal_Variable_XIST, Return_R2)

Regression$pval_Immune_Variable_Mean <- lapply(Res.Immune_Variable_XIST, Return_pval)
Regression$R2_Immune_Variable_Mean <- lapply(Res.Immune_Variable_XIST, Return_R2)

# ______________________________________________________________________________________________________________________
# Table 1; Columns 21:29
# Correlation of all genes reported as incompletely silenced with XIST
# ______________________________________________________________________________________________________________________
# Categories:
# "Incomplete_In_At_Least_One", Incomplete_In_Tukiainen", "Incomplete_In_Balaton", 
# "Immune_Genes_Incomplete_In_At_Least_One"

# Subset list of genes incompletly silenced from each df of x counts in list 
One.Incomplete <- Map(Filter_Func, 
                      x=Ind_XCounts, 
                      y='Incomplete_In_At_Least_One')
Tuk.Incomplete <- Map(Filter_Func, 
                      x=Ind_XCounts, 
                      y='Incomplete_In_Tukiainen')
Bal.Incomplete <- Map(Filter_Func, 
                      x=Ind_XCounts, 
                      y='Incomplete_In_Balaton')
Immune.Incomplete <- Map(Filter_Func, 
                         x=Ind_XCounts, 
                         y='Immune_Genes_Incomplete_In_At_Least_One')

# Get the mean values of putatively incompletely silenced X chm genes
One.Mean_Incomplete <- lapply(One.Incomplete, Mean_Val)
Tuk.Mean_Incomplete <- lapply(Tuk.Incomplete, Mean_Val)
Bal.Mean_Incomplete <- lapply(Bal.Incomplete, Mean_Val)
Immune.Mean_Incomplete <- lapply(Immune.Incomplete, Mean_Val)

# Combine list of mean value of incompletely silenced genes and XIST expression
One.Mean_Incomplete_Vs_XIST <- Combine_Lsts(x='f.One.Mean_Incomplete_Vs_XIST',
                                            y=One.Mean_Incomplete, 
                                            z=XIST_Ind_Counts)
Tuk.Mean_Incomplete_Vs_XIST <- Combine_Lsts(x='f.Tuk.Mean_Incomplete_Vs_XIST', 
                                            y=Tuk.Mean_Incomplete, 
                                            z=XIST_Ind_Counts)
Bal.Mean_Incomplete_Vs_XIST <- Combine_Lsts(x='f.Bal.Mean_Incomplete_Vs_XIST', 
                                            y=Bal.Mean_Incomplete, 
                                            z=XIST_Ind_Counts)
Immune.Mean_Incomplete_Vs_XIST <- Combine_Lsts(x='f.Immune.Mean_Incomplete_Vs_XIST', 
                                               y=Immune.Mean_Incomplete, 
                                               z=XIST_Ind_Counts)

# Rename columns in each df
One.Mean_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                   x=One.Mean_Incomplete_Vs_XIST, 
                                   a='Mean_Incomplete', 
                                   b='XIST') 
Tuk.Mean_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                   x=Tuk.Mean_Incomplete_Vs_XIST, 
                                   a='Mean_Incomplete',
                                   b='XIST') 
Bal.Mean_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                   x=Bal.Mean_Incomplete_Vs_XIST, 
                                   a='Mean_Incomplete', 
                                   b='XIST') 
Immune.Mean_Incomplete_Vs_XIST <- Map(Rename_Col, 
                                      x=Immune.Mean_Incomplete_Vs_XIST, 
                                      a='Mean_Incomplete', 
                                      b='XIST')  

# Apply lm to each df in list
Linear_Model.4 <- function(x) {
  z <- lm(Mean_Incomplete ~ XIST, data = x)
  return(z)
}

lm.One_Incomplete_XIST <- lapply(One.Mean_Incomplete_Vs_XIST, Linear_Model.4)
lm.Tuk_Incomplete_XIST <- lapply(Tuk.Mean_Incomplete_Vs_XIST, Linear_Model.4)
lm.Bal_Incomplete_XIST <- lapply(Bal.Mean_Incomplete_Vs_XIST, Linear_Model.4)
lm.Immune_Incomplete_XIST <- lapply(Immune.Mean_Incomplete_Vs_XIST, Linear_Model.4)

# Apply function to list of dfs
Res.One_Incomplete_XIST <- lapply(lm.One_Incomplete_XIST, Regression_Res)
Res.Tuk_Incomplete_XIST <- lapply(lm.Tuk_Incomplete_XIST, Regression_Res)
Res.Bal_Incomplete_XIST <- lapply(lm.Bal_Incomplete_XIST, Regression_Res)
Res.Immune_Incomplete_XIST <- lapply(lm.Immune_Incomplete_XIST, Regression_Res)

# Add vector as column to df
Regression$pval_One_Incomplete_Mean <- lapply(Res.One_Incomplete_XIST, Return_pval) 
Regression$R2_One_Incomplete_Mean <- lapply(Res.One_Incomplete_XIST, Return_R2) 

Regression$pval_Tuk_Incomplete_Mean <- lapply(Res.Tuk_Incomplete_XIST, Return_pval)
Regression$R2_Tuk_Incomplete_Mean <- lapply(Res.Tuk_Incomplete_XIST, Return_R2)

Regression$pval_Bal_Incomplete_Mean <- lapply(Res.Bal_Incomplete_XIST, Return_pval)
Regression$R2_Bal_Incomplete_Mean <- lapply(Res.Bal_Incomplete_XIST, Return_R2)

Regression$pval_Immune_Incomplete_Mean <- lapply(Res.Immune_Incomplete_XIST, Return_pval)
Regression$R2_Immune_Incomplete_Mean <- lapply(Res.Immune_Incomplete_XIST, Return_R2)

# ______________________________________________________________________________________________________________________
# Table 1; Columns 29:37
# Correlation of all genes/ all genes not evaluated / PAR genes with XIST
# ______________________________________________________________________________________________________________________
# Categories:
# "All_Evaluated_Balaton_Tukiainen", "Not_Evaluated_In_Either", "Immune_Genes_Not_Evaluated", "PAR_In_Balaton"  
# Subset list of genes from each df of x counts in list 
All_Eval <- Map(Filter_Func, 
                x=Ind_XCounts, 
                y='All_Evaluated_Balaton_Tukiainen')
Not_Eval <- Map(Filter_Func,
                x=Ind_XCounts, 
                y='Not_Evaluated_In_Either')
Immune_Not_Eval <- Map(Filter_Func, 
                       x=Ind_XCounts, 
                       y='Immune_Genes_Not_Evaluated')
PAR_Bal <-  Map(Filter_Func, 
                x=Ind_XCounts, 
                y='PAR_In_Balaton')

# Get the mean values X chm genes
All_Eval.Mean <- lapply(All_Eval, Mean_Val)
Not_Eval.Mean <- lapply(Not_Eval, Mean_Val)
Immune_Not_Eval.Mean <- lapply(Immune_Not_Eval, Mean_Val)
PAR.Mean <- lapply(PAR_Bal, Mean_Val)

# Combine list of mean value of genes and XIST expression
All_Eval.Mean_Vs_XIST <- Combine_Lsts(x='All_Eval.Mean_Vs_XIST', 
                                      y=All_Eval.Mean, 
                                      z=XIST_Ind_Counts)
Not_Eval.Mean_Vs_XIST <- Combine_Lsts(x='Not_Eval.Mean_Vs_XIST', 
                                      y=Not_Eval.Mean, 
                                      z=XIST_Ind_Counts)
Immune_Not_Eval.Mean_Vs_XIST <- Combine_Lsts(x='Immune_Not_Eval.Mean_Vs_XIST', 
                                             y=Immune_Not_Eval.Mean, 
                                             z=XIST_Ind_Counts)
PAR.Mean_Vs_XIST <- Combine_Lsts(x='PAR.Mean_Vs_XIST', 
                                 y=PAR.Mean, 
                                 z=XIST_Ind_Counts)

# Rename columns in each df
All_Eval.Mean_Vs_XIST <- Map(Rename_Col, 
                             x=All_Eval.Mean_Vs_XIST,
                             a='Misc', 
                             b='XIST') 
Not_Eval.Mean_Vs_XIST <- Map(Rename_Col, 
                             x=Not_Eval.Mean_Vs_XIST, 
                             a='Misc',
                             b='XIST') 
Immune_Not_Eval.Mean_Vs_XIST <- Map(Rename_Col, 
                                    x=Immune_Not_Eval.Mean_Vs_XIST, 
                                    a='Misc', 
                                    b='XIST') 
PAR.Mean_Vs_XIST <- Map(Rename_Col,
                        x=PAR.Mean_Vs_XIST, 
                        a='Misc', 
                        b='XIST') 

# Apply lm to each df in list
Linear_Model.5 <- function(x) {
  z <- lm(Misc ~ XIST, data = x)
  return(z)
}

lm.All_Eval <- lapply(All_Eval.Mean_Vs_XIST, Linear_Model.5)
lm.Not_Eval <- lapply(Not_Eval.Mean_Vs_XIST, Linear_Model.5)
lm.Immune_Not_Eval <- lapply(Immune_Not_Eval.Mean_Vs_XIST, Linear_Model.5)
lm.PAR <- lapply(PAR.Mean_Vs_XIST, Linear_Model.5)

# Apply function to list of dfs
Res.All_Eval <- lapply(lm.All_Eval, Regression_Res)
Res.Not_Eval <- lapply(lm.Not_Eval, Regression_Res)
Res.Immune_Not_Eval <- lapply(lm.Immune_Not_Eval, Regression_Res)
Res.PAR <- lapply(lm.PAR, Regression_Res)

# Add vector as column to df
Regression$pval_All_Eval <- lapply(Res.All_Eval, Return_pval)
Regression$R2_All_Eval <- lapply(Res.All_Eval, Return_R2)

Regression$pval_Not_Eval <- lapply(Res.Not_Eval, Return_pval)
Regression$R2_Not_Eval <- lapply(Res.Not_Eval, Return_R2)

Regression$pval_Immune_Not_Eval <- lapply(Res.Immune_Not_Eval, Return_pval)
Regression$R2_Immune_Not_Eval <- lapply(Res.Immune_Not_Eval, Return_R2)

Regression$pval_PAR <- lapply(Res.PAR, Return_pval)
Regression$R2_PAR <- lapply(Res.PAR, Return_R2)

# ______________________________________________________________________________________________________________________
#  Write table 1 and summary
# ______________________________________________________________________________________________________________________
# Add col with mean(XIST) and sd(XIST)
Regression$Mean_XIST <- lapply(XIST_Ind_Counts, mean)
Regression$sd_XIST <- lapply(XIST_Ind_Counts, sd)

# Add column with number of samples per individual
Count_Rows <- function(x){
  res <- ncol(x)
  return(res)
}
Num_Samples <- as.character(lapply(Ind_Counts, Count_Rows))

# Add column with number of samples per individual
Regression <- cbind(Num_Samples=Num_Samples, Regression) 

# Convert Num_Tissues from factor to numeric
Regression$Num_Samples <- as.numeric(as.character(Regression$Num_Samples))

# Convert rest of cols from list to numeric
Regression[,4:ncol(Regression)] <- lapply(Regression[,4:ncol(Regression)], function(x) unlist(x))

# Sort by number of sampels per individual
Regression <- Regression[order(Regression$Num_Samples, decreasing=TRUE),]

# Write to file
write.csv(Regression, TABLE, row.names=TRUE)

# ______________________________________________________________________________________________________________________
#  Session data
# ______________________________________________________________________________________________________________________
save.image(file=DATA)


