# Perform linear regression on the gene count data (version 7).
# Find correlation between mean X chm expression and XIST per tissue across individuals. 

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

