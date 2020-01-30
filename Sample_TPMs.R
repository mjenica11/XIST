#!/usr/bin/env Rscript

# Check how many genes remain after stict TPM filtering
setwd("~/XIST/")

# Constants
COUNTS <- "~/XIST/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct" # TPM normalized
METRICS <- "~/XIST/Files/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"
PHENOTYPES <- "~/XIST/Files/GTEX_v7_Annotations_SubjectPhenotypesDS.txt"
GENCODE <- "gencode.v19.genes.v7.patched_contigs.gff3" # has to be in working dir
GENE_LST <- "~/XIST/Files/X_Genes_Status.json"

# Results
DATA <- "~/XIST/Tissue/TPM/TPM_Filtering_012320.RData"

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
# Organize count dfs by tissue type, then split by sex. Drop replicates. 
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
# Filter by TPM 
# ______________________________________________________________________________________________________________________
###### TEST #######
# Subset each df from list to make a test object
# f.Subset_Cts <- lapply(f.Tissue_Counts, function(x) x[1:50,1:5])
# m.Subset_Cts <- lapply(m.Tissue_Counts, function(x) x[1:50,1:5])

# Function to filter list of count dfs by tpm value
Filter_Cts <- function(DF, TPM){
  DF <- DF[apply(DF[,3:ncol(DF)], MARGIN=1, function(x) all(x>TPM)),]
  DF <- DF[,1:2] # Keep only gene name and id cols; drop all sample counts
  return(DF)
}

# Keep genes with TPM > 0; counts are seperated by sex
f.Filtered_0 <- Map(Filter_Cts, DF=f.Tissue_Counts, TPM=0)
m.Filtered_0 <- Map(Filter_Cts, DF=m.Tissue_Counts, TPM=0)
# Keep genes with TPM > 1
f.Filtered_1 <- Map(Filter_Cts, DF=f.Tissue_Counts, TPM=1)
m.Filtered_1 <- Map(Filter_Cts, DF=m.Tissue_Counts, TPM=1)
# Keep genes with TPM > 3
f.Filtered_3 <- Map(Filter_Cts, DF=f.Tissue_Counts, TPM=3)
m.Filtered_3 <- Map(Filter_Cts, DF=m.Tissue_Counts, TPM=3)
# Keep genes with TPM > 5
f.Filtered_5 <- Map(Filter_Cts, DF=f.Tissue_Counts, TPM=5)
m.Filtered_5 <- Map(Filter_Cts, DF=m.Tissue_Counts, TPM=5)
# Keep genes with TPM > 7
f.Filtered_7 <- Map(Filter_Cts, DF=f.Tissue_Counts, TPM=7)
m.Filtered_7 <- Map(Filter_Cts, DF=m.Tissue_Counts, TPM=7)
# Keep genes with TPM > 10
f.Filtered_10 <- Map(Filter_Cts, DF=f.Tissue_Counts, TPM=10) 
m.Filtered_10 <- Map(Filter_Cts, DF=m.Tissue_Counts, TPM=10)

# Sanity check; make sure all dfs in list are there and in the same order
identical(names(f.Tissue_Counts), names(f.Filtered_5)) # TRUE
identical(names(m.Tissue_Counts), names(m.Filtered_5)) # TRUE

# Make vector of files names; sex sample tissue type
# First, remove spaces and dashes from the tissue names
Fix_Name <- function(x){
  x <- gsub(" ", "", x)
  x <- gsub("-", "_", x)
  return(x)
}
f.Tissues <- names(f.Filtered_5)
f.Tissues <- sapply(f.Tissues, Fix_Name)

m.Tissues <- names(m.Filtered_5)
m.Tissues <- sapply(m.Tissues, Fix_Name)

# Make character vector listing sex that is equal in length to the number of dfs
Fems <- replicate(length(f.Tissues), "Fem")
Males <- replicate(length(m.Tissues), "Male")

# Function to write each df to csv; Name includes tissue type and sex
Name_Files <- function(DF, SUBDIR, SEX, TISSUE){
  res <- write.csv(DF,
                   file=paste0("/Users/Mollie/XIST/Tissue/TPM/", # Directory
                               SEX,
                               "_",
                               TISSUE,
                               ".csv"),
                   row.names = FALSE) # add sample sex to file name
  return(res)
}
# Only write results for TPM > 5 for now
Map(Name_Files, DF=f.Filtered_5, SEX=Fems, TISSUE=f.Tissues)
Map(Name_Files, DF=m.Filtered_5, SEX=Males, TISSUE=m.Tissues)

# ______________________________________________________________________________________________________________________
# Summarize TPM filtering results
# ______________________________________________________________________________________________________________________
# How many genes before TPM filter?
# lapply(f.Filtered_5, nrow)

 



