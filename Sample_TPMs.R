# Check how many genes remain after stict TPM filtering
setwd("~/XIST/")

# Constants
COUNTS <- "~/XIST/Files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct" # TPM normalized
METRICS <- "~/XIST/Files/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv"
PHENOTYPES <- "~/XIST/Files/GTEX_v7_Annotations_SubjectPhenotypesDS.txt"
GENCODE <- "gencode.v19.genes.v7.patched_contigs.gff3" # has to be in working dir
GENE_LST <- "~/XIST/Files/X_Genes_Status.json"

# Results
DATA <- "~/XIST/XIST_Indiv_121819.RData"

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

# Keep genes with TPM > 5, excluding XIST
# Store XIST count row
tmp <- Gene_Cts[Gene_Cts$gene_name=='XIST',]

tester <- Gene_Cts[1:20,1:20]
tester.1 <- tester[apply(tester[,3:ncol(tester)], MARGIN=1, function(x) all(x>0)),]
View(tester.1)

# Filter for genes with TPM > 5, excluding XIST
Gene_Cts.5 <- Gene_Cts[apply(Gene_Cts[,3:ncol(Gene_Cts)], MARGIN=1, function(x) all(x>5)),] # Skip first 2 cols
Gene_Cts.3 <- Gene_Cts[apply(Gene_Cts[,3:ncol(Gene_Cts)], MARGIN=1, function(x) all(x>3)),] # Skip first 2 cols
Gene_Cts.1 <- Gene_Cts[apply(Gene_Cts[,3:ncol(Gene_Cts)], MARGIN=1, function(x) all(x>1)),] # Skip first 2 cols
 
# How many genes before TPM filter?
nrow(Gene_Cts) # 56,202
 
# How many genes left after TPM filter?
nrow(Gene_Cts.1) # TPM > 1 --> 2,591
nrow(Gene_Cts.3) # TPM > 3 --> 792
nrow(Gene_Cts.5) # TPM > 5 --> 404
 
# Add XIST back to df
Gene_Cts <- rbind(Gene_Cts, tmp)
