# ______________________________________________________________________________________________________________________
# Summary statistics of tissue linear models
# ______________________________________________________________________________________________________________________
setwd("~/XIST/")

# Input
DATA <- "~/XIST/Gene_Tissue_121919.RData"

# Output
AVG <- "~/XIST/Tissue/Tissue_Linear_Model_Averages.csv"
SLOPES <- "~/XIST/Tissue/Tissue_Slopes_Table.csv"
MEANX_WILCOX <- "~/XIST/Tissue/Wilcox_MeanX.csv"
XIST_WILCOX <- "~/XIST/Tissue/Wilcox_XIST.csv"

# Load data
load(DATA)

# ______________________________________________________________________________________________________________________
#  Correlations summary
# ______________________________________________________________________________________________________________________
# Average R^2 of silenced genes reported in both studies for females and males
Summary.df <- data.frame(Female=colMeans(f.Regression[,2:ncol(f.Regression)]),
                         Male=colMeans(m.Regression[,2:ncol(m.Regression)]))
write.csv(Summary.df, AVG)

# ______________________________________________________________________________________________________________________
#  Table of Slopes
# ______________________________________________________________________________________________________________________
# Get list of female and male tissue type samples
f.Tissues <- rownames(f.Regression)
m.Tissues <- rownames(m.Regression)

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

write.csv(Slopes.df, SLOPES)

# ______________________________________________________________________________________________________________________
#  Wilcoxon rank sum test to test for difference in meanX expression between female and male tissues
# ______________________________________________________________________________________________________________________
# Compare MeanX in f/m tissues
Shared <- c("Brain - Cortex", "Brain - Hippocampus", "Brain - Substantia nigra",                 
            "Brain - Anterior cingulate cortex (BA24)", "Brain - Frontal Cortex (BA9)",             
            "Brain - Cerebellar Hemisphere", "Brain - Caudate (basal ganglia)",          
            "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)",          
            "Brain - Hypothalamus", "Brain - Spinal cord (cervical c-1)",       
            "Brain - Amygdala", "Brain - Cerebellum", "Adipose - Subcutaneous", 
            "Muscle - Skeletal", "Artery - Tibial", "Artery - Coronary",                        
            "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)",                                   
            "Breast - Mammary Tissue", "Skin - Not Sun Exposed (Suprapubic)", 
            "Minor Salivary Gland", "Adrenal Gland", "Thyroid", "Lung","Spleen", "Pancreas",                                  
            "Esophagus - Muscularis", "Esophagus - Mucosa","Esophagus - Gastroesophageal Junction",
            "Stomach", "Colon - Sigmoid", "Small Intestine - Terminal Ileum",         
            "Colon - Transverse", "Skin - Sun Exposed (Lower leg)",           
            "Nerve - Tibial", "Heart - Left Ventricle", "Pituitary",                                                        
            "Cells - Transformed fibroblasts", "Whole Blood",                              
            "Artery - Aorta", "Liver", "Kidney - Cortex", "Bladder")

# Keep only shared tissues
f.tmp_MeanX_XIST <- f.MeanX_Vs_XIST[Shared]
m.tmp_MeanX_XIST <- m.MeanX_Vs_XIST[Shared]

# Check lists are in same order
identical(names(f.tmp_MeanX_XIST), names(m.tmp_MeanX_XIST)) # TRUE

# Functions to make df by binding cols of unequal lengths 
# Add lists of unequal lengths as columns to df
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

Combine_Unequal <- function(DF.1, DF.2, COL){
  res <- data.frame(cbind.fill(DF.1[[COL]], DF.2[[COL]]))
  return(res)
}

# Function to rename cols
Rename_Col <- function(x, a, b){
  x <- setNames(x, c(a, b))
  return(x)
}

# Apply funcs to each df in list; Unequal lengths bc diff number of samples per tissue for females/males
# Combine female and male meanX expression into one df per tissue
MeanX <- Map(Combine_Unequal, 
             DF.1=f.tmp_MeanX_XIST, 
             DF.2=m.tmp_MeanX_XIST,
             COL='MeanX')
MeanX <- Map(Rename_Col, 
             x=MeanX, 
             a='f.MeanX', 
             b='m.MeanX')

# Function to perform two-sided Wilcoxon rank sum test
# H0: MeanX is not different b/w f/m
# alpha: 0.05
Wilcox_Func <- function(DF, COL.1, COL.2 ){
  res <- wilcox.test(DF[[COL.1]], DF[[COL.2]], conf.level=0.05) 
  return(res)
}

# Function to get p-value
Extract_pVal <- function(x){
  res <- x[[3]]
  return(res)
}

# Apply functions to list of dfs
Wilcox_MeanX <- Map(Wilcox_Func, DF=MeanX, COL.1='f.MeanX', COL.2='m.MeanX')
pVal_MeanX <- lapply(Wilcox_MeanX, Extract_pVal)

# Convert list to vector then apply Bonferoni correction for multiple testing
pVal_MeanX <- p.adjust(as.character(pVal_MeanX), "bonferroni")

# Convert vector to df
pVal.df <- data.frame(Tissue=names(Wilcox_MeanX), Wilcox_pVal=pVal_MeanX, row.names=NULL)

# Write to table
write.csv(pVal.df, file=MEANX_WILCOX)

# ______________________________________________________________________________________________________________________
#  Wilcoxon rank sum test to test for difference in XIST expression between female and male tissues
# ______________________________________________________________________________________________________________________
# H0: There is no difference in XIST expression between female and male tissues
# Apply funcs to each df in list; Unequal lengths bc diff number of samples per tissue for females/males
# Combine female and male XIST expression into one df per tissue
XIST.lst <- Map(Combine_Unequal, 
             DF.1=f.tmp_MeanX_XIST, 
             DF.2=m.tmp_MeanX_XIST,
             COL='XIST')
XIST.lst <- Map(Rename_Col, 
             x=XIST.lst, 
             a='f.XIST', 
             b='m.XIST')

# Apply functions to list of dfs
Wilcox_XIST <- Map(Wilcox_Func, DF=XIST.lst, COL.1='f.XIST', COL.2='m.XIST')
pVal_XIST <- lapply(Wilcox_MeanX, Extract_pVal)

# Convert list to vector then apply Bonferoni correction for multiple testing
pVal_XIST <- p.adjust(as.character(pVal_MeanX), "bonferroni")

# Convert vector to df
pVal_XIST.df <- data.frame(Tissue=names(Wilcox_XIST), Wilcox_pVal=pVal_XIST, row.names=NULL)

# Write to table
write.csv(pVal.df, file=XIST_WILCOX)
