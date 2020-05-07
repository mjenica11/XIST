#!/usr/bin/env Rscript

# Plots related to tissue linear model results
setwd("~/XIST/")

DATA <- "~/XIST/Tissue/Xpressed/Command/Mean_XIST.RData"
F.QQ <- "~/XIST/Tissue/Xpressed/Command/Fem_QQ_Plots.pdf"
M.QQ <- "~/XIST/Tissue/Xpressed/Command/Male_QQ_Plots.pdf"
F.SCATTER <- "~/XIST/Tissue/Xpressed/Command/Fem_Tissue_Scatter_Plots.pdf"
M.SCATTER <- "~/XIST/Tissue/Xpressed/Command/Male_Tissue_Scatter_Plots.pdf"
R2_SCATTER <- "~/XIST/Tissue/Xpressed/Command/R2_Scatter.tiff"
F.BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Command/Fem_Brain_XIST.tiff"
M.BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Command/Male_Brain_XIST.tiff"
F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Command/Fem_NotBrain_XIST.tiff"
M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Command/Male_NotBrain_XIST.tiff"
SEX_SPECIFIC <- "~/XIST/Tissue/Xpressed/Command/Sex_Specific_Violin.tiff"
SHARED_VIOLIN <- "~/XIST/Tissue/Xpressed/Command/Shared_Violin.tiff"
BRAIN_VIOLIN <- "~/XIST/Tissue/Xpressed/Command/Brain_Violin.tiff"
XIST.R2_SCATTER <- "~/XIST/Tissue/Xpressed/Command/R2_vs_CentralXIST.tiff"
VENN.INCOMPLETE <- "~/XIST/Tissue/Xpressed/Command/Venn_Incomplete.tiff"
VENN.SILENCED <- "~/XIST/Tissue/Xpressed/Command/Venn_Silenced.tiff"
VENN.VARIABLE <- "~/XIST/Tissue/Xpressed/Command/Venn_Variable.tiff"
INTER_VENN <- "~/XIST/Tissue/Xpressed/Command/Intersect_Venn.tiff"
R2_VIOLIN <- "~/XIST/Tissue/Xpressed/Command/R2_Violin.tiff"
BALATON_COR <- "~/XIST/Tissue/Xpressed/Command/Balaton_Scatter.tiff"

# Command line arguments:
# args[1]: 'mean' or 'median'
# arge[2]: gene of interest as response variable (i.e. 'XIST', 'DDX3X', or 'AR')
#args = commandArgs(trailingOnly=TRUE)
#operation = args[1]
#operation.2 = args[2]

# Constants
#if (operation == 'mean' & operation.2 == 'XIST'){
#    DATA <- "~/XIST/Tissue/Command/Mean/XIST/Mean_XIST.RData"
#    F.QQ <- "~/XIST/Tissue/Command/Mean/XIST/Fem_QQ_Plots.pdf"
#    M.QQ <- "~/XIST/Tissue/Command/Mean/XIST/Male_QQ_Plots.pdf"
#    F.SCATTER <- "~/XIST/Tissue/Command/Mean/XIST/Fem_Tissue_Scatter_Plots.pdf"
#    M.SCATTER <- "~/XIST/Tissue/Command/Mean/XIST/Male_Tissue_Scatter_Plots.pdf"
#    R2_SCATTER <- "~/XIST/Tissue/Command/Mean/XIST/R2_Scatter.tiff"
#    F.BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/XIST/Fem_Brain_XIST.tiff"
#    M.BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/XIST/Male_Brain_XIST.tiff"
#    F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/XIST/Fem_NotBrain_XIST.tiff"
#    M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/XIST/Male_NotBrain_XIST.tiff"
#    SEX_SPECIFIC <- "~/XIST/Tissue/Command/Mean/XIST/Sex_Specific_Violin.tiff"
#    SHARED_VIOLIN <- "~/XIST/Tissue/Command/Mean/XIST/Shared_Violin.tiff"
#    BRAIN_VIOLIN <- "~/XIST/Tissue/Command/Mean/XIST/Brain_Violin.tiff"
#    XIST.R2_SCATTER <- "~/XIST/Tissue/Command/Mean/XIST/R2_vs_CentralXIST.tiff"
#    VENN.INCOMPLETE <- "~/XIST/Tissue/Command/Mean/XIST/Venn_Incomplete.tiff"
#    VENN.SILENCED <- "~/XIST/Tissue/Command/Mean/XIST/Venn_Silenced.tiff"
#    VENN.VARIABLE <- "~/XIST/Tissue/Command/Mean/XIST/Venn_Variable.tiff"
#    INTER_VENN <- "~/XIST/Tissue/Command/Mean/XIST/Intersect_Venn.tiff"
#    R2_VIOLIN <- "~/XIST/Tissue/Command/Mean/XIST/R2_Violin.tiff"
#    BALATON_COR <- "~/XIST/Tissue/Command/Mean/XIST/Balaton_Scatter.tiff"
#} else if (operation == 'median' & operation.2 == 'XIST'){
#    DATA <- "~/XIST/Tissue/Command/Median/XIST/Command_040620.RData"
#    F.QQ <- "~/XIST/Tissue/Command/Median/XIST/Fem_QQ_Plots.pdf"
#    M.QQ <- "~/XIST/Tissue/Command/Median/XIST/Male_QQ_Plots.pdf"
#    F.SCATTER <- "~/XIST/Tissue/Command/Median/XIST/Fem_Tissue_Scatter_Plots.pdf"
#    M.SCATTER <- "~/XIST/Tissue/Command/Median/XIST/Male_Tissue_Scatter_Plots.pdf"
#    R2_SCATTER <- "~/XIST/Tissue/Command/Median/XIST/R2_Scatter.tiff"
#    F.BRAIN_XIST <- "~/XIST/Tissue/Command/Median/XIST/Fem_Brain_XIST.tiff"
#    M.BRAIN_XIST <- "~/XIST/Tissue/Command/Median/XIST/Male_Brain_XIST.tiff"
#    F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Median/XIST/Fem_NotBrain_XIST.tiff"
#    M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Median/XIST/Male_NotBrain_XIST.tiff"
#    SEX_SPECIFIC <- "~/XIST/Tissue/Command/Median/XIST/Sex_Specific_Violin.tiff"
#    SHARED_VIOLIN <- "~/XIST/Tissue/Command/Median/XIST/Shared_Violin.tiff"
#    BRAIN_VIOLIN <- "~/XIST/Tissue/Command/Median/XIST/Brain_Violin.tiff"
#    XIST.R2_SCATTER <- "~/XIST/Tissue/Command/Median/XIST/R2_vs_MedianXIST.tiff"
#    VENN.INCOMPLETE <- "~/XIST/Tissue/Command/Median/XIST/Venn_Incomplete.tiff"
#    VENN.SILENCED <- "~/XIST/Tissue/Command/Median/XIST/Venn_Silenced.tiff"
#    VENN.VARIABLE <- "~/XIST/Tissue/Command/Median/XIST/Venn_Variable.tiff"
#    INTER_VENN <- "~/XIST/Tissue/Command/Median/XIST/Intersect_Venn.tiff"
#    R2_VIOLIN <- "~/XIST/Tissue/Command/Median/XIST/R2_Violin.tiff"
#    BALATON_COR <- "~/XIST/Tissue/Command/Median/XIST/Balaton_Scatter.tiff"
#if (operation == 'mean' & operation.2 == 'DDX3X'){
#    DATA <- "~/XIST/Tissue/Command/Mean/DDX3X/Command_040620.RData"
#    F.QQ <- "~/XIST/Tissue/Command/Mean/DDX3X/Fem_QQ_Plots.pdf"
#    M.QQ <- "~/XIST/Tissue/Command/Mean/DDX3X/Male_QQ_Plots.pdf"
#    F.SCATTER <- "~/XIST/Tissue/Command/Mean/DDX3X/Fem_Tissue_Scatter_Plots.pdf"
#    M.SCATTER <- "~/XIST/Tissue/Command/Mean/DDX3X/Male_Tissue_Scatter_Plots.pdf"
#    R2_SCATTER <- "~/XIST/Tissue/Command/Mean/DDX3X/R2_Scatter.tiff"
#    F.BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/DDX3X/Fem_Brain_XIST.tiff"
#    M.BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/DDX3X/Male_Brain_XIST.tiff"
#    F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/DDX3X/Fem_NotBrain_XIST.tiff"
#    M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/DDX3X/Male_NotBrain_XIST.tiff"
#    SEX_SPECIFIC <- "~/XIST/Tissue/Command/Mean/DDX3X/Sex_Specific_Violin.tiff"
#    SHARED_VIOLIN <- "~/XIST/Tissue/Command/Mean/DDX3X/Shared_Violin.tiff"
#    BRAIN_VIOLIN <- "~/XIST/Tissue/Command/Mean/DDX3X/Brain_Violin.tiff"
#    XIST.R2_SCATTER <- "~/XIST/Tissue/Command/Mean/DDX3X/R2_vs_CentralXIST.tiff"
#    VENN.INCOMPLETE <- "~/XIST/Tissue/Command/Mean/DDX3X/Venn_Incomplete.tiff"
#    VENN.SILENCED <- "~/XIST/Tissue/Command/Mean/DDX3X/Venn_Silenced.tiff"
#    VENN.VARIABLE <- "~/XIST/Tissue/Command/Mean/DDX3X/Venn_Variable.tiff"
#    INTER_VENN <- "~/XIST/Tissue/Command/Mean/DDX3X/Intersect_Venn.tiff"
#    R2_VIOLIN <- "~/XIST/Tissue/Command/Mean/DDX3X/R2_Violin.tiff"
#    BALATON_COR <- "~/XIST/Tissue/Command/Mean/DDX3X/Balaton_Scatter.tiff"
#} else if (operation == 'median' & operation.2 == 'DDX3X'){
#    DATA <- "~/XIST/Tissue/Command/Median/DDX3X/Command_040620.RData"
#    F.QQ <- "~/XIST/Tissue/Command/Median/DDX3X/Fem_QQ_Plots.pdf"
#    M.QQ <- "~/XIST/Tissue/Command/Median/DDX3X/Male_QQ_Plots.pdf"
#    F.SCATTER <- "~/XIST/Tissue/Command/Median/DDX3X/Fem_Tissue_Scatter_Plots.pdf"
#    M.SCATTER <- "~/XIST/Tissue/Command/Median/DDX3X/Male_Tissue_Scatter_Plots.pdf"
#    R2_SCATTER <- "~/XIST/Tissue/Command/Median/DDX3X/R2_Scatter.tiff"
#    F.BRAIN_XIST <- "~/XIST/Tissue/Command/Median/DDX3X/Fem_Brain_XIST.tiff"
#    M.BRAIN_XIST <- "~/XIST/Tissue/Command/Median/DDX3X/Male_Brain_XIST.tiff"
#    F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Median/DDX3X/Fem_NotBrain_XIST.tiff"
#    M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Median/DDX3X/Male_NotBrain_XIST.tiff"
#    SEX_SPECIFIC <- "~/XIST/Tissue/Command/Median/DDX3X/Sex_Specific_Violin.tiff"
#    SHARED_VIOLIN <- "~/XIST/Tissue/Command/Median/DDX3X/Shared_Violin.tiff"
#    BRAIN_VIOLIN <- "~/XIST/Tissue/Command/Median/DDX3X/Brain_Violin.tiff"
#    XIST.R2_SCATTER <- "~/XIST/Tissue/Command/Median/DDX3X/R2_vs_MedianXIST.tiff"
#    VENN.INCOMPLETE <- "~/XIST/Tissue/Command/Median/DDX3X/Venn_Incomplete.tiff"
#    VENN.SILENCED <- "~/XIST/Tissue/Command/Median/DDX3X/Venn_Silenced.tiff"
#    VENN.VARIABLE <- "~/XIST/Tissue/Command/Median/DDX3X/Venn_Variable.tiff"
#    INTER_VENN <- "~/XIST/Tissue/Command/Median/DDX3X/Intersect_Venn.tiff"
#    R2_VIOLIN <- "~/XIST/Tissue/Command/Median/DDX3X/R2_Violin.tiff"
#    BALATON_COR <- "~/XIST/Tissue/Command/Median/DDX3X/Balaton_Scatter.tiff"
#if (operation == 'mean' & operation.2 == 'AR'){
#    DATA <- "~/XIST/Tissue/Command/Mean/AR/Command_040620.RData"
#    F.QQ <- "~/XIST/Tissue/Command/Mean/AR/Fem_QQ_Plots.pdf"
#    M.QQ <- "~/XIST/Tissue/Command/Mean/AR/Male_QQ_Plots.pdf"
#    F.SCATTER <- "~/XIST/Tissue/Command/Mean/AR/Fem_Tissue_Scatter_Plots.pdf"
#    M.SCATTER <- "~/XIST/Tissue/Command/Mean/AR/Male_Tissue_Scatter_Plots.pdf"
#    R2_SCATTER <- "~/XIST/Tissue/Command/Mean/AR/R2_Scatter.tiff"
#    F.BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/AR/Fem_Brain_XIST.tiff"
#    M.BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/AR/Male_Brain_XIST.tiff"
#    F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/AR/Fem_NotBrain_XIST.tiff"
#    M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Mean/AR/Male_NotBrain_XIST.tiff"
#    SEX_SPECIFIC <- "~/XIST/Tissue/Command/Mean/AR/Sex_Specific_Violin.tiff"
#    SHARED_VIOLIN <- "~/XIST/Tissue/Command/Mean/AR/Shared_Violin.tiff"
#    BRAIN_VIOLIN <- "~/XIST/Tissue/Command/Mean/AR/Brain_Violin.tiff"
#    XIST.R2_SCATTER <- "~/XIST/Tissue/Command/Mean/AR/R2_vs_CentralXIST.tiff"
#    VENN.INCOMPLETE <- "~/XIST/Tissue/Command/Mean/AR/Venn_Incomplete.tiff"
#    VENN.SILENCED <- "~/XIST/Tissue/Command/Mean/AR/Venn_Silenced.tiff"
#    VENN.VARIABLE <- "~/XIST/Tissue/Command/Mean/AR/Venn_Variable.tiff"
#    INTER_VENN <- "~/XIST/Tissue/Command/Mean/AR/Intersect_Venn.tiff"
#    R2_VIOLIN <- "~/XIST/Tissue/Command/Mean/AR/R2_Violin.tiff"
#    BALATON_COR <- "~/XIST/Tissue/Command/Mean/AR/Balaton_Scatter.tiff"
#} else if (operation == 'median' & operation.2 == 'AR'){
#    DATA <- "~/XIST/Tissue/Command/Median/AR/Command_040620.RData"
#    F.QQ <- "~/XIST/Tissue/Command/Median/AR/Fem_QQ_Plots.pdf"
#    M.QQ <- "~/XIST/Tissue/Command/Median/AR/Male_QQ_Plots.pdf"
#    F.SCATTER <- "~/XIST/Tissue/Command/Median/AR/Fem_Tissue_Scatter_Plots.pdf"
#    M.SCATTER <- "~/XIST/Tissue/Command/Median/AR/Male_Tissue_Scatter_Plots.pdf"
#    R2_SCATTER <- "~/XIST/Tissue/Command/Median/AR/R2_Scatter.tiff"
#    F.BRAIN_XIST <- "~/XIST/Tissue/Command/Median/AR/Fem_Brain_XIST.tiff"
#    M.BRAIN_XIST <- "~/XIST/Tissue/Command/Median/AR/Male_Brain_XIST.tiff"
#    F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Median/AR/Fem_NotBrain_XIST.tiff"
#    M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Command/Median/AR/Male_NotBrain_XIST.tiff"
#    SEX_SPECIFIC <- "~/XIST/Tissue/Command/Median/AR/Sex_Specific_Violin.tiff"
#    SHARED_VIOLIN <- "~/XIST/Tissue/Command/Median/AR/Shared_Violin.tiff"
#    BRAIN_VIOLIN <- "~/XIST/Tissue/Command/Median/AR/Brain_Violin.tiff"
#    XIST.R2_SCATTER <- "~/XIST/Tissue/Command/Median/AR/R2_vs_MedianXIST.tiff"
#    VENN.INCOMPLETE <- "~/XIST/Tissue/Command/Median/AR/Venn_Incomplete.tiff"
#    VENN.SILENCED <- "~/XIST/Tissue/Command/Median/AR/Venn_Silenced.tiff"
#    VENN.VARIABLE <- "~/XIST/Tissue/Command/Median/AR/Venn_Variable.tiff"
#    INTER_VENN <- "~/XIST/Tissue/Command/Median/AR/Intersect_Venn.tiff"
#    R2_VIOLIN <- "~/XIST/Tissue/Command/Median/AR/R2_Violin.tiff"
#    BALATON_COR <- "~/XIST/Tissue/Command/Median/AR/Balaton_Scatter.tiff"
#}

# Load libraries
library(readr) 
library(plyr) # Load plyr before dplyer
library(dplyr) 
library(data.table) 
library(stringr)
library(broom)
library(rjson)
library(plotly)
library(grDevices)
library(gridExtra)
library(reshape2)
library(ggunchained)
library(ggVennDiagram)
library(ggplot2)

# Load session data
load(DATA)

# Vectors for subsetting dfs
Brain_Tissues <- c("Brain - Cortex", "Brain - Hippocampus", "Brain - Substantia nigra",                 
                   "Brain - Anterior cingulate cortex (BA24)", "Brain - Frontal Cortex (BA9)",             
                   "Brain - Cerebellar Hemisphere", "Brain - Caudate (basal ganglia)",          
                   "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)",          
                   "Brain - Hypothalamus", "Brain - Spinal cord (cervical c-1)",       
                   "Brain - Amygdala", "Brain - Cerebellum")

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

# ______________________________________________________________________________________________________________________
#  QQ Plots
# ______________________________________________________________________________________________________________________
# Functon to plot LM residuals
QQ_Func <- function(LM, PREDICTOR, RESPONSE, SEX, TISSUE){
  plot <- qqnorm(y = LM[['model']][['Gene']], # predictor/independent var
                 ylab = paste("Sample", PREDICTOR, "TPMs"),
                 xlab = "Theoretical Quantiles",
                 main = paste("lm(", RESPONSE, " ~", PREDICTOR, ") in", SEX, TISSUE))
          qqline(LM[['model']][['Gene']])
          return(plot)
}
pdf(F.QQ)
Map(QQ_Func, LM=lm_f.CentralX_Gene,
             PREDICTOR="XIST",
             RESPONSE="MeanX",
             SEX=c("Female"), 
             TISSUE=names(lm_f.CentralX_Gene))
dev.off()
pdf(M.QQ)
Map(QQ_Func, LM=lm_m.CentralX_Gene,
             PREDICTOR="XIST",
             RESPONSE="MeanX",
             SEX=c("Male"), 
             TISSUE=names(lm_m.CentralX_Gene))
dev.off()

# ______________________________________________________________________________________________________________________
#  Scatter Plots
# ______________________________________________________________________________________________________________________
# Find the max x/y values to set as x/y lim
Max_Func <- function(x, y){
    lim <- max(x[[y]])
    res <- round(max(unlist(lim)))
    return(res)
}
# X limits
f.xmax <- max(unlist(Map(Max_Func, x=f.CentralX_Vs_Gene, y='Gene')))
m.xmax <- max(unlist(Map(Max_Func, x=m.CentralX_Vs_Gene, y='Gene')))
# Y limits
f.ymax <- max(unlist(Map(Max_Func, x=f.CentralX_Vs_Gene, y='Central')))
m.ymax <- max(unlist(Map(Max_Func, x=m.CentralX_Vs_Gene, y='Central')))

# Return the larger X limit + 5 to leave room near the plot boundaries
x.lim <- max(f.xmax, m.xmax) + 5
y.lim <- max(f.ymax, m.ymax) + 5

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, R2, PVAL, XLAB, YLAB){
  plot(LM[['model']][['Gene']], 
       LM[['model']][['Central']], 
       main=TITLE, 
       xlab=XLAB, 
       ylab=YLAB,
       xlim=c(0, x.lim), 
       ylim=c(0, y.lim))
  legend("bottomright", 
         bty="n", 
         legend=paste("R^2: ", 
                      format(R2, digits=3), 
                      "; p_Val: ", 
                      format(PVAL, digits=3)))
  abline(LM)
}

# Print plots
pdf(F.SCATTER)
f.Scatter <- Map(Scatter_Func, 
                 LM=lm_f.CentralX_Gene, 
                 TITLE=names(lm_f.CentralX_Gene), 
                 XLAB='XIST',
                 YLAB='Mean X',
                 R2=f.Regression[['R2_CentralX']],
                 PVAL=f.Regression[['pval_CentralX']])
dev.off()

pdf(M.SCATTER)
m.Scatter <- Map(Scatter_Func, 
                 LM=lm_m.CentralX_Gene, 
                 TITLE=names(lm_m.CentralX_Gene), 
                 XLAB='XIST',
                 YLAB='Mean X',
                 R2=m.Regression[['R2_CentralX']],
                 PVAL=m.Regression[['pval_CentralX']])
dev.off()

# ______________________________________________________________________________________________________________________
# R2 scatter plots of Balaton gene categories-merged 
# ______________________________________________________________________________________________________________________
# Subset regression df to contain only relevant data
Bal.df <- f.Regression[,c('Tissue', 'R2_CentralX', 'R2_Bal_Silenced', 'R2_Bal_Variable', 'R2_Bal_Incomplete')]

# Convert df to 'tall' format; i.e. combine n cols to one column and repeat the col containing the rownames (Tissue) n times 
Bal.df <- melt(Bal.df, id.vars='Tissue')

# Rename cols 
setnames(Bal.df, old=c('variable', 'value'), new=c('LM_Test', 'R2_Value'))

ggplot(Bal.df, aes(x=Tissue, y=R2_Value)) + 
 geom_point(aes(color=LM_Test)) +
# scale_shape_manual(values=c(21,22)) +
 scale_color_manual(values=c('black', 'blue', 'green', 'red')) +
 ggtitle('Scatter plot of R2 of Various Gene Sets in Female Tissues') +
 xlab('Tissue Type') +
 ylab('R2 Value') +
 ylim(c(0,1)) +
 theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(BALATON_COR, device="tiff")

# ______________________________________________________________________________________________________________________
# Prepare df for plotting R2 scatter plots
# ______________________________________________________________________________________________________________________
# Combine relevant cols from dfs
# df of R2_CentralX and Mean_XIST for both females and males
Subset_f.df <- f.Regression[ ,c('Tissue', 'R2_CentralX', 'Mean_XIST', 'pval_CentralX')]
Common_f.df <- Subset_f.df[Subset_f.df$Tissue %in% Shared, ]
Common_f.df$Sex <- 'Female'

# Add values from males
Subset_m.df <- m.Regression[ ,c('Tissue', 'R2_CentralX', 'Mean_XIST', 'pval_CentralX')]
Common_m.df <- Subset_m.df[Subset_m.df$Tissue %in% Shared, ]
Common_m.df$Sex <- 'Male'

# Combine 
Common.df <- rbind(Common_f.df, Common_m.df)

# x and y lims
max(Common.df$Mean_XIST) # 141.7757
max(Common.df$R2_CentralX) # 0.635009; in original version of script it was 0.5555296

# Add col with pval factor (p> or < 0.05)
Common.df$p_Factor <- as.factor(Common.df$pval_CentralX<0.05)
head(Common.df)

# ______________________________________________________________________________________________________________________
#  Scatter plot of Mean X vs XIST R2 values across female/male tissues
# ______________________________________________________________________________________________________________________
# Plot females and males together
ggplot(Common.df, aes(x=Tissue, y=R2_CentralX)) + 
 geom_point(aes(shape=Sex, fill=Sex, alpha=p_Factor)) +
 scale_shape_manual(values=c(21,22)) +
 scale_fill_manual(values=c('blue', 'green')) +
 ggtitle('Scatter plot of R2 (MeanX ~ XIST) in Female and Male Tissues') +
 xlab('Tissue Type') +
 ylab('R2 (MeanX~ XIST)') +
 ylim(c(0,1)) +
 theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(R2_SCATTER, device="tiff")

# ______________________________________________________________________________________________________________________
#  Scatter plot of R^2 of MeanX and XIST vs XIST
# ______________________________________________________________________________________________________________________
# Scatter plot
ggplot(Common.df, aes(x=Mean_XIST, y=R2_CentralX)) +
  geom_point(aes(shape=Sex, fill=Sex, alpha=p_Factor)) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=c('blue', 'green')) +
  ggtitle('R^2 for XIST and Mean X vs Mean XIST') +
  xlab('Mean XIST') +
  ylab('R2 for XIST and Mean X') +
  xlim(0,150) +
  ylim(0,1) 
ggsave(XIST.R2_SCATTER, device="tiff")

# ______________________________________________________________________________________________________________________
#  Prepare dfs of meanX and XIST expression in f + m brain tissues, sex-shared tissues, 
#  and sex-specific tissues for violin plots 
# ______________________________________________________________________________________________________________________
# Prepare df
# Collapse list of dfs into one df of just the meanX counts
f.df <- ldply(f.CentralX_Vs_Gene, data.frame)
f.df$Tissue <- f.df$.id
f.df$.id <- NULL
f.df$Sex <- 'Female'
head(f.df)

m.df <- ldply(m.CentralX_Vs_Gene, data.frame)
m.df$Tissue <- m.df$.id
m.df$.id <- NULL
m.df$Sex <- 'Male'
head(m.df)

# Seperate into brain, and two dfs non-brain tissues
f.Not_Brain <- setdiff(names(lm_f.CentralX_Gene), Brain_Tissues)
m.Not_Brain <- setdiff(names(lm_m.CentralX_Gene), Brain_Tissues)

# Subset MeanX dfs for plots
Brain_f.df <- f.df[f.df$Tissue %in% Brain_Tissues,]
Not_Brain_f.df <- f.df[f.df$Tissue %in% f.Not_Brain,]

Brain_m.df <- m.df[m.df$Tissue %in% Brain_Tissues,]
Not_Brain_m.df <- m.df[m.df$Tissue %in% m.Not_Brain,]

# Combine MeanX dfs into one df
# Brain tissues only
Brain.df <- rbind(Brain_f.df, Brain_m.df)
head(Brain.df); tail(Brain.df)
nrow(Brain.df) == nrow(Brain_f.df) + nrow(Brain_m.df)

# Sex-specific tissue types list
f.Only <- setdiff(f.Not_Brain, m.Not_Brain)
m.Only <- setdiff(m.Not_Brain, f.Not_Brain)

# Subset MeanX dfs and combine to get df of one sex, not brain tissues
Only_f.df <- f.df[f.df$Tissue %in% f.Only,]
Only_m.df <- m.df[m.df$Tissue %in% m.Only,]

# Add col indicating sex to MeanX and XIST dfs
Only_f.df$Sex <- 'Female'
Only_m.df$Sex <- 'Male'

# Join for both
Sex_Specific.df <- rbind(Only_f.df, Only_m.df)
head(Sex_Specific.df);tail(Sex_Specific.df)

# Create df of tissues excluding brain and sex-specific tissues for both
f.Shared <- f.df[f.df$Tissue %in% Shared,]
m.Shared <- m.df[m.df$Tissue %in% Shared,]

# Join for both
Shared.df <- rbind(f.Shared, m.Shared)
head(Shared.df); tail(Shared.df)

# ______________________________________________________________________________________________________________________
#  MeanX Split-violin plots 
# ______________________________________________________________________________________________________________________
# ylims
# Set to 90 TPM to keep yaxis proportional between plots
max(Brain.df$Central) # 76.295
max(Shared.df$Central) # 86.90617

Split_Violin <- function(DF, TITLE){
   ggplot(DF, aes(x = Tissue, y = Central, color = Sex, fill=Sex)) +
      geom_split_violin(color='black') +
      scale_fill_manual(values=c('blue', 'green')) +
      ylim(c(0,90)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(TITLE)
}
# Split violin plot of brain tissues
Split_Violin(DF=Brain.df, TITLE="Split violin plot of mean X expression in brain tissues")
ggsave(BRAIN_VIOLIN)

# Split violin plot of sex-shared tissues
Split_Violin(DF=Shared.df, TITLE="Split violin plot of mean X expression in sex-shared tissues")
ggsave(SHARED_VIOLIN, device="tiff")

# ______________________________________________________________________________________________________________________
#  Sex-specific Tissues Violin Plot
# ______________________________________________________________________________________________________________________
# ylim
max(Sex_Specific.df$Central) # 76.0642
levels(as.factor(Sex_Specific.df$Tissue))

# Plot
# Explicitely set the order of samples; default is alphabetical
ggplot(Sex_Specific.df, 
       aes(x = factor(Tissue, level=c("Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", 
                                      "Ovary", "Uterus", "Vagina", "Prostate", "Testis")),
           y = Central,
           color=Sex,
           fill=Sex)) +
   geom_violin(color='black') +
   scale_fill_manual(values=c('blue', 'green')) +
   ylim(c(0,80)) + 
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
   ggtitle("Violin plot of mean X expression in sex-specific tissues")
ggsave(SEX_SPECIFIC, width=5, height=5, device="tiff")

# ______________________________________________________________________________________________________________________
#  XIST violin plots 
# ______________________________________________________________________________________________________________________
# Log transform just female XIST, since high variance
head(Brain_f.df); tail(Brain_f.df)
Brain_f.df[,1:2] <- log10(Brain_f.df[,1:2]) 
Not_Brain_f.df[,1:2] <- log10(Not_Brain_f.df[,1:2]) 

#Violin plot function
XIST_Violin <- function(DF, YLIM, TITLE){
    ggplot(DF, 
           aes(x = Tissue,
               y = Gene,
               color = Tissue,
               fill = Tissue)) +
        geom_violin(color='black') +
        ylim(YLIM) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(TITLE) +
        guides(fill=FALSE)
}
XIST_Violin(DF=Brain_f.df, 
            YLIM=c(0, ceiling(max(Brain_f.df$Gene))), 
            TITLE="Violin plot of log10 XIST expression across female brain tissues")
ggsave(F.BRAIN_XIST, device="tiff")
XIST_Violin(DF=Not_Brain_f.df, 
            YLIM=c(0, ceiling(max(Not_Brain_f.df$Gene))), 
            TITLE="Violin plot of log10 XIST expression across female tissues")
ggsave(F.NOT_BRAIN_XIST, device="tiff")
XIST_Violin(DF=Brain_m.df, 
            YLIM=c(0, ceiling(max(Brain_m.df$Gene))), 
            TITLE="Violin plot of XIST expression across male brain tissues")
ggsave(M.BRAIN_XIST, device="tiff")
XIST_Violin(DF=Not_Brain_m.df, 
            YLIM=c(0, ceiling(max(Not_Brain_m.df$Gene))), 
            TITLE="Violin plot of XIST expression across male tissues")
ggsave(M.NOT_BRAIN_XIST, device="tiff")

# ______________________________________________________________________________________________________________________
#  Violin plot: R2 by gene classification scheme per tissue
# ______________________________________________________________________________________________________________________
# Reshape correlation df 
Category_Lst <- c('Silenced In Both', 'Silenced In Tukainen', 'Silenced In Balaton', 'Silenced In At Least One', 
                  'Silenced Immune Genes', 'Variable In At Least One', 'Variable In Tukainen',
                  'Variable In Balaton', 'Variable Immune Genes', 'Incomplete In At Least One', 'Incomplete in Tukainen',
                  'Invariable In Balaton', 'Invariable Immune Genes', 'All Evaluated In Both', 'Not Evaluated In Either', 
                  'Immune Genes Not Evaluated', 'PAR Genes In Balaton')

Columns_Lst <- c('R2_Silenced', 'R2_Tuk_Silenced', 'R2_Bal_Silenced', 'R2_One_Silenced',
                 'R2_Immune_Silenced', 'R2_One_Variable', 'R2_Tuk_Variable', 'R2_Bal_Variable', 
                 'R2_Immune_Variable', 'R2_One_Incomplete', 'R2_Tuk_Incomplete', 'R2_Bal_Incomplete', 
                 'R2_Immune_Incomplete', 'R2_All_Eval', 'R2_Not_Eval', 'R2_Immune_Not_Eval', 'R2_PAR')

# Function to subset and reshape df 
Fem_Subset <- function(SUBSET, CAT){
  res.df <- data.frame(Tissue=f.Regression[['Tissue']], R2=f.Regression[[SUBSET]], Category=CAT)
  return(res.df)
}
Sub_Lst.f <- Map(Fem_Subset, 
                 SUBSET=Columns_Lst, 
                 CAT=Category_Lst)

Male_Subset <- function(SUBSET, CAT){
  res.df <- data.frame(Tissue=m.Regression[['Tissue']], R2=m.Regression[[SUBSET]], Category=CAT)
  return(res.df)
}
Sub_Lst.m <- Map(Fem_Subset, 
                 SUBSET=Columns_Lst, 
                 CAT=Category_Lst)

# rbind list of dfs and add column indicating sex
Combined.f <- rbindlist(Sub_Lst.f)
Combined.f$Sex <- 'Female'

Combined.m <- rbindlist(Sub_Lst.m)
Combined.m$Sex <- 'Male'

# bind dfs
Correlations.df <- rbind(Combined.f, Combined.m)
head(Correlations.df); tail(Correlations.df)

# Plot
ggplot(Correlations.df, aes(x = Category, y = R2, color = Sex, fill=Sex)) +
   geom_split_violin(color='black') +
   scale_fill_manual(values=c('blue', 'green')) +
   ylim(c(0,1)) +
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
   ggtitle("Split violin plot of R2 per gene subset")
ggsave(R2_VIOLIN, width=5, height=5, device="tiff")

# ______________________________________________________________________________________________________________________
#  Venn diagram; overlap between gene classification schemes
# ______________________________________________________________________________________________________________________
summary(Gene_Lst)

# Plot
tiff(VENN.INCOMPLETE)
ggVennDiagram(Gene_Lst[c('Incomplete_In_Tukiainen', 'Incomplete_In_Balaton')]) +
    scale_fill_gradient(low='blue', high='red') +
    ggtitle('Venn Diagram of Incomplete Genes')
dev.off()
tiff(VENN.VARIABLE)
ggVennDiagram(Gene_Lst[c('Variable_In_Tukiainen', 'Variable_In_Balaton')]) +
    scale_fill_gradient(low='blue', high='red') +
    ggtitle('Venn Diagram of Variable Genes')
dev.off()
tiff(VENN.SILENCED)
ggVennDiagram(Gene_Lst[c('Silenced_In_Tukiainen', 'Silenced_In_Balaton')]) +
    scale_fill_gradient(low='blue', high='red') +
    ggtitle('Venn Diagram of Silenced Genes')
dev.off()
# ______________________________________________________________________________________________________________________
# Venn diagram; intersection of gene intersections 
# ______________________________________________________________________________________________________________________
# Make a list of the intersection of the silenced, variable, and incomplete gene sets in both studies
Incomp <- intersect(Gene_Lst[[12]], Gene_Lst[[15]])
Var <- intersect(Gene_Lst[[13]], Gene_Lst[[17]])
Sil <- intersect(Gene_Lst[[11]], Gene_Lst[[14]])

# Combine into one list
Inter <- list(Incomp, Var, Sil)
names(Inter) <- c("Incomplete", "Variable", "Silenced")
summary(Inter)

# Plot
ggVennDiagram(Inter[1:3]) +
    scale_fill_gradient(low='blue', high='red') +
    ggtitle('Venn Diagram of Intersection of Intersection of Gene Categories')
ggsave(INTER_VENN, device="tiff")

