#!/usr/bin/env Rscript

# Plots related to tissue linear model results
setwd("~/XIST/")

# Constants
DATA <- "~/XIST/Tissue/Xpressed/Median/Xpressed_012720.RData"
F.QQ <- "~/XIST/Tissue/Xpressed/Median/Fem_QQ_Plots.pdf"
M.QQ <- "~/XIST/Tissue/Xpressed/Median/Male_QQ_Plots.pdf"
F.SCATTER <- "~/XIST/Tissue/Xpressed/Median/Fem_Tissue_Scatter_Plots.pdf"
M.SCATTER <- "~/XIST/Tissue/Xpressed/Median/Male_Tissue_Scatter_Plots.pdf"
R2_SCATTER <- "~/XIST/Tissue/Xpressed/Median/R2_Scatter.tiff"
F.BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Median/Fem_Brain_XIST.tiff"
M.BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Median/Male_Brain_XIST.tiff"
F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Median/Fem_NotBrain_XIST.tiff"
M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Median/Male_NotBrain_XIST.tiff"
SEX_SPECIFIC <- "~/XIST/Tissue/Xpressed/Median/Sex_Specific_Violin.tiff"
SHARED_VIOLIN <- "~/XIST/Tissue/Xpressed/Median/Shared_Violin.tiff"
NOT_BRAIN_VIOLIN <- "~/XIST/Tissue/Xpressed/Median/NotBrain_Violin.tiff"
BRAIN_VIOLIN <- "~/XIST/Tissue/Xpressed/Median/Brain_Violin.tiff"
XIST.R2_SCATTER <- "~/XIST/Tissue/Xpressed/Median/R2_XIST_and_MedX_Vs_MedXIST.tiff"
VENN_INCOMPLETE <- "~/XIST/Tissue/Xpressed/Median/Venn_Incomplete.tiff"
VENN_SILENCED <- "~/XIST/Tissue/Xpressed/Median/Venn_Silenced.tiff"
VENN_VARIABLE <- "~/XIST/Tissue/Xpressed/Median/Venn_Variable.tiff"
R2_VIOLIN <- "~/XIST/Tissue/Xpressed/Median/R2_Violin.tiff"
BALATON_COR <- "~/XIST/Tissue/Xpressed/Median/Balaton_Scatter.tiff"
INTER_VENN <- "~/XIST/Tissue/Xpressed/Median/Intersection_Venn.tiff"

# Load libraries
library(readr) 
library(refGenome) # Parse .gff file
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
QQ_Func <- function(LM, SEX, TISSUE){
  plot <- qqnorm(y = LM[['model']][['XIST']], # predictor/independent var
                 ylab = "Sample XIST TPMs",
                 xlab = "Theoretical Quantiles",
                 main = paste("lm(MedX ~ XIST) in", SEX, TISSUE))
  qqline(LM[['model']][['XIST']])
  return(plot)
}

pdf(F.QQ)
Map(QQ_Func, LM=lm_f.MedX_XIST, SEX=c("Female"), TISSUE=names(lm_f.MedX_XIST))
dev.off()
pdf(M.QQ)
Map(QQ_Func, LM=lm_m.MedX_XIST, SEX=c("Male"), TISSUE=names(lm_f.MedX_XIST))
dev.off()

# ______________________________________________________________________________________________________________________
#  Scatter Plots
# ______________________________________________________________________________________________________________________
# Find the max x/y values to set as x/y lim
Max_Func <- function(x, y){
  lim <- max(x[[y]])
  return(lim)
}

# x limits
f.xmax <- round(max(unlist(Map(Max_Func, x=f.MedX_Vs_XIST, y='XIST'))))
m.xmax <- max(unlist(Map(Max_Func, x=m.MedX_Vs_XIST, y='XIST')))
# y limits
f.ymax <- round(max(unlist(Map(Max_Func, x=f.MedX_Vs_XIST, y='MedX'))))
m.ymax <- max(unlist(Map(Max_Func, x=m.MedX_Vs_XIST, y='MedX')))

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, XMAX, YMAX, R2, PVAL){
  plot(LM[['model']][['XIST']],
       LM[['model']][['MedX']], 
       main=TITLE, 
       xlab='XIST', 
       ylab='Median X chromosome',
       xlim=c(0, XMAX), 
       ylim=c(0, YMAX))
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
                 LM=lm_f.MedX_XIST, 
                 TITLE=names(lm_f.MedX_XIST), 
                 XMAX=f.xmax, 
                 YMAX=f.ymax, 
                 R2=f.Regression[['R2_MedX']],
                 PVAL=f.Regression[['pval_MedX']])
dev.off()

# set x lim to male max(XIST) but keep y lim as female max(median X)
pdf(M.SCATTER)
m.Scatter <- Map(Scatter_Func, 
                 LM=lm_m.MedX_XIST, 
                 TITLE=names(lm_m.MedX_XIST), 
                 XMAX=m.xmax, 
                 YMAX=f.ymax, 
                 R2=m.Regression[['R2_MedX']],
                 PVAL=m.Regression[['pval_MedX']])
dev.off()

# ______________________________________________________________________________________________________________________
# R2 scatter plots of Balaton gene categories-merged 
# ______________________________________________________________________________________________________________________
# Subset regression df to contain only relevant data
Bal.df <- f.Regression[,c('Tissue', 'R2_MedX', 'R2_Bal_Silenced_Median', 'R2_Bal_Variable_Median', 'R2_Bal_Incomplete_Median')]

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
# df of R2_median X and median_XIST for both females and males
Subset_f.df <- f.Regression[ ,c('Tissue', 'R2_MedX', 'Median_XIST', 'pval_MedX')]
Common_f.df <- Subset_f.df[Subset_f.df$Tissue %in% Shared, ]
Common_f.df$Sex <- 'Female'

# Add values from males
Subset_m.df <- m.Regression[ ,c('Tissue', 'R2_MedX', 'Median_XIST', 'pval_MedX')]
Common_m.df <- Subset_m.df[Subset_m.df$Tissue %in% Shared, ]
Common_m.df$Sex <- 'Male'

# Combine 
Common.df <- rbind(Common_f.df, Common_m.df)

# x and y lims
max(Common.df$Median_XIST) # 135.85
max(Common.df$R2_MedX) # 0.8402523

# Add col with pval factor (p> or < 0.05)
Common.df$p_Factor <- as.factor(Common.df$pval_MedX<0.05)
head(Common.df)

# ______________________________________________________________________________________________________________________
#  Scatter plot of median X vs XIST R2 values across female/male tissues
# ______________________________________________________________________________________________________________________
# Plot females and males together
ggplot(Common.df, aes(x=Tissue, y=R2_MedX)) + 
 geom_point(aes(shape=Sex, fill=Sex, alpha=p_Factor)) +
 scale_shape_manual(values=c(21,22)) +
 scale_fill_manual(values=c('blue', 'green')) +
 ggtitle('Scatter plot of R2 (MedianX ~ XIST) in Female and Male Tissues') +
 xlab('Tissue Type') +
 ylab('R2 (MedianX~ XIST)') +
 ylim(c(0,1)) +
 theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(R2_SCATTER, device="tiff")

# ______________________________________________________________________________________________________________________
#  Scatter plot of R^2 of median X and XIST vs XIST
# ______________________________________________________________________________________________________________________
# Scatter plot
ggplot(Common.df, aes(x=Median_XIST, y=R2_MedX)) +
  geom_point(aes(shape=Sex, fill=Sex, alpha=p_Factor)) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=c('blue', 'green')) +
  ggtitle('R^2 for XIST and Median X vs Median XIST') +
  xlab('Median XIST') +
  ylab('R2 for XIST and Median X') +
  xlim(0,150) +
  ylim(0,1) 
ggsave(XIST.R2_SCATTER, device="tiff")

# ______________________________________________________________________________________________________________________
#  Prepare dfs of median X and XIST expression in f + m brain tissues, sex-shared tissues, 
#  and sex-specific tissues for violin plots 
# ______________________________________________________________________________________________________________________
# Prepare df
# Collapse list of dfs into one df of just the median X counts
f.df <- ldply(f.MedX_Vs_XIST, data.frame)
f.df$Tissue <- f.df$.id
f.df$.id <- NULL
f.df$Sex <- 'Female'
head(f.df)

m.df <- ldply(m.MedX_Vs_XIST, data.frame)
m.df$Tissue <- m.df$.id
m.df$.id <- NULL
m.df$Sex <- 'Male'
head(m.df)

# Seperate into brain, and two dfs non-brain tissues
f.Not_Brain <- setdiff(names(lm_f.MedX_XIST), Brain_Tissues)
m.Not_Brain <- setdiff(names(lm_m.MedX_XIST), Brain_Tissues)

# Subset median X dfs for plots
Brain_f.df <- f.df[f.df$Tissue %in% Brain_Tissues,]
Not_Brain_f.df <- f.df[f.df$Tissue %in% f.Not_Brain,]

Brain_m.df <- m.df[m.df$Tissue %in% Brain_Tissues,]
Not_Brain_m.df <- m.df[m.df$Tissue %in% m.Not_Brain,]

# Combine median X dfs into one df
# Brain tissues only
Brain.df <- rbind(Brain_f.df, Brain_m.df)
head(Brain.df); tail(Brain.df)
nrow(Brain.df) == nrow(Brain_f.df) + nrow(Brain_m.df)

# Sex-specific tissue types list
f.Only <- setdiff(f.Not_Brain, m.Not_Brain)
m.Only <- setdiff(m.Not_Brain, f.Not_Brain)

# Subset median dfs and combine to get df of one sex, not brain tissues
Only_f.df <- f.df[f.df$Tissue %in% f.Only,]
Only_m.df <- m.df[m.df$Tissue %in% m.Only,]

# Add col indicating sex to median X and XIST dfs
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
#  Median X Split-violin plots 
# ______________________________________________________________________________________________________________________
# ylims
# Set to 90 TPM to keep yaxis proportional between plots
max(Brain.df$MedX) # 16.81
max(Shared.df$MedX) # 17.38

Split_Violin <- function(DF, TITLE){
   ggplot(DF, aes(x = Tissue, y = MedX, color = Sex, fill=Sex)) +
      geom_split_violin(color='black') +
      scale_fill_manual(values=c('blue', 'green')) +
      ylim(c(0,90)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(TITLE)
}
# Split violin plot of brain tissues
Split_Violin(DF=Brain.df, TITLE="Split violin plot of median X expression in brain tissues")
ggsave(BRAIN_VIOLIN, device="tiff")

# Split violin plot of sex-shared tissues
Split_Violin(DF=Shared.df, TITLE="Split violin plot of median X expression in sex-shared tissues")
ggsave(SHARED_VIOLIN, device="tiff")

# ______________________________________________________________________________________________________________________
#  Sex-specific Tissues Violin Plot
# ______________________________________________________________________________________________________________________
# ylim
max(Sex_Specific.df$MedX) # 16.525
levels(as.factor(Sex_Specific.df$Tissue))

# Plot
# Explicitely set the order of samples; default is alphabetical
ggplot(Sex_Specific.df, 
       aes(x = factor(Tissue, level=c("Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", 
                                      "Ovary", "Uterus", "Vagina", "Prostate", "Testis")),
           y = MedX,
           color=Sex,
           fill=Sex)) +
   geom_violin(color='black') +
   scale_fill_manual(values=c('blue', 'green')) +
   ylim(c(0,80)) + 
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
   ggtitle("Violin plot of median X expression in sex-specific tissues")
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
               y = XIST,
               color = Tissue,
               fill = Tissue)) +
        geom_violin(color='black') +
        ylim(YLIM) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(TITLE) +
        guides(fill=FALSE)
}
XIST_Violin(DF=Brain_f.df, 
            YLIM=c(0, ceiling(max(Brain_f.df$XIST))), 
            TITLE="Violin plot of log10 XIST expression across female brain tissues")
ggsave(F.BRAIN_XIST, device="tiff")
XIST_Violin(DF=Not_Brain_f.df, 
            YLIM=c(0, ceiling(max(Not_Brain_f.df$XIST))), 
            TITLE="Violin plot of log10 XIST expression across female tissues")
ggsave(F.NOT_BRAIN_XIST, device="tiff")
XIST_Violin(DF=Brain_m.df, 
            YLIM=c(0, ceiling(max(Brain_m.df$XIST))), 
            TITLE="Violin plot of XIST expression across male brain tissues")
ggsave(M.BRAIN_XIST, device="tiff")
XIST_Violin(DF=Not_Brain_m.df, 
            YLIM=c(0, ceiling(max(Not_Brain_m.df$XIST))), 
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

Columns_Lst <- c('R2_Silenced_Median', 'R2_Tuk_Silenced_Median', 'R2_Bal_Silenced_Median', 'R2_One_Silenced_Median',
                 'R2_Immune_Silenced_Median', 'R2_One_Variable_Median', 'R2_Tuk_Variable_Median', 'R2_Bal_Variable_Median', 
                 'R2_Immune_Variable_Median', 'R2_One_Incomplete_Median', 'R2_Tuk_Incomplete_Median', 'R2_Bal_Incomplete_Median', 
                 'R2_Immune_Incomplete_Median', 'R2_All_Eval', 'R2_Not_Eval', 'R2_Immune_Not_Eval', 'R2_PAR')

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
ggsave(R2_VIOLIN, width=5, height=5)

# ______________________________________________________________________________________________________________________
#  Venn diagram; overlap between gene classification schemes
# ______________________________________________________________________________________________________________________
summary(Gene_Lst)

# Plot
ggVennDiagram(Gene_Lst[c('Incomplete_In_Tukiainen', 'Incomplete_In_Balaton')]) +
    scale_fill_gradient(low='blue', high='red') +
    ggtitle('Venn Diagram of Incomplete Genes')
ggsave(VENN_INCOMPLETE, device="tiff")
ggVennDiagram(Gene_Lst[c('Variable_In_Tukiainen', 'Variable_In_Balaton')]) +
    scale_fill_gradient(low='blue', high='red') +
    ggtitle('Venn Diagram of Variable Genes')
ggsave(VENN_VARIABLE, device="tiff")
ggVennDiagram(Gene_Lst[c('Silenced_In_Tukiainen', 'Silenced_In_Balaton')]) +
    scale_fill_gradient(low='blue', high='red') +
    ggtitle('Venn Diagram of Silenced Genes')
ggsave(VENN_SILENCED, device="tiff")

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
ggsave(INTER_VENN)

