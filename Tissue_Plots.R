#!/usr/bin/env Rscript

# Plots related to tissue linear model results
setwd("~/XIST/")

# Constants
DATA <- "~/XIST/Tissue/Xpressed/Mean/Xpressed_012720.RData"
F.QQ <- "~/XIST/Tissue/Xpressed/Mean/Fem_QQ_Plots.pdf"
M.QQ <- "~/XIST/Tissue/Xpressed/Mean/Male_QQ_Plots.pdf"
F.SCATTER <- "~/XIST/Tissue/Xpressed/Mean/Fem_Tissue_Scatter_Plots.pdf"
M.SCATTER <- "~/XIST/Tissue/Xpressed/Mean/Male_Tissue_Scatter_Plots.pdf"
R2_SCATTER <- "~/XIST/Tissue/Xpressed/Mean/R2_Scatter.pdf"
F.BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Mean/Fem_Brain_XIST.pdf"
M.BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Mean/Male_Brain_XIST.pdf"
F.NOT_BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Mean/Fem_NotBrain_XIST.pdf"
M.NOT_BRAIN_XIST <- "~/XIST/Tissue/Xpressed/Mean/Male_NotBrain_XIST.pdf"
SEX_SPECIFIC <- "~/XIST/Tissue/Xpressed/Mean/Sex_Specific_Violin.pdf"
NOT_BRAIN_VIOLIN <- "~/XIST/Tissue/Xpressed/Mean/NotBrain_Violin.pdf"
BRAIN_VIOLIN <- "~/XIST/Tissue/Xpressed/Mean/Brain_Violin.pdf"
XIST.R2_SCATTER <- "~/XIST/Tissue/Xpressed/Mean/R2_vs_MeanXIST.pdf"
VENN <- "~/XIST/Tissue/Xpressed/Mean/Venn_GeneCategories.pdf"
R2_VIOLIN <- "~/XIST/Tissue/Xpressed/Mean/R2_Violin.pdf"

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
library(VennDiagram)
library(grDevices)
library(grid)

# Load session data
load(DATA)

# ______________________________________________________________________________________________________________________
#  QQ Plots
# ______________________________________________________________________________________________________________________
# Functon to plot LM residuals
QQ_Func <- function(LM, SEX, TISSUE){
  plot <- qqnorm(y = LM[['model']][['XIST']], # predictor/independent var
                 ylab = "Sample XIST TPMs",
                 xlab = "Theoretical Quantiles",
                 main = paste("lm(MeanX ~ XIST) in", SEX, TISSUE))
  qqline(LM[['model']][['XIST']])
  return(plot)
}

pdf(F.QQ)
Map(QQ_Func, LM=lm_f.MeanX_XIST, SEX=c("Female"), TISSUE=names(lm_f.MeanX_XIST))
dev.off()
pdf(M.QQ)
Map(QQ_Func, LM=lm_m.MeanX_XIST, SEX=c("Male"), TISSUE=names(lm_m.MeanX_XIST))
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
f.xmax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_XIST, y='XIST'))))
m.xmax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_XIST, y='XIST')))
# y limits
f.ymax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_XIST, y='MeanX'))))
m.ymax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_XIST, y='MeanX')))

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, XMAX, YMAX, R2, PVAL){
  plot(LM[['model']][['XIST']], 
       LM[['model']][['MeanX']], 
       main=TITLE, 
       xlab='XIST', 
       ylab='Mean X chromosome',
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
                 LM=lm_f.MeanX_XIST, 
                 TITLE=names(lm_f.MeanX_XIST), 
                 XMAX=f.xmax, 
                 YMAX=f.ymax, 
                 R2=f.Regression[['R2_MeanX']],
                 PVAL=f.Regression[['pval_MeanX']])
dev.off()

# set x lim to male max(XIST) but keep y lim as female max(MeanX)
pdf(M.SCATTER)
m.Scatter <- Map(Scatter_Func, 
                 LM=lm_m.MeanX_XIST, 
                 TITLE=names(lm_m.MeanX_XIST), 
                 XMAX=m.xmax, 
                 YMAX=f.ymax, 
                 R2=m.Regression[['R2_MeanX']],
                 PVAL=m.Regression[['pval_MeanX']])
dev.off()

# ______________________________________________________________________________________________________________________
# Prepare df for plotting R2 scatter plots
# ______________________________________________________________________________________________________________________
# Combine relevant cols from dfs
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

# df of R2_MeanX and Mean_XIST for both females and males
Subset_f.df <- f.Regression[ ,c('Tissue', 'R2_MeanX', 'Mean_XIST', 'pval_MeanX')]
Common_f.df <- Subset_f.df[Subset_f.df$Tissue %in% Shared, ]
Common_f.df$Sex <- 'Female'

# Add values from males
Subset_m.df <- m.Regression[ ,c('Tissue', 'R2_MeanX', 'Mean_XIST', 'pval_MeanX')]
Common_m.df <- Subset_m.df[Subset_m.df$Tissue %in% Shared, ]
Common_m.df$Sex <- 'Male'

# Combine 
Common.df <- rbind(Common_f.df, Common_m.df)

# x and y lims
max(Common.df$Mean_XIST) # 141.7757
max(Common.df$R2_MeanX) # 0.5555296

# Add col with pval factor (p> or < 0.05)
Common.df$p_Factor <- as.factor(Common.df$pval_MeanX<0.05)
head(Common.df)


# ______________________________________________________________________________________________________________________
#  Scatter plot of Mean X vs XIST R2 values across female/male tissues
# ______________________________________________________________________________________________________________________
# Plot females and males together
pdf(R2_SCATTER)
ggplot(Common.df, aes(x=Tissue, y=R2_MeanX)) + 
 geom_point(aes(shape=Sex, fill=Sex, alpha=p_Factor)) +
 scale_shape_manual(values=c(21,22)) +
 scale_fill_manual(values=c('blue', 'green')) +
 ggtitle('Scatter plot of R2 (MeanX ~ XIST) in Female and Male Tissues') +
 xlab('Tissue Type') +
 ylab('R2 (MeanX~ XIST') +
 ylim(c(0,1)) +
 theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# If I want to keep fe/male plots seperate
#ggplot(m.Regression, aes(x=Tissue, y=R2_MeanX)) + # refer directly to df columns
#  geom_point() +
#  ggtitle(paste('Scatter plot of Mean X vs XIST R2 in', 'Male', 'Tissues')) +
#  xlab('Tissue Type') +
#  ylab('R^2 Mean X vs XIST') +
#  ylim(c(0,1)) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ______________________________________________________________________________________________________________________
#  Scatter plot of R^2 of MeanX and XIST vs XIST
# ______________________________________________________________________________________________________________________

# Scatter plot
pdf(XIST.R2_SCATTER)
ggplot(Common.df, aes(x=Mean_XIST, y=R2_MeanX)) +
  geom_point(aes(shape=Sex, fill=Sex, alpha=p_Factor)) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=c('blue', 'green')) +
  ggtitle('R^2 for XIST and Mean X vs Mean XIST') +
  xlab('Mean XIST') +
  ylab('R2 for XIST and Mean X') +
  xlim(0,150) +
  ylim(0,1) 
dev.off()

# ______________________________________________________________________________________________________________________
#  Single-Sex Violin Plots 
# ______________________________________________________________________________________________________________________
# For future reference, don't use plotly; You have to pay for a subscription to use pdf()/tiff() with plotly objects 
# Printed plots using viewer

# Collapse list of dfs into one with just XIST column
f.df <- ldply(f.MeanX_Vs_XIST, data.frame)
f.df$Tissue <- f.df$.id
f.df$.id <- NULL

m.df <- ldply(m.MeanX_Vs_XIST, data.frame)
m.df$Tissue <- m.df$.id
m.df$.id <- NULL

# Seperate into brain, and two dfs non-brain tissues
Brain_Tissues <- c("Brain - Cortex", "Brain - Hippocampus", "Brain - Substantia nigra",                 
                   "Brain - Anterior cingulate cortex (BA24)", "Brain - Frontal Cortex (BA9)",             
                   "Brain - Cerebellar Hemisphere", "Brain - Caudate (basal ganglia)",          
                   "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)",          
                   "Brain - Hypothalamus", "Brain - Spinal cord (cervical c-1)",       
                   "Brain - Amygdala", "Brain - Cerebellum")

f.Not_Brain <- setdiff(names(lm_f.MeanX_XIST), Brain_Tissues)
m.Not_Brain <- setdiff(names(lm_m.MeanX_XIST), Brain_Tissues)

# Subset df for plots
Brain_f.df <- f.df[f.df$Tissue %in% Brain_Tissues,]
Not_Brain_f.df <- f.df[f.df$Tissue %in% f.Not_Brain,]

Brain_m.df <- m.df[m.df$Tissue %in% Brain_Tissues,]
Not_Brain_m.df <- m.df[m.df$Tissue %in% m.Not_Brain,]

# Log transform 
Brain_f.df[,1:2] <- log10(Brain_f.df[1:2])
Brain_m.df[,1:2] <- log10(Brain_m.df[1:2])
Not_Brain_f.df[,1:2] <- log10(Not_Brain_f.df[1:2])
Not_Brain_m.df[,1:2] <- log10(Not_Brain_m.df[1:2])

# Whole violin function
Violin_Func <- function(DF, RANGE, TITLE){
  plot <- DF %>%
    plot_ly(
      x = ~Tissue,
      y = ~XIST,
      split = ~Tissue,
      type = 'violin',
      box = list(visible = T),
      meanline = list(visible = T)
    ) %>% 
    layout(
      title = TITLE,
      showlegend = FALSE,
      xaxis = list(title = "Tissue Type"),
      yaxis = list(title = "XIST", zeroline = F, range = RANGE)
    )
  return(plot)
}
pdf(F.BRAIN_XIST)
Violin_Func(DF=Brain_f.df, RANGE=c(0, 5), TITLE='Violin Plots of log10 XIST Expression in Female Brain Tissues')
dev.off()
pdf(M.BRAIN_XIST)
Violin_Func(DF=Brain_m.df, RANGE=c(0,5), TITLE='Violin Plots of log10 XIST Expression in Male Brain Tissues')
dev.off()
pdf(F.NOT_BRAIN_XIST)
Violin_Func(DF=Not_Brain_m.df, RANGE=c(0,5), TITLE='Violin Plots of log10 XIST Expression in Male Tissues')
dev.off()
pdf(M.NOT_BRAIN_XIST)
Violin_Func(DF=Not_Brain_f.df, RANGE=c(0, 5), TITLE='Violin Plots of log10 XIST Expression in Female Tissues')
dev.off()

# ______________________________________________________________________________________________________________________
#  Grouped Violin Plot
# ______________________________________________________________________________________________________________________
# Sex-specific tissue types list
f.Only <- setdiff(f.Not_Brain, m.Not_Brain)
m.Only <- setdiff(m.Not_Brain, f.Not_Brain)

# Subset dfs and combine to get df of one sex, not brain tissues
Only_f.df <- f.df[f.df$Tissue %in% f.Only,]
Only_m.df <- m.df[m.df$Tissue %in% m.Only,]

# Add col indicating sex
Only_f.df$Sex <- 'Female'
Only_m.df$Sex <- 'Male'

# Join
Sex_Specific <- rbind(Only_f.df, Only_m.df)

# Log transform
Sex_Specific[,1:2] <- log10(Sex_Specific[1:2])

# Plot
p <- Sex_Specific %>%
  plot_ly(type = 'violin') %>%
  add_trace(
    x = ~Tissue[Sex_Specific$Sex == 'Female'],
    y = ~MeanX[Sex_Specific$Sex == 'Female'],
    legendgroup = 'Female',
    scalegroup = 'Female',
    name = 'Female',
    box = list(visible = T),
    meanline = list(visible = T),
    line = list(color = 'darkblue'),
    fillcolor = 'lightblue',
    marker = list(color='darkblue')
  ) %>%
  add_trace(
    x = ~Tissue[Sex_Specific$Sex == 'Male'],
    y = ~MeanX[Sex_Specific$Sex == 'Male'],
    legendgroup = 'Male',
    scalegroup = 'Male',
    name = 'Male',
    box = list(visible = T),
    meanline = list(visible = T),
    line = list(color = 'darkgreen'),
    fillcolor = 'lightgreen'
  ) %>% 
  layout(
    title = 'Violin Plot of log10 Mean X Expression in Sex-Specific Tissues',
    xaxis = list(title = "Tissue Types"),
    yaxis = list(title = "Mean X Chm Read Count", zeroline = F, range=c(0,5))
  ) 
pdf(SEX_SPECIFIC)
p
dev.off()

# ______________________________________________________________________________________________________________________
#  Split Violin Plots 
# ______________________________________________________________________________________________________________________
# Brain; Female/Male split violin plots
# Subset female and male dfs to only sex-shared tissues
Shared_f.df <- f.df[f.df$Tissue %in% Shared, ]
Shared_m.df <- m.df[m.df$Tissue %in% Shared, ]

# Check that both dfs have same factor levels is Tissue column
all(levels(factor(Shared_m.df$Tissue)) == levels(factor(Shared_f.df$Tissue))) # TRUE

# Add column indicating sex 
Shared_f.df$Sex <- 'Female' # nrow = 3704
Shared_m.df$Sex <- 'Male'   # nrow = 7040

# Merge female/male dfs containing shared tissues
Shared.df <- rbind(Shared_f.df, Shared_m.df) # nrow = 10744; expected 
Brain_Shared.df <- Shared.df[Shared.df$Tissue %in% Brain_Tissues,]

# Not-Brain tissues shared between sexes; Split violin plot
all.Not_Brain <- c("Adipose - Subcutaneous", "Muscle - Skeletal","Artery - Tibial", "Artery - Coronary",                        
                   "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)", "Breast - Mammary Tissue", 
                   "Skin - Not Sun Exposed (Suprapubic)", "Minor Salivary Gland", "Adrenal Gland", "Thyroid",                                  
                   "Lung", "Spleen", "Pancreas", "Esophagus - Muscularis",                   
                   "Esophagus - Mucosa", "Esophagus - Gastroesophageal Junction",
                   "Stomach", "Colon - Sigmoid","Small Intestine - Terminal Ileum", "Colon - Transverse",                       
                   "Skin - Sun Exposed (Lower leg)", "Nerve - Tibial",                           
                   "Heart - Left Ventricle", "Pituitary", "Cells - Transformed fibroblasts",          
                   "Whole Blood", "Artery - Aorta","Cells - EBV-transformed lymphocytes", "Liver",                                    
                   "Kidney - Cortex", "Bladder")    

all.Not_Brain <- setdiff(Shared, Brain_Tissues)

# Subset df for plots
Not_Brain_Shared.df <- Shared.df[Shared.df$Tissue %in% all.Not_Brain,]

# log transform 
Not_Brain_Shared.df[,1:2] <- log10(Not_Brain_Shared.df[1:2])
Brain_Shared.df[,1:2] <- log10(Brain_Shared.df[1:2])

# ylim
max(Not_Brain_Shared.df$MeanX) # untransformed = 21.47722, log10 = 1.331978

# Split violin function
Split_Violin <- function(DF, TITLE, RANGE){
  plot <- plot_ly(data = DF, type = 'violin') %>%
    add_trace(
      x = ~Tissue[DF$Sex == 'Female'],
      y = ~MeanX[DF$Sex == 'Female'],
      legendgroup = 'Female',
      scalegroup = 'Female',
      name = 'Female',
      side = 'negative',
      box = list(visible = T),
      meanline = list(visible = T),
      line = list(color = 'darkblue'),
      fillcolor = 'darkblue',
      marker = list(color='darkblue')
    ) %>%
    add_trace(
      x = ~Tissue[DF$Sex == 'Male'],
      y = ~MeanX[DF$Sex == 'Male'],
      legendgroup = 'Male',
      scalegroup = 'Male',
      name = 'Male',
      side = 'positive',
      box = list(visible = T),
      meanline = list(visible = T),
      line = list(color = 'darkgreen'),
      fillcolor = 'lightgreen',
      marker = list(color='darkgreen')
    ) %>% 
    layout(
      title = TITLE,
      xaxis = list(title = "Tissue Types"),
      yaxis = list(title = "Mean X Chm Read Count", zeroline = F, range=RANGE)
    ) 
  return(plot)
}
pdf(NOT_BRAIN_VIOLIN)
Split_Violin(DF=Not_Brain_Shared.df, 
             TITLE="log10 Mean X Chromosome Expression in Tissues Common to Both Sexes",
             RANGE=c(0,5))
dev.off()
pdf(BRAIN_VIOLIN)
Split_Violin(DF=Brain_Shared.df, 
             TITLE="log10 Mean X Chromosome Expression in Female and Male Brain Tissues",
             RANGE=c(0,5))
dev.off()

# ______________________________________________________________________________________________________________________
#  Venn diagram; overlap between gene classification schemes
# ______________________________________________________________________________________________________________________
# Venn diagram function for 2 groups
Venn_2 <- function(a, b, TITLE){
  venn.plot <- venn.diagram(
    x = list('Tukiainen' = Gene_Lst[[a]], 'Balaton'= Gene_Lst[[b]]),
    label=TRUE,
    filename = NULL,
    scaled = TRUE,
    col = "transparent",
    fill = c("darkblue", "darkgreen"),
    main.pos = c(0.5, 1.0),
    cex = 1.5,
    cat.cex = 1.5,
    main.cex = 2,
    cat.default.pos = "outer",
    cat.just= list(c(0, 0.5), c(0, 0.5)), 
    cat.dist = c(0.01, 0.01),
    cat.fontfamily = "sans",
    main = TITLE,
    fontfamily = "sans",
    na = "remove",
    inverted = FALSE)
  grid.newpage()
  return(grid.draw(venn.plot))
}
pdf(VENN)
Venn_2(a='Incomplete_In_Tukiainen', b='Incomplete_In_Balaton', TITLE='Venn Diagram of Incomplete Genes')
Venn_2(a='Variable_In_Tukiainen', b='Variable_In_Balaton', TITLE='Venn Diagram of Variable Genes')
Venn_2(a='Silenced_In_Tukiainen', b='Silenced_In_Balaton', TITLE='Venn Diagram of Silenced Genes')
dev.off()

# ______________________________________________________________________________________________________________________
#  Violin plot: R2 by gene classification scheme per tissue
# ______________________________________________________________________________________________________________________
# Reshape correlation df 
Category_Lst <- c('Silenced In Both', 'Silenced In Tukainen', 'Silenced In Balaton', 'Silenced In At Least One', 
                  'Silenced Immune Genes', 'Variable In At Least One', 'Variable In Tukainen',
                  'Variable In Balaton', 'Variable Immune Genes', 'Incomplete In At Least One', 'Incomplete in Tukainen',
                  'Invariable In Balaton', 'Invariable Immune Genes', 'All Evaluated In Both', 'Not Evaluated In Either', 
                  'Immune Genes Not Evaluated', 'PAR Genes In Balaton')

Columns_Lst <- c('R2_Silenced_Mean', 'R2_Tuk_Silenced_Mean', 'R2_Bal_Silenced_Mean', 'R2_One_Silenced_Mean',
                 'R2_Immune_Silenced_Mean', 'R2_One_Variable_Mean', 'R2_Tuk_Variable_Mean', 'R2_Bal_Variable_Mean', 
                 'R2_Immune_Variable_Mean', 'R2_One_Incomplete_Mean', 'R2_Tuk_Incomplete_Mean', 'R2_Bal_Incomplete_Mean', 
                 'R2_Immune_Incomplete_Mean', 'R2_All_Eval', 'R2_Not_Eval', 'R2_Immune_Not_Eval', 'R2_PAR')

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

# Grouped violin plot
plot <- plot_ly(data = Correlations.df, type = 'violin') %>%
    add_trace(
      x = ~Category[Correlations.df$Sex == 'Female'],
      y = ~R2[Correlations.df$Sex == 'Female'],
      legendgroup = 'Female',
      scalegroup = 'Female',
      name = 'Female',
      side = 'negative',
      box = list(visible = T),
      meanline = list(visible = T),
      line = list(color = 'darkblue'),
      fillcolor = 'darkblue',
      marker = list(color='darkblue')
    ) %>%
    add_trace(
      x = ~Category[Correlations.df$Sex == 'Male'],
      y = ~R2[Correlations.df$Sex == 'Male'],
      legendgroup = 'Male',
      scalegroup = 'Male',
      name = 'Male',
      side = 'positive',
      box = list(visible = T),
      meanline = list(visible = T),
      line = list(color = 'darkgreen'),
      fillcolor = 'lightgreen',
      marker = list(color='darkgreen')
    ) %>% 
    layout(
      title = 'Violin Plot of R^2 per Gene Subset',
      xaxis = list(title = "Gene Subset"),
      yaxis = list(title = "R^2", zeroline = F, range=c(0,1))
    ) 
pdf(R2_VIOLIN)
plot
dev.off()



