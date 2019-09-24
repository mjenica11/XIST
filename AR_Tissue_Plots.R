# Plots related to tissue linear model results
setwd("~/XIST_Vs_TSIX/Files")

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
load('AR_Tissue_092319.RData')

# _________________________________________________________________________________________________________________________________
#  Scatter Plots
# _________________________________________________________________________________________________________________________________
# Find the max x/y values to set as x/y lim
Max_Func <- function(x, y){
  lim <- max(x[[y]])
  return(lim)
}

# x limits
f.xmax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_AR, y='AR'))))
m.xmax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_AR, y='AR')))
# y limits
f.ymax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_AR, y='MeanX'))))
m.ymax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_AR, y='MeanX')))

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, XMAX, YMAX, RESULTS){
  plot(LM$model$AR, LM$model$MeanX, main=TITLE, xlab='AR', ylab='Mean X chromosome',
       xlim=c(0, XMAX), ylim=c(0, YMAX))
  legend("bottomright", bty="n", legend=paste("R^2: ", format(RESULTS$r_2, digits=3), "; p_Val: ", format(RESULTS$p_val, digits=3)))
  abline(LM)
}

# Print plots
pdf('~/XIST_Vs_TSIX/Files/AR/AR_Fem_Tissue_Scatter_Plots.pdf')
f.Scatter <- Map(Scatter_Func, LM=lm_f.MeanX_AR, TITLE=names(lm_f.MeanX_AR), XMAX=f.xmax, YMAX=f.ymax, RESULTS=Res_f.MeanX_AR)
dev.off()

# set x lim to male max(AR) but keep y lim as female max(MeanX)
pdf('~/XIST_Vs_TSIX/Files/AR/AR_Male_Tissue_Scatter_Plots.pdf')
m.Scatter <- Map(Scatter_Func, LM=lm_m.MeanX_AR, TITLE=names(lm_m.MeanX_AR), XMAX=m.xmax, YMAX=f.ymax, RESULTS=Res_m.MeanX_AR)
dev.off()

# _________________________________________________________________________________________________________________________________
#  Scatter plot of Mean X vs AR R2 values across female tissues
# _________________________________________________________________________________________________________________________________
# Drop transformed cell lines (known to re-activate X)
Drop <- c('Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts')
which(rownames(f.Regression) %in% Drop) # 13, 35

f.Regression <- f.Regression[-c(13,35), ]

ggplot(f.Regression, aes(x=Tissue, y=R2_MeanX)) +
  geom_point() +
  ggtitle('Scatter plot of Mean X vs AR R2 in Female Tissues') +
  xlab('Tissue Type') +
  ylab('R^2 Mean X vs AR') +
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# _________________________________________________________________________________________________________________________________
#  Scatter plot of R^2 of MeanX and AR vs AR
# _________________________________________________________________________________________________________________________________
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

# df of R2_MeanX and Mean_AR for both females and males
Subset_f.df <- f.Regression[ ,c('Tissue', 'R2_MeanX', 'Mean_AR')]
Common_f.df <- Subset_f.df[Subset_f.df$Tissue %in% Shared, ]
Common_f.df$Sex <- 'Female'

# Add values from males
Subset_m.df <- m.Regression[ ,c('Tissue', 'R2_MeanX', 'Mean_AR')]
Common_m.df <- Subset_m.df[Subset_m.df$Tissue %in% Shared, ]
Common_m.df$Sex <- 'Male'

# Combine 
Common.df <- rbind(Common_f.df, Common_m.df)

# x and y lims
max(Common.df$Mean_AR) # 18.31298
max(Common.df$R2_MeanX) # 0.7070481

# Highest R2 in male tissues 
m.Regression[which.max(m.Regression$R2_MeanX), ] # Testis

# Scatter plot
pdf('~/XIST_Vs_TSIX/Files/AR/R2_AR_and_MeanX_Vs_MeanAR.pdf')
ggplot(Common.df, aes(x=Mean_AR, y=R2_MeanX)) +
  geom_point(aes(shape=Sex, fill=Sex), size=2) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=c('blue', 'green')) +
  ggtitle('R^2 for AR ~ Mean X vs Mean AR') +
  xlab('Mean AR') +
  ylab('R2 for AR and Mean X') +
  xlim(0,150) +
  ylim(0,1)
dev.off()

# _________________________________________________________________________________________________________________________________
#  Single-Sex Violin Plots 
# _________________________________________________________________________________________________________________________________
# For future reference, don't use plotly; You have to pay for a subscription to use pdf()/tiff() with plotly objects 
# Printed plots using viewer

# Collapse list of dfs into one with just AR column
f.df <- ldply(f.MeanX_Vs_AR, data.frame)
f.df$Tissue <- f.df$.id
f.df$.id <- NULL

m.df <- ldply(m.MeanX_Vs_AR, data.frame)
m.df$Tissue <- m.df$.id
m.df$.id <- NULL

# Seperate into brain, and two dfs non-brain tissues
Brain_Tissues <- c("Brain - Cortex", "Brain - Hippocampus", "Brain - Substantia nigra",                 
                   "Brain - Anterior cingulate cortex (BA24)", "Brain - Frontal Cortex (BA9)",             
                   "Brain - Cerebellar Hemisphere", "Brain - Caudate (basal ganglia)",          
                   "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)",          
                   "Brain - Hypothalamus", "Brain - Spinal cord (cervical c-1)",       
                   "Brain - Amygdala", "Brain - Cerebellum")

f.Not_Brain <- c("Adipose - Subcutaneous", "Muscle - Skeletal", "Artery - Tibial", "Artery - Coronary",                        
                 "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)", "Ovary", "Uterus",                                   
                 "Vagina", "Breast - Mammary Tissue", "Skin - Not Sun Exposed (Suprapubic)", 
                 "Minor Salivary Gland", "Adrenal Gland", "Thyroid", "Lung","Spleen", "Pancreas",                                  
                 "Esophagus - Muscularis", "Esophagus - Mucosa","Esophagus - Gastroesophageal Junction",
                 "Stomach", "Colon - Sigmoid", "Small Intestine - Terminal Ileum",         
                 "Colon - Transverse", "Skin - Sun Exposed (Lower leg)",           
                 "Nerve - Tibial", "Heart - Left Ventricle", "Pituitary",                                                        
                 "Cells - Transformed fibroblasts", "Whole Blood",                              
                 "Artery - Aorta", "Cells - EBV-transformed lymphocytes",       
                 "Liver", "Kidney - Cortex", "Fallopian Tube",                           
                 "Bladder", "Cervix - Ectocervix", "Cervix - Endocervix")

m.Not_Brain <- c("Adipose - Subcutaneous", "Muscle - Skeletal","Artery - Tibial", "Artery - Coronary",                        
                 "Heart - Atrial Appendage", "Adipose - Visceral (Omentum)", "Breast - Mammary Tissue", 
                 "Skin - Not Sun Exposed (Suprapubic)", "Minor Salivary Gland", "Adrenal Gland", "Thyroid",                                  
                 "Lung", "Spleen", "Pancreas", "Esophagus - Muscularis",                   
                 "Esophagus - Mucosa", "Esophagus - Gastroesophageal Junction",
                 "Stomach", "Colon - Sigmoid","Small Intestine - Terminal Ileum", "Colon - Transverse",                       
                 "Prostate", "Testis", "Skin - Sun Exposed (Lower leg)", "Nerve - Tibial",                           
                 "Heart - Left Ventricle", "Pituitary", "Cells - Transformed fibroblasts",          
                 "Whole Blood", "Artery - Aorta","Cells - EBV-transformed lymphocytes", "Liver",                                    
                 "Kidney - Cortex", "Bladder")  

# Subset df for plots
Brain_f.df <- f.df[f.df$Tissue %in% Brain_Tissues,]
Not_Brain_f.df <- f.df[f.df$Tissue %in% f.Not_Brain,]

Brain_m.df <- m.df[m.df$Tissue %in% Brain_Tissues,]
Not_Brain_m.df <- m.df[m.df$Tissue %in% m.Not_Brain,]

# Whole violin function
Violin_Func <- function(DF, RANGE, TITLE){
  plot <- DF %>%
    plot_ly(
      x = ~Tissue,
      y = ~AR,
      split = ~Tissue,
      type = 'violin',
      box = list(visible = T),
      meanline = list(visible = T)
    ) %>% 
    layout(
      title = TITLE,
      showlegend = FALSE,
      xaxis = list(title = "Tissue Type"),
      yaxis = list(title = "AR", zeroline = F, range = RANGE)
    )
  return(plot)
}
Violin_Func(DF=Brain_f.df, RANGE=c(0,100), TITLE='Violin Plots of AR Expression in Female Brain Tissues')
Violin_Func(DF=Brain_m.df, RANGE=c(0,100), TITLE='Violin Plots of AR Expression in Male Brain Tissues')
Violin_Func(DF=Not_Brain_m.df, RANGE=c(0,100), TITLE='Violin Plots of AR Expression in Male Tissues')
Violin_Func(DF=Not_Brain_f.df, RANGE=c(0,100), TITLE='Violin Plots of AR Expression in Female Tissues')

# _________________________________________________________________________________________________________________________________
#  Violin plot: R2 by gene classification scheme per tissue
# _________________________________________________________________________________________________________________________________
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
Sub_Lst.f <- Map(Fem_Subset, SUBSET=Columns_Lst, CAT=Category_Lst)

Male_Subset <- function(SUBSET, CAT){
  res.df <- data.frame(Tissue=m.Regression[['Tissue']], R2=m.Regression[[SUBSET]], Category=CAT)
  return(res.df)
}
Sub_Lst.m <- Map(Fem_Subset, SUBSET=Columns_Lst, CAT=Category_Lst)

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
plot






