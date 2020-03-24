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
load('DDX3X_Tissue_092319.RData')

# _________________________________________________________________________________________________________________________________
#  Scatter Plots
# _________________________________________________________________________________________________________________________________
# Find the max x/y values to set as x/y lim
Max_Func <- function(x, y){
  lim <- max(x[[y]])
  return(lim)
}

# x limits
f.xmax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_DDX3X, y='DDX3X'))))
m.xmax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_DDX3X, y='DDX3X')))
# y limits
f.ymax <- round(max(unlist(Map(Max_Func, x=f.MeanX_Vs_DDX3X, y='MeanX'))))
m.ymax <- max(unlist(Map(Max_Func, x=m.MeanX_Vs_DDX3X, y='MeanX')))

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, XMAX, YMAX, RESULTS){
  plot(LM$model$DDX3X, LM$model$MeanX, main=TITLE, xlab='DDX3X', ylab='Mean X chromosome',
       xlim=c(0, XMAX), ylim=c(0, YMAX))
  legend("bottomright", bty="n", legend=paste("R^2: ", format(RESULTS$r_2, digits=3), "; p_Val: ", format(RESULTS$p_val, digits=3)))
  abline(LM)
}

# Print plots
pdf('DDX3X_Fem_Tissue_Scatter_Plots.pdf')
f.Scatter <- Map(Scatter_Func, LM=lm_f.MeanX_DDX3X, TITLE=names(lm_f.MeanX_DDX3X), XMAX=f.xmax, YMAX=f.ymax, RESULTS=Res_f.MeanX_DDX3X)
dev.off()

# set x lim to male max(DDX3X) but keep y lim as female max(MeanX)
pdf('DDX3X_Male_Tissue_Scatter_Plots.pdf')
m.Scatter <- Map(Scatter_Func, LM=lm_m.MeanX_DDX3X, TITLE=names(lm_m.MeanX_DDX3X), XMAX=m.xmax, YMAX=f.ymax, RESULTS=Res_m.MeanX_DDX3X)
dev.off()

# _________________________________________________________________________________________________________________________________
#  Scatter plot of Mean X vs DDX3X R2 values across female tissues
# _________________________________________________________________________________________________________________________________
# Drop transformed cell lines (known to re-activate X)
Drop <- c('Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts')
which(rownames(f.Regression) %in% Drop) # 13, 35

f.Regression <- f.Regression[-c(13,35), ]

ggplot(f.Regression, aes(x=Tissue, y=R2_MeanX)) +
  geom_point() +
  ggtitle('Scatter plot of Mean X vs DDX3X R2 in Female Tissues') +
  xlab('Tissue Type') +
  ylab('R^2 Mean X vs DDX3X') +
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# _________________________________________________________________________________________________________________________________
#  Scatter plot of R^2 of MeanX and DDX3X vs DDX3X
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

# df of R2_MeanX and Mean_DDX3X for both females and males
Subset_f.df <- f.Regression[ ,c('Tissue', 'R2_MeanX', 'Mean_DDX3X')]
Common_f.df <- Subset_f.df[Subset_f.df$Tissue %in% Shared, ]
Common_f.df$Sex <- 'Female'

# Add values from males
Subset_m.df <- m.Regression[ ,c('Tissue', 'R2_MeanX', 'Mean_DDX3X')]
Common_m.df <- Subset_m.df[Subset_m.df$Tissue %in% Shared, ]
Common_m.df$Sex <- 'Male'

# Combine 
Common.df <- rbind(Common_f.df, Common_m.df)

# x and y lims
max(Common.df$Mean_DDX3X) # 192.3521
max(Common.df$R2_MeanX) # 0.8966322

# Male outlier- unusually high R2 
m.Regression[which.max(m.Regression$R2_MeanX), ] # Kidney - Cortex

# Scatter plot
pdf('R2_DDX3X_and_MeanX_Vs_MeanDDX3X.pdf')
ggplot(Common.df, aes(x=Mean_DDX3X, y=R2_MeanX)) +
  geom_point(aes(shape=Sex, fill=Sex), size=2) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=c('blue', 'green')) +
  ggtitle('R^2 for DDX3X and Mean X vs Mean DDX3X') +
  xlab('Mean DDX3X') +
  ylab('R2 for DDX3X and Mean X') +
  xlim(0,150) +
  ylim(0,1)
dev.off()

# _________________________________________________________________________________________________________________________________
#  Single-Sex Violin Plots 
# _________________________________________________________________________________________________________________________________
# For future reference, don't use plotly; You have to pay for a subscription to use pdf()/tiff() with plotly objects 
# Printed plots using viewer

# Collapse list of dfs into one with just DDX3X column
f.df <- ldply(f.MeanX_Vs_DDX3X, data.frame)
f.df$Tissue <- f.df$.id
f.df$.id <- NULL

m.df <- ldply(m.MeanX_Vs_DDX3X, data.frame)
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
      y = ~DDX3X,
      split = ~Tissue,
      type = 'violin',
      box = list(visible = T),
      meanline = list(visible = T)
    ) %>% 
    layout(
      title = TITLE,
      showlegend = FALSE,
      xaxis = list(title = "Tissue Type"),
      yaxis = list(title = "DDX3X", zeroline = F, range = RANGE)
    )
  return(plot)
}
Violin_Func(DF=Brain_f.df, RANGE=c(0,400), TITLE='Violin Plots of DDX3X Expression in Female Brain Tissues')
Violin_Func(DF=Brain_m.df, RANGE=c(0,400), TITLE='Violin Plots of DDX3X Expression in Male Brain Tissues')
Violin_Func(DF=Not_Brain_m.df, RANGE=c(0,400), TITLE='Violin Plots of DDX3X Expression in Male Tissues')
Violin_Func(DF=Not_Brain_f.df, RANGE=c(0,400), TITLE='Violin Plots of DDX3X Expression in Female Tissues')

# _________________________________________________________________________________________________________________________________
#  Grouped Violin Plot
# _________________________________________________________________________________________________________________________________
# Sex-specific tissue types list
f.Only <- setdiff(f.Not_Brain, m.Not_Brain)
m.Only <- setdiff(m.Not_Brain, f.Not_Brain)

# Subset dfs and combine
Only_f.df <- f.df[f.df$Tissue %in% f.Only,]
Only_m.df <- m.df[m.df$Tissue %in% m.Only,]

# Add col indicating sex
Only_f.df$Sex <- 'Female'
Only_m.df$Sex <- 'Male'

# Join
Sex_Specific <- rbind(Only_f.df, Only_m.df)

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
    title = 'Violin Plot of Mean X Expression in Sex-Specific Tissues',
    xaxis = list(title = "Tissue Types"),
    yaxis = list(title = "Mean X Chm Read Count", zeroline = F, range=c(0,25))
  ) 
p

# _________________________________________________________________________________________________________________________________
#  Split Violin Plots 
# _________________________________________________________________________________________________________________________________
# Brain; Female/Male split violin plots
# Subset female and male dfs to only sex-shared tissues
Shared_f.df <- f.df[f.df$Tissue %in% Shared, ]
Shared_m.df <- m.df[m.df$Tissue %in% Shared, ]

# Check that both dfs have same factor levels is Tissue column
levels(factor(Shared_m.df$Tissue)) == levels(factor(Shared_f.df$Tissue)) # all TRUE

# Add column indicating sex 
Shared_f.df$Sex <- 'Female' # nrow = 3704
Shared_m.df$Sex <- 'Male'   # nrow = 7040

# Merge female/male dfs containing shared tissues
all.df <- rbind(Shared_f.df, Shared_m.df) # nrow = 10744; expected 
Brain_all.df <- all.df[all.df$Tissue %in% Brain_Tissues,]

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

# Subset df for plots
Not_Brain_all.df <- all.df[all.df$Tissue %in% all.Not_Brain,]

# ylim
max(Not_Brain_all.df$MeanX) # 21.4484

# Split violin function
Split_Violin <- function(DF, TITLE){
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
      yaxis = list(title = "Mean X Chm Read Count", zeroline = F, range=c(0,25))
    ) 
  return(plot)
}
Split_Violin(DF=Not_Brain_all.df, TITLE="Mean X Chromosome Expression in Tissues Common to Both Sexes")
Split_Violin(DF=Brain_all.df, TITLE="Mean X Chromosome Expression in Female and Male Brain Tissues")

# _________________________________________________________________________________________________________________________________
#  Venn diagram; overlap between gene classification schemes
# _________________________________________________________________________________________________________________________________
# Venn diagram function for 2 groups
Venn_2 <- function(a, b, TITLE){
  venn.plot <- venn.diagram(
    x = list('Tukainen' = Gene_Lst[[a]], 'Balaton'= Gene_Lst[[b]]),
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
pdf('RD_Venn_GeneCategories.pdf')
Map(Venn_2, a='Incomplete_In_Tukiainen', b='Incomplete_In_Balaton', TITLE='Venn Diagram of Incomplete Genes')
Map(Venn_2, a='Variable_In_Tukiainen', b='Variable_In_Balaton', TITLE='Venn Diagram of Variable Genes')
Map(Venn_2, a='Silenced_In_Tukiainen', b='Silenced_In_Balaton', TITLE='Venn Diagram of Silenced Genes')
dev.off()

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






