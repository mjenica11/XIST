# Plots related to tissue linear model results
setwd("~/XIST")

# Constants
SCATTER.f <- '~/XIST/Tissue/DDX3X_XIST/Fem_Scatter.tiff'
SCATTER.m <- '~/XIST/Tissue/DDX3X_XIST/Male_Scatter.tiff'
CORRELATION.f <- "~/XIST/Tissue/DDX3X_XIST/Fem_Tissue_Correlations.csv"
CORRELATION.m <- "~/XIST/Tissue/DDX3X_XIST/Male_Tissue_Correlations.csv"
SUMMARY <- "~/XIST/Tissue/DDX3X_XIST/Tissue_Regression_Averages.csv"
SLOPES <- "~/XIST/Tissue/DDX3X_XIST/Tissue_Slopes_Table.csv"

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
load('DDX3X_Gene_Tissue_011320.RData')

# ______________________________________________________________________________________________________________________
#  Correlate DDX3X ~ XIST 
# ______________________________________________________________________________________________________________________
# Get XIST values from each data frame
Get_XIST <- function(x){
  res <- as.numeric(x[which(x['gene_name'] == 'XIST'),][,-c(1:2)]) # drop gene_id and gene_name cols
  return(res)
}
f.XIST_Tissue_Counts <- lapply(f.Tissue_Counts, Get_XIST)
m.XIST_Tissue_Counts <- lapply(m.Tissue_Counts, Get_XIST)

# Cmbine vectors into df
f.DDX3X_XIST <- Combine_Lsts(x='f.DDX3X_XIST', 
                             y=f.DDX3X_Tissue_Counts, 
                             z=f.XIST_Tissue_Counts)
m.DDX3X_XIST <- Combine_Lsts(x='m.DDX3X_XIST', 
                             y=m.DDX3X_Tissue_Counts, 
                             z=m.XIST_Tissue_Counts)

# Rename col names in each df
f.DDX3X_XIST <- Map(Rename_Col, 
                    x=f.DDX3X_XIST,
                    a='DDX3X', 
                    b='XIST')
m.DDX3X_XIST <- Map(Rename_Col, 
                    x=m.DDX3X_XIST, 
                    a='DDX3X', 
                    b='XIST')

# Apply lm to each df in list
Linear_Model <- function(x) {
  z <- lm(DDX3X ~ XIST, data = x)
  return(z)
}
lm_f.DDX3X_XIST <- lapply(f.DDX3X_XIST, Linear_Model)
lm_m.DDX3X_XIST <- lapply(m.DDX3X_XIST, Linear_Model)

# Apply function to list of dfs
Res_f.DDX3X_XIST <- lapply(lm_f.DDX3X_XIST, Regression_Res)
Res_m.DDX3X_XIST <- lapply(lm_m.DDX3X_XIST, Regression_Res)

# Make table summarizing results
f.Regression <- as.data.frame(do.call(rbind, Res_f.DDX3X_XIST))
m.Regression <- as.data.frame(do.call(rbind, Res_m.DDX3X_XIST))

# Label columns
colnames(f.Regression) <- c("pval_DDX3X", "R2_DDX3X")
colnames(m.Regression) <- c("pval_DDX3X","R2_DDX3X")

# ______________________________________________________________________________________________________________________
#  Plot DDX3X ~ XIST
# ______________________________________________________________________________________________________________________
# Find the max x/y values to set as x/y lim
Max_Func <- function(x, y){
  lim <- max(x[[y]])
  return(lim)
}

# x limits
f.xmax <- round(max(unlist(Map(Max_Func, x=f.DDX3X_XIST, y='XIST'))))
m.xmax <- round(max(unlist(Map(Max_Func, x=m.DDX3X_XIST, y='XIST'))))
# y limits
f.ymax <- round(max(unlist(Map(Max_Func, x=f.DDX3X_XIST, y='DDX3X'))))
m.ymax <- round(max(unlist(Map(Max_Func, x=m.DDX3X_XIST, y='DDX3X'))))

# Scatter plots: Correlation by tissue
Scatter_Func <- function(LM, TITLE, XMAX, YMAX, RESULTS){
  plot(LM$model$XIST, LM$model$DDX3X, main=TITLE, xlab='XIST', ylab='DDX3X',
       xlim=c(0, XMAX), ylim=c(0, YMAX))
  legend("bottomright", 
         bty="n", 
         legend=paste("R^2: ", 
                      format(RESULTS$r_2, digits=3), 
                      "; p_Val: ", 
                      format(RESULTS$p_val, digits=3)))
  abline(LM)
}

# Print plots
tiff(SCATTER.f, width = 500, heigh = 500)
f.Scatter <- Map(Scatter_Func, 
                 LM=lm_f.DDX3X_XIST, 
                 TITLE=names(lm_f.DDX3X_XIST), 
                 XMAX=f.xmax, 
                 YMAX=f.ymax, 
                 RESULTS=Res_f.DDX3X_XIST)
dev.off()

tiff(SCATTER.m, width = 500, heigh = 500)
m.Scatter <- Map(Scatter_Func, 
                 LM=lm_m.DDX3X_XIST, 
                 TITLE=names(lm_m.DDX3X_XIST), 
                 XMAX=m.xmax, 
                 YMAX=m.ymax, 
                 RESULTS=Res_m.DDX3X_XIST)
dev.off()

# ______________________________________________________________________________________________________________________
#  Summary table
# ______________________________________________________________________________________________________________________
# Add col with mean(DDX3X) and sd(DDX3X)
f.Regression$Mean_DDX3X <- lapply(f.DDX3X_Tissue_Counts, mean)
m.Regression$Mean_DDX3X <- lapply(m.DDX3X_Tissue_Counts, mean)

f.Regression$sd_DDX3X <- lapply(f.DDX3X_Tissue_Counts, sd)
m.Regression$sd_DDX3X <- lapply(m.DDX3X_Tissue_Counts, sd)

# Add column with number of tissues
Count_Rows <- function(x){
  res <- nrow(x)
  return(res)
}

Num_Fem <- as.character(lapply(f.Tissue_Lst, Count_Rows))
Num_Male <- as.character(lapply(m.Tissue_Lst, Count_Rows))

# Add column with number of samples per tissue
f.Regression <- cbind(Num_Tissues=Num_Fem, f.Regression) 
m.Regression <- cbind(Num_Tissues=Num_Male, m.Regression)

# Convert Num_Tissues from factor to numeric
f.Regression$Num_Tissues <- as.numeric(as.character(f.Regression$Num_Tissues))
m.Regression$Num_Tissues <- as.numeric(as.character(m.Regression$Num_Tissues))

# Convert rest of cols from list to numeric
f.Regression[,4:ncol(f.Regression)] <- lapply(f.Regression[,4:ncol(f.Regression)], function(x) unlist(x))
m.Regression[,4:ncol(m.Regression)] <- lapply(m.Regression[,4:ncol(m.Regression)], function(x) unlist(x))

# Add column indicating sample tissue type as first column
f.Regression <- cbind(Tissue=names(f.Tissue_Lst), f.Regression) 
m.Regression <- cbind(Tissue=names(m.Tissue_Lst), m.Regression) 

# Sort by number of tissues
f.Regression <- f.Regression[order(f.Regression$Num_Tissues, decreasing=TRUE),]
m.Regression <- m.Regression[order(m.Regression$Num_Tissues, decreasing=TRUE),]

# Write to file
write.csv(f.Regression, CORRELATION.f, row.names=FALSE)
write.csv(m.Regression, CORRELATION.m, row.names=FALSE)

# ______________________________________________________________________________________________________________________
#  Correlations summary
# ______________________________________________________________________________________________________________________
# Average R^2 of silenced genes reported in both studies for females and males
Summary.df <- data.frame(Female=colMeans(f.Regression[,2:ncol(f.Regression)]),
                         Male=colMeans(m.Regression[,2:ncol(m.Regression)]))
write.csv(Summary.df, SUMMARY)

# ______________________________________________________________________________________________________________________
#  Table of Slopes
# ______________________________________________________________________________________________________________________
# Get list of female and male tissue type samples
f.Tissues <- rownames(f.Regression)
m.Tissues <- rownames(m.Regression)

# Extract list of slopes
f.Slopes <- lapply(lm_f.DDX3X_XIST, function(x){
  res <- x[[1]][[2]]
  return(res)
})

m.Slopes <- lapply(lm_m.DDX3X_XIST, function(x){
  res <- x[[1]][[2]]
  return(res)
})

# Make and write df
l <- list(f.Slopes, m.Slopes)
Slopes.df <- rbindlist(l, use.names=TRUE, fill=TRUE, idcol="Sex")
Slopes.df$Sex <- c("Female", "Male")

write.csv(Slopes.df, SLOPES)

