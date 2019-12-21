# To:Do: Need to add sample sex to plot titles for this to be useful
# Split list (lm.MeanX_XIST) by sex and plot seperately

# Constants
DATA <- "~/XIST/XIST_Indiv_121819.RData"
QQ <- "~/XIST/Indiv/Indiv_QQPlots.pdf"

# # Function to extract resduals from list of LM objects 
# Res_Func <- function(x){
#   tmp <- resid(x)
#   return(tmp)
# }
# Res.lst <- lapply(lm.MeanX_XIST, Res_Func)

# Functon to plot LM residuals
QQ_Func <- function(LM, INDIV){
  plot <- qqnorm(LM[['residuals']],
                 ylab = "Residuals",
                 xlab = "Normal Scores",
                 main = paste("lm(MeanX ~ XIST) in", INDIV))
  qqline(LM[['residuals']])
  return(plot)
}
#pdf(QQ)
Map(QQ_Func, LM=lm.MeanX_XIST, INDIV=names(lm.MeanX_XIST))
#dev.off()




