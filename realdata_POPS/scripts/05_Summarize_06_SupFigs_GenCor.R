#' ---
#' title: "Get Correlation Plots"
#' subtitle: "Longitudinal MVMR in POPS"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' **Correlation Plot between raw - log - Z-scores and mean - slope - variability**
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
load("../results/_tables/SupplementalTables_EGG.RData")
myNames = names(tab2)[grepl("POPS_EFW_beta",names(tab2))]
data_wide = dcast(tab2, rsID ~ flag, value.var=myNames)
data_wide
names(data_wide) = gsub("POPS_EFW_beta_","",names(data_wide))
data_wide_matrix = as.matrix(data_wide[,-1])

#' # Get correlation plot ####
#' ***
CorTab = cor(data_wide_matrix,use = "pairwise.complete.obs")
filt1 = rownames(CorTab) %in% c("slope_SENS_noSlope","var_SENS_noVar") 
CorTab = CorTab[!filt1,!filt1]
colnames(CorTab) = rep(" ",13)
corrplot(CorTab,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/SupFigs/GenCorPlot_EWF.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)
dev.off()

#' # Venn diagram ####
#' ***
source("../../helperfunctions/Venn_hk.R")
mean = tab2[flag=="MAIN" & POPS_EFW_pval_mean<0.05,rsID]
slope = tab2[flag=="MAIN" & POPS_EFW_pval_slope<0.05,rsID]
var = tab2[flag=="MAIN" & POPS_EFW_pval_var<0.05,rsID]

venn3(x1=mean,y1=slope,z1=var)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
