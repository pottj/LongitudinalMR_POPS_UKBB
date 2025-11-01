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

source("../../SourceFile.R")

#' # Load data ####
#' ***
load("../results/_tables/SupplementalTables_GLGC.RData")
tab2 = tab2[flag %in% c("0 - MAIN","2 - SENS - no statins","1A - SENS - no slope","1B - SENS - no variability")]
myNames = names(tab2)[grepl("UKB_TC_beta",names(tab2))]

data_wide = dcast(tab2, rsID ~ flag, value.var=myNames)
data_wide
names(data_wide) = gsub("UKB_TC_beta_","",names(data_wide))
data_wide_matrix = as.matrix(data_wide[,-1])

#' # Get correlation plot ####
#' ***
CorTab = cor(data_wide_matrix,use = "pairwise.complete.obs")
filt1 = rownames(CorTab) %in% c("slope_1A - SENS - no slope","var_1B - SENS - no variability") 
CorTab = CorTab[!filt1,!filt1]
filt2 = grepl("ageCorr",rownames(CorTab))
CorTab = CorTab[!filt2,!filt2]

colnames(CorTab) = rep(" ",10)
corrplot(CorTab,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/SupFigs/GenCorPlot_TC.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)
dev.off()

#' # Venn diagram ####
#' ***
source("../../helperfunctions/Venn_hk.R")
mean = tab2[flag=="0 - MAIN" & UKB_TC_pval_mean<5e-8,rsID]
slope = tab2[flag=="0 - MAIN" & UKB_TC_pval_slope<5e-8,rsID]
var = tab2[flag=="0 - MAIN" & UKB_TC_pval_var<5e-8,rsID]

venn3(x1=mean,y1=slope,z1=var)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
