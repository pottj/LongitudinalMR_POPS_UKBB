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
source("../../helperfunctions/Venn_hk.R")

#' # Load data ####
#' ***
load("../results/_tables/SupplementalTables_EGG.RData")
myNames = names(tab2)[grepl("POPS_EFW_beta",names(tab2))]

data_wide1 = dcast(tab2[POPS_EFW_phenotype == "logefwcomb"], rsID ~ flag, value.var=myNames)
names(data_wide1) = gsub("POPS_EFW_beta_","",names(data_wide1))
data_wide_matrix1 = as.matrix(data_wide1[,-1])

data_wide2 = dcast(tab2[POPS_EFW_phenotype == "efwcombZv2"], rsID ~ flag, value.var=myNames)
names(data_wide2) = gsub("POPS_EFW_beta_","",names(data_wide2))
data_wide_matrix2 = as.matrix(data_wide2[,-1])

#' # Get correlation plot ####
#' ***
CorTab1 = cor(data_wide_matrix1,use = "pairwise.complete.obs")
filt1 = rownames(CorTab1) %in% c("slope_1A - SENS - no slope","var_1B - SENS - no variability") 
CorTab1 = CorTab1[!filt1,!filt1]
filt2 = grepl("ageCorrected",rownames(CorTab1))
CorTab1 = CorTab1[!filt2,!filt2]

colnames(CorTab1) = rep(" ",10)
corrplot(CorTab1,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/SupFigs/GenCorPlot_EWFL.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab1,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)
dev.off()

CorTab2 = cor(data_wide_matrix2,use = "pairwise.complete.obs")
filt1 = rownames(CorTab2) %in% c("slope_1A - SENS - no slope","var_1B - SENS - no variability") 
CorTab2 = CorTab2[!filt1,!filt1]
filt2 = grepl("ageCorrected",rownames(CorTab2))
CorTab2 = CorTab2[!filt2,!filt2]

colnames(CorTab2) = rep(" ",10)
corrplot(CorTab2,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)

png(file=paste0("../results/_figures/SupFigs/GenCorPlot_EWFZ.png"),
    width=1500,height=900,res = 200)
corrplot(CorTab2,order = "hclust",
         #col= colorRampPalette(c("#7D0722","white", "#053061"))(10),
         tl.col = "black", tl.srt = 45)
dev.off()

#' # Venn diagram ####
#' ***
mean = tab2[flag=="0 - MAIN" & POPS_EFW_pval_mean<0.05 & POPS_EFW_phenotype == "logefwcomb",rsID]
slope = tab2[flag=="0 - MAIN" & POPS_EFW_pval_slope<0.05 & POPS_EFW_phenotype == "logefwcomb",rsID]
var = tab2[flag=="0 - MAIN" & POPS_EFW_pval_var<0.05 & POPS_EFW_phenotype == "logefwcomb",rsID]

venn3(x1=mean,y1=slope,z1=var)

mean = tab2[flag=="0 - MAIN" & POPS_EFW_pval_mean<0.05 & POPS_EFW_phenotype == "efwcombZv2",rsID]
slope = tab2[flag=="0 - MAIN" & POPS_EFW_pval_slope<0.05 & POPS_EFW_phenotype == "efwcombZv2",rsID]
var = tab2[flag=="0 - MAIN" & POPS_EFW_pval_var<0.05 & POPS_EFW_phenotype == "efwcombZv2",rsID]

venn3(x1=mean,y1=slope,z1=var)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
