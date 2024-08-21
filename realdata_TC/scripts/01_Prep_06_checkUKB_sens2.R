#' ---
#' title: "Get nice subset of TC data"
#' subtitle: "Longitudinal MVMR"
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
#' Here, I want to get a slightly different SNP set for my MVMR approach. When playing around with my SNP selection process, I found one setting in which I only had a significant effect of the variability on CAD. I want to recreate that selection. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load UKB data ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_03_SNPData.RData"))
SNPList = fread("../temp/01_Prep_06_SNPList_TC_test.txt",header = F)
table(is.element(SNPList$V1,pvar$ID))
SNPList[!is.element(V1,pvar$ID)]

pvar2 = copy(pvar)
pvar2 = pvar2[ID %in% SNPList$V1,]
pvar2[,comment_LD := "LD OK"]
filt = is.element(pvar$ID,pvar2$ID)
geno_mat2 = geno_mat[,filt]
dim(geno_mat2)

test = pvar[,.N,CHR]
test = test[N>1,]
myCHRs = test$CHR

LDTab2 = foreach(i = 1:length(myCHRs))%do%{
  #i=20
  myCHR = myCHRs[i]
  filt = pvar2$CHR == myCHR
  dumMatrix = geno_mat2[,filt]
  if(sum(filt)>1){
    CorTab = cor(dumMatrix)^2
    heatmap(CorTab,Rowv = NA,Colv = NA, main =paste0("Chromosome ",myCHRs[i]))
    CorTab2 = as.data.table(CorTab)
    CorTab2[,SNP1 := rownames(CorTab)]
    CorTab_long = melt(CorTab2,id.vars="SNP1",measure.vars=rownames(CorTab),variable.name="SNP2")
    CorTab_long = CorTab_long[SNP1 != SNP2,]
    CorTab2_long = CorTab_long[value>0.05,]
    pvar2[ID %in% CorTab2_long$SNP1 | ID %in% CorTab2_long$SNP2, comment_LD := "LD >0.05 with at least 1 SNP"]
    CorTab_long
  }
  
}
LDTab2 = rbindlist(LDTab2)
dev.off()

table(pvar2$comment_LD)

save(pvar2, psam, geno_mat2, file = paste0(UKB_phenotypes_filtered,"/01_Prep_06_SNPData_filtered_sens.RData"))
save(LDTab2, file = paste0("../results/01_Prep_06_LD_filtered.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

