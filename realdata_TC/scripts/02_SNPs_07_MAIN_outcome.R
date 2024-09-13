#' ---
#' title: "Prep PLINK output for CAD"
#' subtitle: "Longitudinal MR"
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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load data ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_03_SNPData_filtered.RData"))
pvar_main = copy(pvar2)
load(paste0(UKB_phenotypes_filtered,"/01_Prep_06_SNPData_filtered_sens.RData"))
pvar_sens1 = copy(pvar2)

myList = data.table(gene = c("APOE","LDLR","CELSR2","APOC1","PCSK9",
                             "APOB","TM6SF2","CETP","HMGCR","DOCK7",
                             "TOMM40","APOC4","ABCG8","LIPC","ABO",
                             "BUD13","BACE1","APOC3","APOA4","LPL"), 
                    SNP = c("rs7412","rs73015024","rs12740374","rs141622900","rs11591147",
                            "rs6548010","rs58542926","rs183130","rs12916","rs995000",
                            "rs11668327","rs12721109","rs4245791","rs633695","rs2519093",
                            "rs964184","rs116987336","rs12718462","rs61905132","rs6993414"), 
                    varEffect = c(T,F,F,T,F,
                                  F,F,F,F,F,
                                  T,F,F,F,F,
                                  T,T,T,T,T),
                    cytoband = c("19q13.32","19p13.2" ,"01p13.3","19q13.32","01p32.3",
                                 "02p24.1" ,"19p13.11","16q13"  ,"05q13.3" ,"01p31.3",
                                 "19q13.32","19q13.32","02p21"  ,"15q21.3" ,"09q34.2",
                                 "11q23.3" ,"11q23.3", "11q23.3","11q23.3", "08p21.3"))
load(paste0(UKB_phenotypes_filtered,"/01_Prep_03_SNPData.RData"))
filt = is.element(pvar$ID,myList$SNP)
table(filt)
pvar_sens2 = pvar[filt,]

pvar = rbind(pvar_main,pvar_sens1,pvar_sens2)
pvar = pvar[!duplicated(ID),]
setorder(pvar,CHR,POS)

#' load PLINK outfiles
myFiles = list.files(path = UKB_genotypes_filtered, pattern = "01_Prep_09")
myFiles
myFiles = myFiles[-4]

dumTab = foreach(i = 1:3)%do%{
  #i=1
  erg1 = fread(paste0(UKB_genotypes_filtered,myFiles[i]))
  erg1 = erg1[ID %in% pvar$ID,]
  erg1[,BETA := log(OR)]
  setnames(erg1,"LOG(OR)_SE","SE")
  
  # check alleles
  stopifnot(erg1$ID == pvar$ID)
  stopifnot(erg1$REF == pvar$REF)
  stopifnot(erg1$ALT == pvar$ALT)
  table(erg1$A1 == pvar$ALT)
  
  plot(pvar$EAF,erg1$A1_FREQ)
  plot(pvar$MAF,erg1$A1_FREQ)
  table(erg1$A1 == pvar$ALT, pvar$EAF>0.5)
  filt = pvar$EAF>0.5
  erg1[filt,BETA := BETA * (-1)]
  erg1[filt,A1_FREQ:= 1-A1_FREQ]
  erg1[filt,Z_STAT:= Z_STAT * (-1)]
  
  # filter for relevant columns
  erg1[,model := gsub(".*_","",myFiles[i])]
  erg1[,model := gsub(".glm.logistic.hybrid","",model)]

  erg1 = erg1[,c(3,21,15,20,17,18,19)]
  names(erg1) = c("SNP","model","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean")
  erg1
}
myAssocs_Y = rbindlist(dumTab)

save(myAssocs_Y,file=paste0("../results/02_SNPs_07_MAIN_outcome.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
