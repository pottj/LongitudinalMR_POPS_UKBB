#' ---
#' title: "Try 2-sample MVMR"
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
#' Here, I want to run the MVMR using the statistics from all samples with at least 2 exposure measurements on the EGG and UKB outcomes. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
source("../../helperfunctions/MVMR_jp_POPS.R")
source("../../helperfunctions/MVMR_jp_POPS_top20.R")

#' # Get data ####
#' ***
#' ## SNP info
#' Get high quality SNPs only
load("../results/01_Prep_02_LD_filtered_EGG.RData")
LDTab = copy(LDTab2)
load("../results/01_Prep_02_SNPList_filtered_EGG.RData")
SNPList = copy(SNPList_filtered)

#' Filter LD table for good SNPs only
LDTab[,SNP2 := as.character(SNP2)]

#' ## Exposure
load("../results/02_SNPs_01_MAIN.RData")
myAssocs_X = copy(myAssocs_X_gamlssIA)

#' Check the GX associations
myAssocs_X[,table(pval_mean==0,phenotype)]
myAssocs_X[,table(pval_mean<1e-50,phenotype)]

#' Transform into wide format
data_long1 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("beta_mean", "beta_slope", "beta_var"),
                  variable.name="type",
                  value.name="beta")
data_long2 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("SE_mean", "SE_slope", "SE_var"),
                  variable.name="type",
                  value.name="SE")
data_long3 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("tval_mean", "tval_slope", "tval_var"),
                  variable.name="type",
                  value.name="tval")
data_long4 = melt(myAssocs_X,
                  id.vars=names(myAssocs_X)[1:7],
                  measure.vars=c("pval_mean", "pval_slope", "pval_var"),
                  variable.name="type",
                  value.name="pval")

myAssocs_X_long = cbind(data_long1,data_long2[,9],
                        data_long3[,9],data_long4[,9])
myAssocs_X_long[,type := gsub("beta_","",type)]
myAssocs_X_long = myAssocs_X_long[!is.na(beta),]
# myAssocs_X_long[type == "slope",beta:=beta*36]
# myAssocs_X_long[type == "slope",SE:=SE*36]

#' ## Outcomes
#' 
#' 1) Sakaue et al. (2021): Cesarian section (meta-analyses)
#' 2) Neale lab: Delivery methods: Elective caesarean section
#' 3) Neale lab: Delivery methods: Emergency caesarean section
#' 
SNPList[,chrPosEAOA_b37:= paste(chr,pos_b37,effect_allele,other_allele,sep=":")]
SNPList[,chrPosOAEA_b37:= paste(chr,pos_b37,other_allele,effect_allele,sep=":")]

#' ### Sakaue
Sakaue = fread(paste0(pathData,"34594039-GCST90018810-EFO_0009636.h.tsv.gz"),
               header=T, select = c(1:7,11,19,20))

Sakaue = Sakaue[hm_rsid %in% SNPList$rsid,]
matched = match(SNPList$rsid,Sakaue$hm_rsid)
Sakaue = Sakaue[matched,]
stopifnot(Sakaue$hm_rsid == SNPList$rsid)
stopifnot(Sakaue$hm_pos == SNPList$pos_b38)
plot(Sakaue$hm_pos,SNPList$pos_b38)
stopifnot(Sakaue$hm_effect_allele == SNPList$effect_allele)
plot(Sakaue$hm_effect_allele_frequency,SNPList$EAF2)

myAssocs_Y1 = copy(Sakaue)
myAssocs_Y1[,phenotype := "CS_Sakaue"]
myAssocs_Y1[,SNP:= SNPList$ID]
myAssocs_Y1[,rsID:= SNPList$rsid]
myAssocs_Y1[,sampleSize:= 243071]
myAssocs_Y1[,tval:= hm_beta/standard_error]
myAssocs_Y1 = myAssocs_Y1[,c(12,13,11,14,7,9,15,10)]
names(myAssocs_Y1) = c("SNP","rsID","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' ### UKB Elective caesarean section
#' SNPList: dumID2 = chr : pos_b38 : OA : EA
#' UKB: variant = chr : pos_b37 : OA : EA
UKB_elec = fread(paste0(pathData,"41221_7.gwas.imputed_v3.female.tsv.bgz"))
UKB_elec = UKB_elec[variant %in% SNPList$chrPosOAEA_b37,]
stopifnot(UKB_elec$variant == SNPList$chrPosOAEA_b37)

UKB_elec[,EA := gsub(".*:","",variant)]
stopifnot(UKB_elec$EA == SNPList$effect_allele)
UKB_elec[EA==minor_allele, EAF := minor_AF]
UKB_elec[EA!=minor_allele, EAF := 1-minor_AF]
plot(UKB_elec$EAF,SNPList$EAF2)

myAssocs_Y2 = copy(UKB_elec)
myAssocs_Y2[,phenotype := "CS_UKB_elective"]
myAssocs_Y2[,SNP:= SNPList$ID]
myAssocs_Y2[,rsID:= SNPList$rsid]
myAssocs_Y2 = myAssocs_Y2[,c(16,17,15,6,9:12)]
names(myAssocs_Y2) = c("SNP","rsID","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' ### UKB Emergency caesarean section
UKB_emer = fread(paste0(pathData,"41221_8.gwas.imputed_v3.female.tsv.bgz"))
UKB_emer = UKB_emer[variant %in% SNPList$chrPosOAEA_b37,]
stopifnot(UKB_emer$variant == SNPList$chrPosOAEA_b37)

UKB_emer[,EA := gsub(".*:","",variant)]
stopifnot(UKB_emer$EA == SNPList$effect_allele)
UKB_emer[EA==minor_allele, EAF := minor_AF]
UKB_emer[EA!=minor_allele, EAF := 1-minor_AF]
plot(UKB_emer$EAF,SNPList$EAF2)

myAssocs_Y3 = copy(UKB_emer)
myAssocs_Y3[,phenotype := "CS_UKB_emergency"]
myAssocs_Y3[,SNP:= SNPList$ID]
myAssocs_Y3[,rsID:= SNPList$rsid]
myAssocs_Y3 = myAssocs_Y3[,c(16,17,15,6,9:12)]
names(myAssocs_Y3) = c("SNP","rsID","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' ### Combine data
#myAssocs_Y = rbind(myAssocs_Y2,myAssocs_Y3)
myAssocs_Y = rbind(myAssocs_Y1,myAssocs_Y2,myAssocs_Y3)

#' ## save as temporary files
save(myAssocs_X_long,myAssocs_Y, file = paste0("../temp/07_MVMRInput_2SampleMR_CS.RData"))

#' # Do MVMR ####
#' ***
#' create a to do list
#' 
myExposures = unique(myAssocs_X_long$phenotype)
myOutcomes = unique(myAssocs_Y$phenotype)
myFlag = "main_2SMR_CS"

#registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
#registerDoParallel(4)

#dumTab2 = foreach(j = 1:length(myExposures))%dorng%{
dumTab2 = foreach(j = 1:length(myExposures))%do%{
  #j=1
  source("../../SourceFile_HPC.R")
  source("../../helperfunctions/MVMR_jp_POPS.R")
  source("../../helperfunctions/MVMR_jp_POPS_top20.R")
  
  myExposure = myExposures[j]
  
  # filter data
  myAssocs_X_long2 = copy(myAssocs_X_long)
  myAssocs_X_long2 = myAssocs_X_long2[phenotype == myExposure,]
  myAssocs_X_long2[,dumID := myFlag]
  
  dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
    #k=1
    myOutcome = myOutcomes[k]
    myAssocs_Y2 = copy(myAssocs_Y)
    myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
    
    message("Working on exposure ",myExposure," and outcome ",myOutcome," ...")
    
    # do MVMRs
    MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = myFlag,
                         GX_pval_treshold = 1,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 2996,random = F,getCondF = T,getUni = T)
    
    MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = myFlag,
                         GX_pval_treshold = 0.05,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = 2996,random = F,getCondF = T,getUni = T)
    
    MVMR4 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2,
                               data_outcome = myAssocs_Y2,
                               exposure_name = myExposure, 
                               outcome_name = myOutcome,
                               flag = myFlag,
                               SNPSets = "overlap",
                               getPlot = F,
                               corTab = LDTab,
                               corTab_threshold = 0.1,sampleSize_GX = 2996,random = F,getCondF = T,getUni = T)
    
    MVMR5 = MVMR_jp_POPS_top20(data_exposure = myAssocs_X_long2,
                               data_outcome = myAssocs_Y2,
                               exposure_name = myExposure, 
                               outcome_name = myOutcome,
                               flag = myFlag,
                               SNPSets = "distinct",
                               getPlot = F,
                               corTab = LDTab,
                               corTab_threshold = 0.1,sampleSize_GX = 2996,random = F,getCondF = T,getUni = T)
    
    MVMR0[,threshold := "all_SNPs"]
    MVMR2[,threshold := "nominal_SNPs"]
    MVMR4[,threshold := "top20_overlap"]
    MVMR5[,threshold := "top20_distinct"]
    MVMR = rbind(MVMR0,MVMR2,MVMR4,MVMR5,fill=T)
    MVMR
    
  }
  MVMR_Tab1 = rbindlist(dumTab3)
  MVMR_Tab1
  
}

MVMR_results = rbindlist(dumTab2,fill = T)
save(MVMR_results,file = paste0("../results/07_MVMR_2SampleMR_CS.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
