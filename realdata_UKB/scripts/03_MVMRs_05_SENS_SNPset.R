#' ---
#' title: "MVMR - sensitivity analysis: SNP selection"
#' subtitle: "Longitudinal MVMR in UKB"
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

source("../../SourceFile.R")
source("../../helperfunctions/MVMR_jp_POPS.R")

#' # Get data ####
#' ***
#' ## SNP info
#' Get high quality SNPs only
load("../results/01_Prep_04_LD.RData")
LDTab[,SNP2 := as.character(SNP2)]
load("../results/01_Prep_03_SNPList.RData")

#' ## Exposure
load("../results/02_SNPs_01_MAIN.RData")
length(unique(myAssocs_X$SNP))

#' SNP set 1: distinct SNP sets for mean and variability
mySNPlist1_distinct = myAssocs_X[(pval_mean<5e-8 & pval_var>0.05) | (pval_mean>0.05 & pval_var<5e-8),SNP]

#' SNP set 2: top20, best associated SNPs per type
myAssocs_X[,NR := 1:346]
setorder(myAssocs_X,pval_mean)
myAssocs_X[1:20,flag:=T]
setorder(myAssocs_X,pval_slope)
myAssocs_X[1:20,flag:=T]
setorder(myAssocs_X,pval_var)
myAssocs_X[1:20,flag:=T]
setorder(myAssocs_X,NR)
mySNPlist2_top20 = myAssocs_X[flag==T,SNP]

#' SNP set 3: FUMA annotation with other CAD risk factors in the GWAS Catalog trait (potential pleiotropy)
FUMA = fread("../temp/gwascatalog.txt")
filt1 = grepl("alcohol",FUMA$Trait)
filt2 = grepl("smoking",FUMA$Trait)
filt3 = grepl("pressure",FUMA$Trait)
filt4 = grepl("diabetes",FUMA$Trait)
filt5 = grepl("obes",FUMA$Trait) | grepl("BMI",FUMA$Trait) | grepl("body mass index",FUMA$Trait)
filt6 = grepl("physical activity",FUMA$Trait)
FUMA = FUMA[filt1 | filt2 | filt3 | filt4 | filt5 | filt6,]

mySNPlist3_FUMA = myAssocs_X[SNP %in% FUMA[,IndSigSNP],SNP]

#' SNP set 4: known lipid biology (drug targets: PCSK9, HMGCR, NPC1L1, and CETP; other obvious hits: APOB, LPA, FADS2, LDLR, TM6SF2, APOE)
#' 
mySNPlist4_bio = c("rs11591147","rs12916","rs17725246","rs247616",
                   "rs934197","rs10455872","rs174564","rs964184","rs73015024","rs58542926","rs7412")
mySNPlist4_bioGenes = c("PCSK9","HMGCR","NPC1L1","CETP",
                        "APOB","LPA","FADS2","APOA5","LDLR","TM6SF2","APOE")

#' SNP set 4: sex interaction using Kanoni et al. data
TC_male = fread(GLGC_TC_males)
TC_male = TC_male[rsID %in% SNPList$rsID]
TC_female = fread(GLGC_TC_females)
TC_female = TC_female[rsID %in% SNPList$rsID]

TC_male[duplicated(rsID),]
TC_male[rsID == "rs9832727"]
TC_female[rsID == "rs9832727"]
SNPList[rsID == "rs9832727",]

TC_male = TC_male[!duplicated(rsID),]
TC_female = TC_female[!duplicated(rsID),]

stopifnot(TC_male$rsID == TC_female$rsID)
table(TC_male$REF == TC_female$REF,
      TC_male$REF == TC_female$ALT)

stopifnot(TC_male$rsID == SNPList$rsID)
table(TC_male$ALT == SNPList$effect_allele, TC_male$REF == SNPList$effect_allele)
table(TC_male$ALT == SNPList$other_allele, TC_male$REF == SNPList$other_allele)

filt = TC_male$ALT != SNPList$effect_allele
TC_male[,EAF := POOLED_ALT_AF]
TC_male[filt,EAF := 1-EAF]
TC_male[,beta := EFFECT_SIZE]
TC_male[filt,beta := beta*(-1)]
plot(SNPList$UKB_TC_EAF, TC_male$EAF)

stopifnot(TC_female$rsID == SNPList$rsID)
table(TC_female$ALT == SNPList$effect_allele, TC_female$REF == SNPList$effect_allele)
table(TC_female$ALT == SNPList$other_allele, TC_female$REF == SNPList$other_allele)

filt = TC_female$ALT != SNPList$effect_allele
TC_female[,EAF := POOLED_ALT_AF]
TC_female[filt,EAF := 1-EAF]
TC_female[,beta := EFFECT_SIZE]
TC_female[filt,beta := beta*(-1)]
plot(SNPList$UKB_TC_EAF, TC_female$EAF)

plot(TC_male$EAF,TC_female$EAF)
min(TC_female$EAF)
max(TC_female$EAF)
min(TC_male$EAF)
max(TC_male$EAF)

plot(TC_male$beta,TC_female$beta)
cor.test(TC_male$beta,TC_female$beta)

# build information
TC_sexIA = copy(SNPList)[,c(1:6)]

TC_sexIA[,trait1 := "TC_males"]
TC_sexIA[,trait1_beta := TC_male[,beta]]
TC_sexIA[,trait1_SE := TC_male[,SE]]
TC_sexIA[,trait1_pval := TC_male[,pvalue]]
TC_sexIA[,trait1_n := TC_male[,N]]

TC_sexIA[,trait2 := "TC_females"]
TC_sexIA[,trait2_beta := TC_female[,beta]]
TC_sexIA[,trait2_SE := TC_female[,SE]]
TC_sexIA[,trait2_pval := TC_female[,pvalue]]
TC_sexIA[,trait2_n := TC_female[,N]]

TC_sexIA[,IA_diff := trait2_beta - trait1_beta]
TC_sexIA[,IA_SE := sqrt(trait1_SE^2 + trait2_SE^2) ]
TC_sexIA[,IA_Zscore := IA_diff/IA_SE ]
TC_sexIA[,IA_pval := pnorm(abs(IA_Zscore), lower.tail = F)*2 ]

mySNPlist5_sexIA = TC_sexIA[IA_pval<1e-6,rsID]

load(paste0("../temp/03_MVMRInput_MAIN.RData"))

#' # Do MVMR ####
#' ***
mySampleSize = c(68467)
myOutcomes = unique(myAssocs_Y_long$type)
myFlag = paste0("sens_SNPset_",c("distinct","top20","pleiotropy","bioknown","sexIA"))
mySet = c("mySNPlist1_distinct","mySNPlist2_top20","mySNPlist3_FUMA","mySNPlist4_bio","mySNPlist5_sexIA")

names(myAssocs_Y_long) = c("SNP", "phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean" )

# filter data
myAssocs_X_long[,phenotype := paste0("TC_",model)]

dumTab2 = foreach(j = 1:length(myFlag))%do%{
  #j=1
  myAssocs_X_long2 = copy(myAssocs_X_long)
  myAssocs_X_long2[,dumID := myFlag[j]]
  myAssocs_X_long2 = myAssocs_X_long2[SNP %in% get(mySet[j])]
  
  dumTab3 = foreach(k = 1:length(myOutcomes))%do%{
    #k=1
    myOutcome = myOutcomes[k]
    myExposure = unique(myAssocs_X_long2$phenotype)
    myAssocs_Y2 = copy(myAssocs_Y_long)
    myAssocs_Y2 = myAssocs_Y2[phenotype == myOutcome,]
    
    # do MVMR
    MVMR0 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = myFlag[j],
                         GX_pval_treshold = 1,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = mySampleSize,
                         random = F,getCondF = T,getUni = T)
    
    MVMR2 = MVMR_jp_POPS(data_exposure = myAssocs_X_long2,
                         data_outcome = myAssocs_Y2,
                         exposure_name = myExposure, 
                         outcome_name = myOutcome,
                         flag = myFlag[j],
                         GX_pval_treshold = 5e-8,
                         getPlot = F,
                         corTab = LDTab,
                         corTab_threshold = 0.1,sampleSize_GX = mySampleSize,
                         random = F,getCondF = T,getUni = T)
    
    MVMR0[,threshold := "all_SNPs"]
    MVMR2[,threshold := "gw_SNPs"]
    MVMR = rbind(MVMR0,MVMR2,fill=T)
    MVMR
    
  }
  MVMR_results2 = rbindlist(dumTab3,fill = T)
  MVMR_results2
}

MVMR_results = rbindlist(dumTab2,fill = T)
save(MVMR_results,file = paste0("../results/03_MVMR_05_SENS_SNPsets.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
