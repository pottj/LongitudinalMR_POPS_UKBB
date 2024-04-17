#' ---
#' title: "Get Supplemental Tables"
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
#' **Supplemental Tables**
#' 
#' 1) Simulation results - main
#' 2) Simulation results - sensitivity
#' 3) POPS - Cohort description (mean and SD of all covariables used in the linear regression model - restricted to setting all2 and phenotype EFW_C)
#' 4) POPS - Overview of sample size per exposure and sample filter setting
#' 5) POPS - MVMR results - main model
#' 6) POPS - MVMR results - sens 1: GBR3
#' 7) POPS - MVMR results - sens 2: LMM (no var)
#' 8) POPS - MVMR results - sens 3: noTimeIA (no slope)
#' 9) POPS - MVMR results - sens 4: sigma function time independent
#' 10) POPS - MVMR results - sens 5: no covariables 
#' 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = paste0("S",c(1:9)),
                    Title = c("Simulation results - main",
                              "Simulation results - senstivity",
                              "POPS - Study description",
                              "POPS - MVMR main",
                              "POPS - MVMR sens1: GBR3 (only British ancestry and 3 measurements per sample available)",
                              "POPS - MVMR sens2: LMM (no var effect estimated)",
                              "POPS - MVMR sens3: noTimeIA (no slope effect estimated)",
                              "POPS - MVMR sens4: SigmaTimeIndep (no adjustment for time in sigma function)",
                              "POPS - MVMR sens5: noCovar (gamlss without any further adjustment, just SNP and time)"),
                    Source = c("simulation/results/_tables/",
                               "simulation_v2/results/_tables/",
                               "done in this script",
                               "realdata/results/04_MVMR_01_MAIN_240414.RData",
                               "realdata/results/04_MVMR_02_SENS_GBR3_240414.RData",
                               "realdata/results/04_MVMR_03_SENS_LMM_240414.RData",
                               "realdata/results/04_MVMR_05_SENS_noTimeIA_240414.RData",
                               "realdata/results/04_MVMR_04_SENS_SigmaTimeIndep_240414.RData",
                               "realdata/results/04_MVMR_07_SENS_noCovars_240414.RData"))
  
  tab0
  
  abbreviations = data.table(abbrev_paper = c("BPD_R","HC_R","AC_R","FL_R","EFW_L",
                                           "EFW_R","EFW_Z","EFW_C",
                                           "AC_Z","FL_Z","HCAC_Z","ACFL_Z",
                                           "BW_R","BW_Z","BW_C","eCS",
                                           "GA","mAge","mHeight","mSmoke","fSex"),
                             phenotype_POPS = c("bpd_cm","hc_cm","ac_cm","fl_cm","logefwcomb",
                                                "efwcomb","efwcombZv2","efwcombv2_cent",
                                                "acZv2","flZv2","hcacZv2","acflZv2",
                                                "pn_bw","BW_SDS_Br1990","BW_Centile_Br1990","pn_emcsall",
                                                "pn_ga_wk","an_est_age","an_height",
                                                "an_smokstat","pn_sex"),
                             comment = c("biparietal diameter (in cm)",
                                         "head circumference (in cm)",
                                         "abdominal circumference (in cm)",
                                         "femur length (in cm)",
                                         "estimated fetal weight (log-transformed)",
                                         "estimated fetal weight (in kg)",
                                         "estimated fetal weight (GA-adjusted Z-scores)",
                                         "estimated fetal weight (GA-adjusted percentiles)",
                                         "abdominal circumference (GA-adjusted Z-scores)",
                                         "femur length (GA-adjusted Z-scores)",
                                         "ratio of head to abdominal circumference (GA-adjusted Z-scores)",
                                         "ration of abdominal circumference to femur length (GA-adjusted Z-scores)", 
                                         "birth weight (in kg)",
                                         "birth weight (GA-adjusted Z-score)",
                                         "birth weight (GA-adjusted percentiles)",
                                         "emergency Cesarean Section (0=no; 1=yes)",
                                         "gestational age (in weeks)",
                                         "maternal age (in years)",
                                         "maternal height (in cm)",
                                         "maternal smoking status",
                                         "fetal sex (1=male; 2=female)"))
  
}

#' # Get Sup Tab 1 ####
#' ***
#' Simulation results - old
#' 
{
  loaded = load("simulation/results/_tables/Simulation_complete_F0.RData")
  tab1 = get(loaded)
  names(tab1) = gsub("X1","mean",names(tab1))
  names(tab1) = gsub("X2","slope",names(tab1))
  names(tab1) = gsub("X3","var",names(tab1))
  tab1
}

#' # Get Sup Tab 2 ####
#' ***
#' Simulation results - old
#' 
{
  loaded = load("simulation_v2/results/_tables/Simulation_complete_F0.RData")
  tab2 = get(loaded)
  names(tab2) = gsub("X1","mean",names(tab2))
  names(tab2) = gsub("X2","slope",names(tab2))
  names(tab2) = gsub("X3","var",names(tab2))
  tab2
}

#' # Get Sup Tab 3 ####
#' ***
#' Study description
#' 
#' a) for the exposures (per X complete): mean and SD per time point for the abs values only (Z scores should all be 0 and 1, is that meaningful? yes, because subgroup? no, because not too many samples filtered?)
#' b) for the covariables (for all2 and EFW complete): 
#'    - gestational age at birth
#'    - birth weight (abs, Z score, centile)
#'    - mother age 
#'    - mother height
#'    - mother smoking status 
#'    - gestational diabetes diet (yes/no)
#'    - babys sex (male/female)
#'    - emergency cesarian section (yes/no)
#' 
{
  newpath = gsub("../../","",POPS_phenotypes)
  myGD_files = list.files(path = newpath,pattern = "01_Prep_04")
  loaded1 = load(paste0(newpath,myGD_files[1]))
  loaded1
  loaded2 = load(paste0(newpath,myGD_files[2]))
  loaded2
  
  # get table for exposure
  myTab_X[,logefwcomb := log(efwcomb)]
  myTab_X[,efwcomb := efwcomb/1000]
  myTab_X[,efwcombv2_cent := efwcombv2_cent/100]
  
  myExposures=names(myTab_X)[c(18,38,23,24, 12:15, 25,26,29,30)]
  
  stopifnot(myExposures %in% abbreviations$phenotype_POPS)
  
  dumTab1 = foreach(i = 1:length(myExposures))%do%{
    #i=1
    meanTab=myTab_X[,round(mean(get(myExposures[i]),na.rm=T),3),by=scan]
    sdTab=myTab_X[,round(sd(get(myExposures[i]),na.rm=T),3),by=scan]
    nTab=myTab_X[!is.na(get(myExposures[i])),.N,by=scan]
    nTab2 =myTab_X[!is.na(get(myExposures[i])),.N,by=POPSID]
    nTab3 =myTab_X[!is.na(get(myExposures[i])) & ancestry=="GBR",.N,by=POPSID]
    n1 = dim(nTab2[N>=2])[1]
    n2 = dim(nTab2[N>2])[1]
    n3 = dim(nTab3[N>=2])[1]
    n4 = dim(nTab3[N>2])[1]
    filt = abbreviations$phenotype_POPS == myExposures[i]
    res = data.table(phenotype = abbreviations[filt,abbrev_paper],
                     variable = abbreviations[filt,phenotype_POPS],
                     description = abbreviations[filt,comment],
                     GA20_n = nTab[1,2],
                     GA20_mean = meanTab[1,2],
                     GA20_SD = sdTab[1,2],
                     GA28_n = nTab[2,2],
                     GA28_mean = meanTab[2,2],
                     GA28_SD = sdTab[2,2],
                     GA36_n = nTab[3,2],
                     GA36_mean = meanTab[3,2],
                     GA36_SD = sdTab[3,2],
                     sampleSize_all2 = n1,
                     sampleSize_all3 = n2,
                     sampleSize_GBR2 = n3,
                     sampleSize_GBR3 = n4)
    names(res)=gsub(".V1","",names(res))
    names(res)=gsub(".N","",names(res))
    res
  }
  tab3a = rbindlist(dumTab1)
  tab3a
  
  # get table for covariables and outcome (filtered for 2 or more time points for efwcomb 
  myTab_Y[,pn_bw := pn_bw/1000]
  myTab_Y[,BW_Centile_Br1990 := BW_Centile_Br1990/100]
  myTab_Y = myTab_Y[efwcomb>=2,]
  
  myVars1 = names(myTab_Y)[c(10,2,3,12,15,16)]
  dumTab2 = foreach(i = 1:length(myVars1))%do%{
    #i=1
    meanTab=myTab_Y[,round(mean(get(myVars1[i]),na.rm=T),3)]
    sdTab=myTab_Y[,round(sd(get(myVars1[i]),na.rm=T),3)]
    nTab1 = myTab_Y[!is.na(get(myVars1[i])),.N]
    nTab2 = myTab_Y[!is.na(get(myVars1[i])) & ancestry=="GBR",.N]
    
    filt = abbreviations$phenotype_POPS == myVars1[i]
    res = data.table(phenotype = abbreviations[filt,abbrev_paper],
                     variable = abbreviations[filt,phenotype_POPS],
                     description = abbreviations[filt,comment],
                     mean = meanTab,
                     SD = sdTab,
                     sampleSize_all = nTab1,
                     sampleSize_GBR = nTab2)
    res
  }
  tab3b = rbindlist(dumTab2)
  tab3b
  
  myVars2 = names(myTab_Y)[c(8,9,18)]
  dumTab3 = foreach(i = 1:length(myVars2))%do%{
    #i=1
    test=myTab_Y[,table(get(myVars2[i]))]
    test2 = test/dim(myTab_Y)[1]
    test2 = round(test2,3)
    nTab1 = myTab_Y[!is.na(get(myVars2[i])),.N]
    nTab2 = myTab_Y[!is.na(get(myVars2[i])) & ancestry=="GBR",.N]
    
    filt = abbreviations$phenotype_POPS == myVars2[i]
    res = data.table(phenotype = abbreviations[filt,abbrev_paper],
                     variable = abbreviations[filt,phenotype_POPS],
                     description = abbreviations[filt,comment],
                     number = test,
                     percent = test2,
                     sampleSize_all = nTab1,
                     sampleSize_GBR = nTab2)
    res[,percent.V1:=NULL]
    setnames(res,"number.V1","group")
    names(res)=gsub(".N","",names(res))
    res
  }
  tab3c = rbindlist(dumTab3)
  tab3c[group==1,group:="cases"]
  tab3c = tab3c[c(2,4,3,1,5,8),]
  tab3c
}

#' # Get Sup Tab 4 ####
#' ***
#' POPS MVMR tables
#' 
{
  myMVMR_results = list.files(path = "realdata/results/",pattern = "04_MVMR_0")
  myMVMR_results = myMVMR_results[grepl("240414",myMVMR_results)]
  
  dumTab4 = foreach(i = 1:length(myMVMR_results))%do%{
    #i=3
    load(paste0("realdata/results/",myMVMR_results[i]))
    names(MVMR_results)
    MVMR_results = MVMR_results[setting == "multivariate",]
    
    MVMR_tab_wide = dcast(MVMR_results, exposure + outcome + ID + threshold + NR_SNPs_total + HeteroStat + HeteroStat_pval ~ exposure_type, 
                          value.var=c("NR_SNPs_type","beta_IVW","SE_IVW","pval_IVW","condFstat"))
    names(MVMR_tab_wide) = gsub("IVW_","",names(MVMR_tab_wide))
    names(MVMR_tab_wide) = gsub("NR_SNPs_type","SNPs",names(MVMR_tab_wide))
    matched = match(MVMR_tab_wide$exposure,abbreviations$phenotype_POPS)
    MVMR_tab_wide[,exposure := abbreviations[matched,abbrev_paper]]
    
    matched = match(MVMR_tab_wide$outcome,abbreviations$phenotype_POPS)
    MVMR_tab_wide[,outcome := abbreviations[matched,abbrev_paper]]
    
    tab4 = copy(MVMR_tab_wide)
    
  }
  tab4 = rbindlist(dumTab4,fill=T)
  names(tab4)
  tab4 = tab4[,c(1:7, 
                                   8,11,14,17,20,
                                   9,12,15,18,21, 
                                   10,13,16,19,22)]
  
  
}  

#' # Save tables ###
#' ***

tosave4 = data.table(data = c("tab0", "tab1", "tab2","tab3a","tab3b","tab3c", "tab4"), 
                     SheetNames = c("Content","TableS1", "TableS2", "TableS3a", "TableS3b",
                                    "TableS3c","TableS4"))
excel_fn = paste0("SupplementalTables_",tag,".xlsx")
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(tab0, tab1, tab2, tab3a, tab3b, tab3c, tab4, 
     file = paste0("SupplementalTables_",tag,".RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
