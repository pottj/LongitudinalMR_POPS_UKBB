#' ---
#' title: "Get Supplemental Tables for real data"
#' subtitle: "Longitudinal MVMR in UKBB"
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
#' **Supplemental Tables of real data analyses**
#' 
#' 1) UKBB - Cohort description (mean and SD of all covariables used in the regression model - restricted to setting main setting)
#' 2) UKBB - Overview of instruments
#' 3) UKBB - MVMR results (main + sensitivity)
#' 4) UKBB - MR results (main + sensitivity)
#' 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Get content table (tab0) ####
#' ***
{
  tab0 = data.table(Table = paste0("S",c(1:4)),
                    Title = c("UKBB - Study description",
                              "UKBB - SNP summary statistics",
                              "UKBB - MVMR results",
                              "UKBB - MR results"),
                    Source = c("done in this script",
                               "realdata_GLGC/results/02_SNPs_*.RData",
                               "realdata_GLGC/results/03_MVMR_*.RData",
                               "realdata_GLGC/results/03_MVMR_*.RData"))
  
  tab0
  
}

#' # Get Sup Tab 1 ####
#' ***
#' Study description
#' 
#' I want for the exposure (MAIN and SENS_sampleSet): 
#' 
#' - total sample size
#' - median (IQR) of number of repeated measurements per individual 
#' - mean (SD) of TC measurements
#' - Female %
#' - mean (SD) age
#' - Statin medication % (any time point
#' - mean (SD) age of start of statin treatment
#' 
{
  load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_GP_TC_GLGC.RData"))
  
  test1 = myTab6[,.N,by=BSU_ID]
  test2 = copy(myTab6)
  test2 = test2[lipLowMed == 1,]
  test2 = test2[!duplicated(BSU_ID),]
  
  tab1a = data.table(setting = "MAIN", 
                     sampleSize = dim(myTab7)[1],
                     NR_measurements_median = as.numeric(summary(test1$N)[3]),
                     NR_measurements_1stQ = as.numeric(summary(test1$N)[2]),
                     NR_measurements_3rdQ = as.numeric(summary(test1$N)[5]),
                     TC_mean = mean(myTab6$exposure_value),
                     TC_SD = sd(myTab6$exposure_value),
                     Female_absolute = dim(myTab7[sex==0,])[1],
                     Female_percent = dim(myTab7[sex==0,])[1]/dim(myTab7)[1],
                     Age_mean = mean(myTab6$exposure_age),
                     Age_SD = sd(myTab6$exposure_age),
                     Statin_absolut = dim(test2)[1],
                     Statin_percent = dim(test2)[1]/dim(myTab7)[1],
                     Age_1st_statin_mean = mean(test2$exposure_age),
                     Age_1st_statin_SD = sd(test2$exposure_age) )
  
  load(paste0(UKB_phenotypes_filtered,"/01_Prep_05_UKB_GP_TC_GLGC_sens.RData"))
  
  test1 = myTab6[,.N,by=BSU_ID]
  test2 = copy(myTab6)
  test2 = test2[lipLowMed == 1,]
  test2 = test2[!duplicated(BSU_ID),]
  
  tab1b = data.table(setting = "SENS_SampleSet", 
                     sampleSize = dim(myTab7)[1],
                     NR_measurements_median = as.numeric(summary(test1$N)[3]),
                     NR_measurements_1stQ = as.numeric(summary(test1$N)[2]),
                     NR_measurements_3rdQ = as.numeric(summary(test1$N)[5]),
                     TC_mean = mean(myTab6$exposure_value),
                     TC_SD = sd(myTab6$exposure_value),
                     Female_absolute = dim(myTab7[sex==0,])[1],
                     Female_percent = dim(myTab7[sex==0,])[1]/dim(myTab7)[1],
                     Age_mean = mean(myTab6$exposure_age),
                     Age_SD = sd(myTab6$exposure_age),
                     Statin_absolut = dim(test2)[1],
                     Statin_percent = dim(test2)[1]/dim(myTab7)[1],
                     Age_1st_statin_mean = mean(test2$exposure_age),
                     Age_1st_statin_SD = sd(test2$exposure_age) )
  
  tab1 = rbind(tab1a,tab1b)
  
  # round the numbers
  tab1[,TC_mean := round(TC_mean,2)]
  tab1[,TC_SD := round(TC_SD,2)]
  tab1[,Female_percent := round(Female_percent,3)*100]
  tab1[,Age_mean := round(Age_mean,2)]
  tab1[,Age_SD := round(Age_SD,2)]
  tab1[,Statin_percent := round(Statin_percent,3)*100]
  tab1[,Age_1st_statin_mean := round(Age_1st_statin_mean,2)]
  tab1[,Age_1st_statin_SD := round(Age_1st_statin_SD,2)]
  
  # merge some columns
  tab1[,measurements := paste0(NR_measurements_median," (",NR_measurements_1stQ,", ",NR_measurements_3rdQ,")")]
  tab1[,TC := paste0(TC_mean," (",TC_SD,")")]
  tab1[,Female := paste0(Female_absolute," (",Female_percent,")")]
  tab1[,Age := paste0(Age_mean, " (",Age_SD,")")]
  tab1[,Statins := paste0(Statin_absolut," (",Statin_percent,")")]
  tab1[,Age_1st_statin := paste0(Age_1st_statin_mean, " (",Age_1st_statin_SD,")")]
  
  # restrict to relevant columns
  tab1 = tab1[,c(1,2,16:21)]
  tab1
}

#' # Get Sup Tab 2 ####
#' ***
#' Selected instruments
#' 
#' Here, I want to give all my summary statistics, that I used throughout the analysis. I want I row per SNP and setting (main, sens with no statin data, sens with no variability, sens with no slope, sens with random effect in sigma model). 
#' 
{
  loaded = load("../results/01_Prep_03_SNPList.RData")
  tab2 = get(loaded)
  
  # reorder the table
  names(tab2)
  tab2 = tab2[,c(1:9,14,10:13,15,16,21,22,17:20,23,24,29,30,25:28)]
  
  if(!file.exists("../results/02_SNPs_05_SENS_RandomEffectSigma.RData")){
    input = list.files(path = "../results/",pattern = "02_SNPs_05")
    dumTab = foreach(i=1:length(input))%do%{
      #i=1
      load(paste0("../results/",input[i]))
      myAssocs_X
    }
    myAssocs_X = rbindlist(dumTab)
    save(myAssocs_X,file = "../results/02_SNPs_05_SENS_RandomEffectSigma.RData")
  }

  mySumStats = list.files(path = "../results/",pattern = "02_SNPs_0")
  mySumStats = mySumStats[!grepl("SNPset",mySumStats)]
  
  settings = gsub(".RData","",mySumStats)
  settings = gsub("02_SNPs_0._","",settings)
  settings[grepl("RandomEffect",settings)] = "SENS_RIsigma"
  
  dumTab = foreach(i = 1:length(mySumStats))%do%{
    #i=1
    loaded2 = load(paste0("../results/",mySumStats[i]))
    tab2.2 = get(loaded2)
    tab2.2 = tab2.2[SNP %in% tab2$rsID,]
    stopifnot(tab2.2$SNP == tab2$rsID)
    names(tab2.2) = paste0("UKB_TC_",names(tab2.2))
    
    tab2.3 = copy(tab2)
    x = dim(tab2.2)[2]
    tab2.3 = cbind(tab2.3[,1:8],tab2.2[,3:x],tab2.3[,9:30])
    tab2.3[,flag := settings[i]]
    tab2.3
  }
  tab2 = rbindlist(dumTab,fill=T)
  tab2 = tab2[,c(44,1:43)]
}

#' # Get Sup Tab 3 ####
#' ***
#' POPS MVMR tables
#' 
{
  myMVMR_results = list.files(path = "../results/",pattern = "03_MVMR_0")
  
  dumTab4 = foreach(i = 1:length(myMVMR_results))%do%{
    #i=1
    load(paste0("../results/",myMVMR_results[i]))
    names(MVMR_results)
    MVMR_results = MVMR_results[setting == "multivariate",]
    MVMR_results[outcome == "UKB", setting := "1-sample"]
    MVMR_results[outcome != "UKB", setting := "2-sample"]
    MVMR_results[outcome != "UKB", outcome := "Aragam_CAD"]
    MVMR_results[outcome == "UKB", outcome := "UKB_CAD"]
    
    MVMR_tab_wide = dcast(MVMR_results, ID + setting + exposure + outcome + threshold + NR_SNPs_total + HeteroStat + HeteroStat_pval ~ exposure_type, 
                          value.var=c("NR_SNPs_type","beta_IVW","SE_IVW","pval_IVW","condFstat"))
    names(MVMR_tab_wide) = gsub("IVW_","",names(MVMR_tab_wide))
    names(MVMR_tab_wide) = gsub("NR_SNPs_type","SNPs",names(MVMR_tab_wide))
    setnames(MVMR_tab_wide,"ID","flag")
    MVMR_tab_wide[,exposure := "TC"]
    
    MVMR_tab_wide
    
  }
  tab3 = rbindlist(dumTab4,fill=T)
  names(tab3)
  tab3[,beta_slope_adj := beta_slope/55.654]
  tab3[,SE_slope_adj := SE_slope/55.654]
  x1 = grep("mean",names(tab3))
  x2 = grep("slope",names(tab3))
  x2 = x2[c(1,2,6,3,7,4,5)]
  x3 = grep("var",names(tab3))
  x = c(1:8,x1,x2,x3)
  tab3 = tab3[,x,with=F]
  tab3[flag == "main",flag:="MAIN"]
  tab3[,flag:=gsub("sens","SENS",flag)]
  tab3[,flag:=gsub("randomEffectSigma","RIsigma",flag)]
  
}  

#' # Get Sup Tab 4 ####
#' ***
#' POPS MR tables
#' 
{
  myMVMR_results = list.files(path = "../results/",pattern = "03_MVMR_0")
  
  dumTab4 = foreach(i = 1:length(myMVMR_results))%do%{
    #i=1
    load(paste0("../results/",myMVMR_results[i]))
    names(MVMR_results)
    MVMR_results = MVMR_results[setting == "univariate",]
    MVMR_results[outcome == "UKB", setting := "1-sample"]
    MVMR_results[outcome != "UKB", setting := "2-sample"]
    MVMR_results[outcome != "UKB", outcome := "Aragam_CAD"]
    MVMR_results[outcome == "UKB", outcome := "UKB_CAD"]
    
    names(MVMR_results) = gsub("_IVW","",names(MVMR_results))
    setnames(MVMR_results,"ID","flag")
    setnames(MVMR_results,"NR_SNPs_total","SNPs")
    MVMR_results = MVMR_results[,c(13,1,2,4,14,3,5,7:12)]
    MVMR_results[,exposure := "TC"]
    MVMR_results
    
  }
  tab4 = rbindlist(dumTab4,fill=T)
  names(tab4)
  tab4.1 = copy(tab4)
  tab4.1 = tab4.1[exposure_type == "slope",]
  tab4.1[,beta := beta/55.654]
  tab4.1[,SE := SE/55.654]
  tab4.1[,exposure_type := "slope_adj",]
  tab4 = rbind(tab4,tab4.1)
  tab4[flag == "main",flag:="MAIN"]
  tab4[,flag:=gsub("sens","SENS",flag)]
  tab4[,flag:=gsub("randomEffectSigma","RIsigma",flag)]
}  

#' # Save tables ###
#' ***

tosave4 = data.table(data = c("tab0", "tab1", "tab2", "tab3", "tab4"), 
                     SheetNames = c("Content","TableS1","TableS2","TableS3","TableS4"))
excel_fn = paste0("../results/_tables/SupplementalTables_realdata_GLGC.xlsx")
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(tab0, tab1, tab2, tab3, tab4, 
     file = paste0("../results/_tables/SupplementalTables_GLGC.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
