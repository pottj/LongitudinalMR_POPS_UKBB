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

source("../../SourceFile.R")

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
  load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_BL_FU_GP_merged_filtered.RData"))
  tab1 = fread("../results/01_Prep_01_summary.txt")
  
  # round the numbers
  tab1[,TC_mean := round(TC_mean,2)]
  tab1[,TC_SD := round(TC_SD,2)]
  tab1[,TC_BL_mean := round(TC_BL_mean,2)]
  tab1[,TC_BL_SD := round(TC_BL_SD,2)]
  tab1[,Female_percent := round(Female_percent,3)*100]
  tab1[,Age_mean := round(Age_mean,2)]
  tab1[,Age_SD := round(Age_SD,2)]
  tab1[,Age_BL_mean := round(Age_BL_mean,2)]
  tab1[,Age_BL_SD := round(Age_BL_SD,2)]
  tab1[,Statin_percent := round(Statin_percent,3)*100]
  tab1[,Age_1st_statin_mean := round(Age_1st_statin_mean,2)]
  tab1[,Age_1st_statin_SD := round(Age_1st_statin_SD,2)]
  
  # merge some columns
  tab1[,measurements := paste0(NR_measurements_median," (",NR_measurements_1stQ,", ",NR_measurements_3rdQ,")")]
  tab1[,TC := paste0(TC_mean," (",TC_SD,")")]
  tab1[,TC_BL := paste0(TC_BL_mean," (",TC_BL_SD,")")]
  tab1[,Female := paste0(Female_absolute," (",Female_percent,")")]
  tab1[,Age := paste0(Age_mean, " (",Age_SD,")")]
  tab1[,Age_BL := paste0(Age_BL_mean, " (",Age_BL_SD,")")]
  tab1[,Statins := paste0(Statin_absolut," (",Statin_percent,")")]
  tab1[,Age_1st_statin := paste0(Age_1st_statin_mean, " (",Age_1st_statin_SD,")")]
  
  # restrict to relevant columns
  tab1 = tab1[,c(2,3,21:28)]
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
  
  mySumStats = list.files(path = "../results/",pattern = "02_SNPs_0")
  mySumStats = mySumStats[!grepl("SNPset",mySumStats)]
  
  settings = gsub(".RData","",mySumStats)
  settings = gsub("02_SNPs_0._","",settings)
  
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
  tab2[, UKB_TC_beta_slope_ageCorrected := UKB_TC_beta_slope * 55.654]
  tab2[, UKB_TC_SE_slope_ageCorrected := UKB_TC_SE_slope * 55.654]
  tab2 = tab2[,c(1:16,45,46,17:44)]
  
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
  x1 = grep("mean",names(tab3))
  x2 = grep("slope",names(tab3))
  x3 = grep("var",names(tab3))
  x = c(1:8,x1,x2,x3)
  tab3 = tab3[,x,with=F]
  tab3[flag == "main",flag:="0 - MAIN"]
  tab3[flag == "sens_noSlope",flag:="1B - SENS - no slope"]
  tab3[flag == "sens_noVar",flag:="1A - SENS - no variability"]
  tab3[flag == "sens_sampleSet1",flag:="2A - SENS - no statins"]
  tab3[flag == "sens_sampleSet2",flag:="2B - SENS - after BL"]
  tab3[flag == "sens_sampleSet3",flag:="2C - SENS - before BL"]
  tab3[flag == "sens_SNPset1",flag:="3A - SENS - GxE enriched SNPs"]
  tab3[flag == "sens_SNPset2",flag:="3B - SENS - mean and var indep."]
  tab3[flag == "sens_SNPset3",flag:="3C - SENS - top 20"]
  
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
  tab4[flag == "main",flag:="0 - MAIN"]
  tab4[flag == "sens_noSlope",flag:="1B - SENS - no slope"]
  tab4[flag == "sens_noVar",flag:="1A - SENS - no variability"]
  tab4[flag == "sens_sampleSet1",flag:="2A - SENS - no statins"]
  tab4[flag == "sens_sampleSet2",flag:="2B - SENS - after BL"]
  tab4[flag == "sens_sampleSet3",flag:="2C - SENS - before BL"]
  tab4[flag == "sens_SNPset1",flag:="3A - SENS - GxE enriched SNPs"]
  tab4[flag == "sens_SNPset2",flag:="3B - SENS - mean and var indep."]
  tab4[flag == "sens_SNPset3",flag:="3C - SENS - top 20"]
  
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
