#' ---
#' title: "Get Supplemental Tables for real data"
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
#' **Supplemental Tables of real data analyses**
#' 
#' 1) POPS - Cohort description (mean and SD of all covariables used in the regression model - restricted to setting main setting)
#' 2) POPS - Overview of instruments
#' 3) POPS - MVMR results (main + sensitivity)
#' 4) POPS - MR results (main + sensitivity)
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
                    Title = c("POPS - Study description",
                              "POPS - SNP summary statistics",
                              "POPS - MVMR results",
                              "POPS - MR results"),
                    Source = c("done in this script",
                               "realdata_GLGC/results/02_SNPs_*.RData",
                               "realdata_GLGC/results/04_MVMR_*.RData",
                               "realdata_GLGC/results/04_MVMR_*.RData"))
  
  tab0
  
}

#' # Get Sup Tab 1 ####
#' ***
#' Study description
#' 
#' I want for the exposure (MAIN and GBR_sampleSet): 
#' 
#' - total sample size
#' - sample size per time point
#' - mean (SD) of EFW / BW measurements per time point
#' - mean (SD) of gestational age per time point
#' - Female % 
#' - Smoking status mother % 
#' - mean (SD) of mothers height 
#' 
{
  myGD_files = list.files(path = POPS_phenotypes,pattern = "01_Prep_04")
  loaded1 = load(paste0(POPS_phenotypes,myGD_files[1]))
  loaded1
  loaded2 = load(paste0(POPS_phenotypes,myGD_files[2]))
  loaded2
  
  # filter for relevant samples 
  myTab_X[,logefwcomb := log(efwcomb)]
  myTab_X = myTab_X[!is.na(logefwcomb)]
  test1 = myTab_X[,.N,by=POPSID]
  myTab_X = myTab_X[POPSID %in% test1[N>1,POPSID]]
  myTab_Y = myTab_Y[POPSID %in% myTab_X$POPSID,]
  
  # get table for time point birth 
  tab1a = data.table(timpoint = "time-fixed", 
                     sampleSize = dim(myTab_Y)[1],
                     EFW_mean = mean(myTab_Y$pn_bw),
                     EFW_SD = sd(myTab_Y$pn_bw),
                     EFW_log_mean = mean(log(myTab_Y$pn_bw)),
                     EFW_log_SD = sd(log(myTab_Y$pn_bw)),
                     Age_mean = mean(myTab_Y$pn_ga_wk),
                     Age_SD = sd(myTab_Y$pn_ga_wk),
                     Female_absolute = dim(myTab_Y[pn_sex=="FEMALE",])[1],
                     Female_percent = dim(myTab_Y[pn_sex=="FEMALE",])[1]/dim(myTab_Y)[1],
                     Smoking_status_never = dim(myTab_Y[an_smokstat=="Never smoked",])[1],
                     Smoking_status_Qpre = dim(myTab_Y[an_smokstat=="Quit (prepregnancy)",])[1],
                     Smoking_status_Qduring = dim(myTab_Y[an_smokstat=="Quit (During pregnancy)",])[1],
                     Smoking_status_current = dim(myTab_Y[an_smokstat=="Currently smoking",])[1],
                     Smoking_status_never_percent = dim(myTab_Y[an_smokstat=="Never smoked",])[1]/dim(myTab_Y)[1],
                     Smoking_status_Qpre_percent = dim(myTab_Y[an_smokstat=="Quit (prepregnancy)",])[1]/dim(myTab_Y)[1],
                     Smoking_status_Qduring_percent = dim(myTab_Y[an_smokstat=="Quit (During pregnancy)",])[1]/dim(myTab_Y)[1],
                     Smoking_status_current_percent = dim(myTab_Y[an_smokstat=="Currently smoking",])[1]/dim(myTab_Y)[1],
                     Heigth_mean = mean(myTab_Y$an_height),
                     Heigth_SD = sd(myTab_Y$an_height) )
  
  # get table for first scan 
  tab1b = data.table(timpoint = "scan 1", 
                     sampleSize = dim(myTab_X[scan==1])[1],
                     EFW_mean = myTab_X[scan==1,mean(efwcomb)],
                     EFW_SD = myTab_X[scan==1,sd(efwcomb)],
                     EFW_log_mean = myTab_X[scan==1,mean(logefwcomb)],
                     EFW_log_SD = myTab_X[scan==1,sd(logefwcomb)],
                     Age_mean = myTab_X[scan==1,mean(ga)],
                     Age_SD = myTab_X[scan==1,sd(ga)] )
  
  # get table for second scan 
  tab1c = data.table(timpoint = "scan 2", 
                     sampleSize = dim(myTab_X[scan==2])[1],
                     EFW_mean = myTab_X[scan==2,mean(efwcomb)],
                     EFW_SD = myTab_X[scan==2,sd(efwcomb)],
                     EFW_log_mean = myTab_X[scan==2,mean(logefwcomb)],
                     EFW_log_SD = myTab_X[scan==2,sd(logefwcomb)],
                     Age_mean = myTab_X[scan==2,mean(ga)],
                     Age_SD = myTab_X[scan==2,sd(ga)] )
  
  # get table for third scan 
  tab1d = data.table(timpoint = "scan 3", 
                     sampleSize = dim(myTab_X[scan==3])[1],
                     EFW_mean = myTab_X[scan==3,mean(efwcomb)],
                     EFW_SD = myTab_X[scan==3,sd(efwcomb)],
                     EFW_log_mean = myTab_X[scan==3,mean(logefwcomb)],
                     EFW_log_SD = myTab_X[scan==3,sd(logefwcomb)],
                     Age_mean = myTab_X[scan==3,mean(ga)],
                     Age_SD = myTab_X[scan==3,sd(ga)] )
  
  tab1 = rbind(tab1a,tab1b,tab1c,tab1d,fill=T)
  
  # round the numbers
  tab1[,EFW_mean := round(EFW_mean,1)]
  tab1[,EFW_SD := round(EFW_SD,1)]
  tab1[,EFW_log_mean := round(EFW_log_mean,2)]
  tab1[,EFW_log_SD := round(EFW_log_SD,2)]
  tab1[,Age_mean := round(Age_mean,1)]
  tab1[,Age_SD := round(Age_SD,1)]
  tab1[,Female_percent := round(Female_percent,3)*100]
  tab1[,Smoking_status_never_percent := round(Smoking_status_never_percent,3)*100]
  tab1[,Smoking_status_Qpre_percent := round(Smoking_status_Qpre_percent,3)*100]
  tab1[,Smoking_status_Qduring_percent := round(Smoking_status_Qduring_percent,3)*100]
  tab1[,Smoking_status_current_percent := round(Smoking_status_current_percent,3)*100]
  tab1[,Heigth_mean := round(Heigth_mean,1)]
  tab1[,Heigth_SD := round(Heigth_SD,1)]
  
  # merge some columns
  tab1[,EFW := paste0(EFW_mean," (",EFW_SD,")")]
  tab1[,logEFW := paste0(EFW_log_mean," (",EFW_log_SD,")")]
  tab1[,Age := paste0(Age_mean, " (",Age_SD,")")]
  tab1[,Female := paste0(Female_absolute," (",Female_percent,")")]
  tab1[,M_Smoking0 := paste0(Smoking_status_never," (",Smoking_status_never_percent,")")]
  tab1[,M_Smoking1 := paste0(Smoking_status_Qpre," (",Smoking_status_Qpre_percent,")")]
  tab1[,M_Smoking2 := paste0(Smoking_status_Qduring," (",Smoking_status_Qduring_percent,")")]
  tab1[,M_Smoking3 := paste0(Smoking_status_current," (",Smoking_status_current_percent,")")]
  tab1[,M_Height := paste0(Heigth_mean, " (",Heigth_SD,")")]
  
  # restrict to relevant columns
  tab1 = tab1[,c(1,2,21:29)]
  tab1[timpoint!="time-fixed",Female := ""]
  tab1[timpoint!="time-fixed",M_Smoking0 := ""]
  tab1[timpoint!="time-fixed",M_Smoking1 := ""]
  tab1[timpoint!="time-fixed",M_Smoking2 := ""]
  tab1[timpoint!="time-fixed",M_Smoking3 := ""]
  tab1[timpoint!="time-fixed",M_Height := ""]
  tab1
}

#' # Get Sup Tab 2 ####
#' ***
#' Selected instruments
#' 
#' Here, I want to give all my summary statistics, that I used throughout the analysis. I want I row per SNP and setting (main, sens with no statin data, sens with no variability, sens with no slope, sens with random effect in sigma model). 
#' 
#' Note: I still have to add the summary statistics from the POPS BW GWAS 
{
  loaded = load("../results/01_Prep_05_SNPList.RData")
  tab2 = get(loaded)
  
  # reorder the table
  names(tab2)
  tab2 = tab2[,c(1:9,14,10:13,15,16,21,17:20)]
  
  mySumStats = list.files(path = "../results/",pattern = "02_SNPs_0")
  #mySumStats = mySumStats[c(1:4)]
  
  settings = gsub(".RData","",mySumStats)
  settings = gsub("02_SNPs_0._","",settings)
  settings[grepl("RandomEffect",settings)] = "SENS_RIsigma"
  
  dumTab = foreach(i = 1:length(mySumStats))%do%{
    #i=1
    loaded2 = load(paste0("../results/",mySumStats[i]))
    tab2.2 = get(loaded2)
    tab2.2 = tab2.2[rsID %in% tab2$rsID,]
    tab2.2 = tab2.2[phenotype == "logefwcomb",]
    stopifnot(tab2.2$rsID == tab2$rsID)
    names(tab2.2) = paste0("POPS_EFW_",names(tab2.2))
    
    tab2.3 = copy(tab2)
    x = dim(tab2.2)[2]
    tab2.3 = cbind(tab2.3[,1:8],tab2.2[,8:x],tab2.3[,9:21])
    tab2.3[,flag := settings[i]]
    
    load("../results/03_SNPs_01_MAIN_Assocs_outcome_nTIA.RData")
    if(settings[i]!="SENS_GBR3"){
      myAssocs_Y = myAssocs_Y[population == "all",]
    }else{
      myAssocs_Y = myAssocs_Y[population == "GBR",]
    }
    myAssocs_Y = myAssocs_Y[rsID %in% tab2.3$rsID,]
    myAssocs_Y = myAssocs_Y[phenotype == "pn_bw",]
    stopifnot(myAssocs_Y$rsID == tab2.3$rsID)
    
    tab2.3[,POPS_beta := myAssocs_Y$beta_mean]
    tab2.3[,POPS_SD := myAssocs_Y$SE_mean]
    tab2.3[,POPS_tval := myAssocs_Y$tval_mean]
    tab2.3[,POPS_pval := myAssocs_Y$pval_mean]
    tab2.3
  }
  tab2 = rbindlist(dumTab,fill=T)
  tab2 = tab2[,c(35,1:34,36:39)]
  tab2[, POPS_EFW_beta_slope_ageCorrected := POPS_EFW_beta_slope * 40.316]
  tab2[, POPS_EFW_SE_slope_ageCorrected := POPS_EFW_SE_slope * 40.316]
  tab2 = tab2[,c(1:16,40,41,17:39)]
  
}

#' # Get Sup Tab 3 ####
#' ***
#' POPS MVMR tables
#' 
{
  myMVMR_results = list.files(path = "../results/",pattern = "04_MVMR_0")
  
  dumTab4 = foreach(i = 1:length(myMVMR_results))%do%{
    #i=1
    load(paste0("../results/",myMVMR_results[i]))
    names(MVMR_results)
    MVMR_results = MVMR_results[setting == "multivariate",]
    MVMR_results = MVMR_results[exposure == "logefwcomb",]
    MVMR_results = MVMR_results[outcome != "EGG_BW",]
    MVMR_results[outcome == "POPS_BW", setting := "1-sample"]
    MVMR_results[outcome != "POPS_BW", setting := "2-sample"]
    
    MVMR_tab_wide = dcast(MVMR_results, ID + setting + exposure + outcome + threshold + NR_SNPs_total + HeteroStat + HeteroStat_pval ~ exposure_type, 
                          value.var=c("NR_SNPs_type","beta_IVW","SE_IVW","pval_IVW","condFstat"))
    names(MVMR_tab_wide) = gsub("IVW_","",names(MVMR_tab_wide))
    names(MVMR_tab_wide) = gsub("NR_SNPs_type","SNPs",names(MVMR_tab_wide))
    setnames(MVMR_tab_wide,"ID","flag")
    MVMR_tab_wide[,exposure := "logEFW"]
    
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
  tab3[flag == "sens_GBR",flag:="2 - SENS - GBR3"]
  tab3[flag == "sens_noSlope",flag:="1B - SENS - no slope"]
  tab3[flag == "sens_noVar",flag:="1A - SENS - no variability"]

}  

#' # Get Sup Tab 4 ####
#' ***
#' POPS MR tables
#' 
{
  myMVMR_results = list.files(path = "../results/",pattern = "04_MVMR_0")
  
  dumTab4 = foreach(i = 1:length(myMVMR_results))%do%{
    #i=1
    load(paste0("../results/",myMVMR_results[i]))
    names(MVMR_results)
    MVMR_results = MVMR_results[setting == "univariate",]
    MVMR_results = MVMR_results[exposure == "logefwcomb",]
    MVMR_results = MVMR_results[outcome != "EGG_BW",]
    MVMR_results[outcome == "POPS_BW", setting := "1-sample"]
    MVMR_results[outcome != "POPS_BW", setting := "2-sample"]
    
    names(MVMR_results) = gsub("_IVW","",names(MVMR_results))
    setnames(MVMR_results,"ID","flag")
    setnames(MVMR_results,"NR_SNPs_total","SNPs")
    MVMR_results = MVMR_results[,c(13,1,2,4,14,3,5,7:12)]
    MVMR_results[,exposure := "logEFW"]
    MVMR_results
    
  }
  tab4 = rbindlist(dumTab4,fill=T)
  names(tab4)
  tab4[flag == "main",flag:="0 - MAIN"]
  tab4[flag == "sens_GBR",flag:="2 - SENS - GBR3"]
  tab4[flag == "sens_noSlope",flag:="1B - SENS - no slope"]
  tab4[flag == "sens_noVar",flag:="1A - SENS - no variability"]
  
}  

#' # Save tables ###
#' ***

tosave4 = data.table(data = c("tab0", "tab1", "tab2", "tab3", "tab4"), 
                     SheetNames = c("Content","TableS1","TableS2","TableS3","TableS4"))
excel_fn = paste0("../results/_tables/SupplementalTables_realdata_EGG.xlsx")
WriteXLS(tosave4$data, 
         ExcelFileName=excel_fn, 
         SheetNames=tosave4$SheetNames, 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(tab0, tab1, tab2, tab3, tab4, 
     file = paste0("../results/_tables/SupplementalTables_EGG.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
