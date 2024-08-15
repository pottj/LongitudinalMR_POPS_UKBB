#' ---
#' title: "Get Sensitivity Forest Plots"
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
#' **Forest Plots with sensitivity data**
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

#' # Get data ####
#' ***
load("../results/04_MVMR_01_MAIN_AF.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
filt0 = grepl("_combined",MVMR_results$exposure) & grepl("_combined",MVMR_results$outcome)
filt1 = grepl("_men",MVMR_results$exposure) & !grepl("_females",MVMR_results$outcome)
filt2 = grepl("_women",MVMR_results$exposure) & !grepl("_males",MVMR_results$outcome)
filt = filt0 | filt1 | filt2
MVMR_results = MVMR_results[filt,]
MVMR_main = copy(MVMR_results)

load("../results/04_MVMR_03_SENS_noVar_AF.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
filt0 = grepl("_combined",MVMR_results$exposure) & grepl("_combined",MVMR_results$outcome)
filt1 = grepl("_men",MVMR_results$exposure) & !grepl("_females",MVMR_results$outcome)
filt2 = grepl("_women",MVMR_results$exposure) & !grepl("_males",MVMR_results$outcome)
filt = filt0 | filt1 | filt2
MVMR_results = MVMR_results[filt,]
MVMR_noVar = copy(MVMR_results)
MVMR_noVar[,ID := "sens1_noVar"]

load("../results/04_MVMR_04_SENS_noSlope_AF.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
filt0 = grepl("_combined",MVMR_results$exposure) & grepl("_combined",MVMR_results$outcome)
filt1 = grepl("_men",MVMR_results$exposure) & !grepl("_females",MVMR_results$outcome)
filt2 = grepl("_women",MVMR_results$exposure) & !grepl("_males",MVMR_results$outcome)
filt = filt0 | filt1 | filt2
MVMR_results = MVMR_results[filt,]
MVMR_noSlope = copy(MVMR_results)
MVMR_noSlope[,ID := "sens2_noSlope"]

load("../results/04_MVMR_02_SENS_rampUp.RData")
MVMR_results = MVMR_results[setting == "multivariate"]
filt0 = grepl("_combined",MVMR_results$exposure) & grepl("_combined",MVMR_results$outcome)
filt1 = grepl("_men",MVMR_results$exposure) & !grepl("_females",MVMR_results$outcome)
filt2 = grepl("_women",MVMR_results$exposure) & !grepl("_males",MVMR_results$outcome)
filt = filt0 | filt1 | filt2
MVMR_results = MVMR_results[filt,]
MVMR_rampUp = copy(MVMR_results)
MVMR_rampUp[,ID := "sens3_rampUp"]

MVMR = rbind(MVMR_main,MVMR_noVar,MVMR_noSlope,MVMR_rampUp)
MVMR[,table(ID)]

#' # Create Forest Plots ####
#' ***
#' I need to do this separately for each phenotype and each exposure type. Let's try a loop... 
#' 
mySex = c("combined","men","women")
myExposure_Types = unique(MVMR$exposure_type)

dumTab0 = foreach(k = 1:length(mySex))%do%{
  #k=1
  message("Working on sex: ",mySex[k])
  MVMR0 = copy(MVMR)
  MVMR0 = MVMR0[grepl(paste0("_",mySex[k]),exposure),]
  myOutcomes = unique(MVMR0$outcome)
  
  dumTab1 = foreach(i = 1:length(myOutcomes))%do%{
    #i=2
    message("   and outcome ",myOutcomes[i])
    MVMR1 = copy(MVMR0)
    MVMR1 = MVMR1[outcome == myOutcomes[i]]
    
    dumTab2 = foreach(j = 1:length(myExposure_Types))%do%{
      #j=1
      myExposure_Type = myExposure_Types[j]
      message("         and type ",myExposure_Type)
      
      MVMR2 = copy(MVMR1)
      MVMR2 = MVMR2[exposure_type == myExposure_Type]
      
      # start plotting 
      MVMR2[,rank := 2]
      dummy = data.table(ID = unique(MVMR2$ID),rank=1)
      data4 = rbind(MVMR2,dummy, fill=T)
      data4[,subgroup := paste0("   ", threshold)]
      data4[is.na(threshold),subgroup := ID]
      data4[,subgroup := gsub("_", " - ",subgroup)]
      setnames(data4,"subgroup", "   Instrument selection")
      
      data4[,lowerCI95 := beta_IVW-1.96*SE_IVW]
      data4[,upperCI95 := beta_IVW+1.96*SE_IVW]
      data4$` ` <- paste(rep(" ", 50), collapse = " ")
      data4$`Estimate [95% CI]` <- ifelse(is.na(data4$SE_IVW), "",
                                          sprintf("%.2f [%.2f, %.2f]",
                                                  data4$beta_IVW, data4$lowerCI95, data4$upperCI95))
      
      min_data4 = signif(min(data4$lowerCI95,na.rm = T),2)
      max_data4 = signif(max(data4$upperCI95,na.rm = T),2)
      min_data4 = min_data4 + 0.2*min_data4
      max_data4 = max_data4 + 0.2*max_data4
      # range = max_data4 - min_data4
      # range = round(range,0)
      # myTicks = seq(min_data4, max_data4, by = round(range/5,2))
      
      setorder(data4,ID,rank)
      # data4 = data4[c(15:21,8:14,1:7,22:28),]
      #data4[threshold == "nominal",threshold := "p<0.05, all SNPs"]
      
      dummy = data4$threshold
      dummy[is.na(dummy)] = "white"
      dummy[grepl("all_SNPs",dummy)] = "lightgrey"
      dummy[grepl("nominal_SNPs",dummy)] = "lightsteelblue2"
      dummy = gsub("top20_overlap","lightgreen",dummy)
      dummy = gsub("top20_distinct","lightcoral",dummy)
      tm1<- forest_theme(core=list(bg_params=list(fill = dummy)))
      
      if(grepl("UKB",myOutcomes[i])){
        myXlab = paste0(myExposure_Type, " of HR on AF (1-sample, UKB, sex: ",mySex[k],")")
      }else{
        myXlab = paste0(myExposure_Type, " of HR on AF (2-sample, Miyazawa, sex: ",mySex[k],")")
      }
      
      p2<- forest(data4[,c(16,19,20)],
                  est = data4$beta_IVW,
                  lower = data4$lowerCI95, 
                  upper = data4$upperCI95,
                  sizes = 0.5,
                  #ticks_at = myTicks,ticks_digits = 2,
                  #ticks_at = c(0,1,2,3),ticks_digits = 2,
                  ci_column = 2,
                  ref_line = 0,
                  xlim = c(min_data4, max_data4),
                  #xlim = c(-0.5, 3.5),
                  title = myXlab,
                  theme = tm1)
      
      plot(p2)
      
      filename = paste0("../results/_figures/05_ForestPlots1_AF/HR_",mySex[k],"_",myOutcomes[i],"_",myExposure_Type,".png")
      png(filename = filename,width = 2000, height = 1200, res=200)
      plot(p2)
      dev.off()
      
      
    }
    
  }
  
}
  

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))
