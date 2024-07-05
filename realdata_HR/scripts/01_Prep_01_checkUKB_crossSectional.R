#' ---
#' title: "Check UK Biobank data - part 1"
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
#' I want to check if the traits I want as outcome for my real data application are available in our UK Biobank project using the heart rate in the ECG bike test as longitudinal exposure. 
#' 
#' Standard columns: 
#' 
#' - sex and age
#' - genetic data available (ethnic background, PCs available, kinship)
#' 
#' ECG columns: 
#' 
#' - cross-sectional ones (there should only be one column for these)
#' - longitudinal ones (either per phase or per trend)
#' 
#' Here, I load the ECG cross-sectional and phase-wise data. 
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load data ####
#' ***
myTab_20 = fread(UKB_phenotypes , 
                 header=TRUE, sep="\t",nrows = 20)
myAnnot = data.table(colNm = names(myTab_20))
myAnnot[,colNR := 1:dim(myTab_20)[2]]

#' # Get standard columns ####
#' *** 
myAnnot[colNm == "f.eid"]
myAnnot[colNm == "f.eid", parameter := "ID"]
myAnnot[colNm == "f.eid", comment := "Encoded anonymised participant ID"]

myAnnot[colNm == "f.31.0.0"]
myAnnot[colNm == "f.31.0.0", parameter := "Sex"]
myAnnot[colNm == "f.31.0.0", comment := "coding females as 0 and males as 1"]

myAnnot[colNm %in% c("f.21022.0.0")]
myAnnot[colNm %in% c("f.21022.0.0"), parameter := "Age"]
myAnnot[colNm %in% c("f.21022.0.0"), comment := c("at recuitment in years")]

myAnnot[colNm %in% c("f.21000.0.0")]
myAnnot[colNm %in% c("f.21000.0.0"), parameter := "Ethnicity"]
myAnnot[colNm %in% c("f.21000.0.0"), comment := c("Ethnic background, 1001 codes for White British")]

myAnnot[colNm %in% c(paste0("f.22009.0.",1:10))]
myAnnot[colNm %in% c(paste0("f.22009.0.",1:10)), parameter := paste0("PC",1:10)]
myAnnot[colNm %in% c(paste0("f.22009.0.",1:10)), comment := c("first 10 (out of 40)")]

myAnnot[colNm %in% c(paste0("f.22021.0.",0))]
myAnnot[colNm %in% c(paste0("f.22021.0.",0)), parameter := "Kinship"]
myAnnot[colNm %in% c(paste0("f.22021.0.",0)), comment := "0 - no kinship; 1 - at least one kinship; 10 - more than 10"]

#' # Get ECG single columns ####
#' ***  
myAnnot[colNm %in% c(paste0("f.6014.0.",0))]
myAnnot[colNm %in% c(paste0("f.6014.0.",0)), parameter := "SQ1"]
myAnnot[colNm %in% c(paste0("f.6014.0.",0)), comment := "Safety question 1: Doctor restricts physical activity due to heart condition"]

myAnnot[colNm %in% c(paste0("f.6015.0.",0))]
myAnnot[colNm %in% c(paste0("f.6015.0.",0)), parameter := "SQ2"]
myAnnot[colNm %in% c(paste0("f.6015.0.",0)), comment := "Safety question 2: Chest pain felt during physical activity"]

myAnnot[colNm %in% c(paste0("f.6016.0.",0))]
myAnnot[colNm %in% c(paste0("f.6016.0.",0)), parameter := "SQ3"]
myAnnot[colNm %in% c(paste0("f.6016.0.",0)), comment := "Safety question 3: Chest pain felt outside physical activity"]

myAnnot[colNm %in% c(paste0("f.6017.0.",0))]
myAnnot[colNm %in% c(paste0("f.6017.0.",0)), parameter := "SQ4"]
myAnnot[colNm %in% c(paste0("f.6017.0.",0)), comment := "Safety question 4: Able to walk or cycle unaided for 10 minutes"]

myAnnot[colNm %in% c(paste0("f.6019.0.",0))]
myAnnot[colNm %in% c(paste0("f.6019.0.",0)), parameter := "Method"]
myAnnot[colNm %in% c(paste0("f.6019.0.",0)), comment := "ECG/bike method for fitness test"]

myAnnot[colNm %in% c(paste0("f.6020.0.",0))]
myAnnot[colNm %in% c(paste0("f.6020.0.",0)), parameter := "Status"]
myAnnot[colNm %in% c(paste0("f.6020.0.",0)), comment := "Completion status of test"]

myAnnot[colNm %in% c(paste0("f.6023.0.",0))]
myAnnot[colNm %in% c(paste0("f.6023.0.",0)), parameter := "Protocol"]
myAnnot[colNm %in% c(paste0("f.6023.0.",0)), comment := "Description of exercise protocol recommended"]

myAnnot[colNm %in% c(paste0("f.6024.0.",0))]

#' This is bad, we do not have the information what category the proband was in. But I will be able to extract that from the protocol data entries. 
#' 
myAnnot[colNm %in% c(paste0("f.6032.0.",0))]
myAnnot[colNm %in% c(paste0("f.6032.0.",0)), parameter := "MaxLoad"]
myAnnot[colNm %in% c(paste0("f.6032.0.",0)), comment := "Maximum workload during fitness test"]

myAnnot[colNm %in% c(paste0("f.6033.0.",0))]
myAnnot[colNm %in% c(paste0("f.6033.0.",0)), parameter := "MaxHR"]
myAnnot[colNm %in% c(paste0("f.6033.0.",0)), comment := "Maximum heart rate during fitness test"]

myAnnot[colNm %in% c(paste0("f.6034.0.",0))]
myAnnot[colNm %in% c(paste0("f.6034.0.",0)), parameter := "TargetHR"]
myAnnot[colNm %in% c(paste0("f.6034.0.",0)), comment := "Target heart rate achieved"]

myAnnot[colNm %in% c(paste0("f.6038.0.",0))]
myAnnot[colNm %in% c(paste0("f.6038.0.",0)), parameter := "NR_Trend"]
myAnnot[colNm %in% c(paste0("f.6038.0.",0)), comment := "Number of trend entries"]

myAnnot[colNm %in% c(paste0("f.6039.0.",0))]
myAnnot[colNm %in% c(paste0("f.6039.0.",0)), parameter := "Duration"]
myAnnot[colNm %in% c(paste0("f.6039.0.",0)), comment := "Duration of fitness test"]

#' # Get ECG multiple columns ####
#' *** 
#' In this script, I just want to check the tree data fields with information per phase (3 phases).
#' 
myAnnot[grepl("f.5991.0",colNm),]
myAnnot[grepl("f.5991.0",colNm), parameter := paste0("Phase",1:3)]
myAnnot[grepl("f.5991.0",colNm), comment := "This field contains phase name"]

myAnnot[grepl("f.5992.0",colNm),]
myAnnot[grepl("f.5992.0",colNm), parameter := paste0("Phase",1:3,"Duration")]
myAnnot[grepl("f.5992.0",colNm), comment := "This field contains phase name"]

myAnnot[grepl("f.5993.0",colNm),]
myAnnot[grepl("f.5993.0",colNm), parameter := paste0("Phase",1:3,"NR_Stages")]
myAnnot[grepl("f.5993.0",colNm), comment := "This field shows total number of stages within a phase."]


#' # Load all samples ####
#' ***
x = myAnnot[!is.na(parameter),colNR]
myTab_cross <- fread(UKB_phenotypes , 
                     header=TRUE, sep="\t",select = x)
names(myTab_cross)

matched = match(names(myTab_cross),myAnnot$colNm)
myAnnot[matched,parameter]
names(myTab_cross) = myAnnot[matched,parameter]

#' # Filter samples ####
#' ***
#' Okay, now I filter my data set for 
#' 
#' - white British people with no kinship (genetic data available)
#' - ECG bike data available
#' 
myTab_cross = myTab_cross[Ethnicity == 1001,]
myTab_cross = myTab_cross[Kinship == 0,]
myTab_cross[,table(is.na(PC1))]

#' Okay, 292,154 white British people with no kinship (genetic data available)
#' 
#' Let's exclude all people flag during the safety questions (I want "no" as answers for SQ1, SQ2, and SW3, coded as 0, and "yes for SQ4, coded as 1)
myTab_cross[,table(is.na(SQ1))]
myTab_cross = myTab_cross[SQ1 == 0,]
myTab_cross = myTab_cross[SQ2 == 0,]
myTab_cross = myTab_cross[SQ3 == 0,]
myTab_cross = myTab_cross[SQ4 == 1,]

#' Okay, 39,268 people passed the safety questions. Let's exclude all people who did not do the bicycle (vs resting, bicycle is coded as 1 in **Method**) and everyone who stopped early (**Status**: 1 = "Fully completed"; 31 = "Participant wanted to stop early"; 32 = "Participant reported chest-pain and/or other discomfort"; 33 = "Heart rate reached safety level"; 34 = "Incomplete - other reason")
myTab_cross = myTab_cross[Method == 1,]
myTab_cross[,table(Status)]
myTab_cross = myTab_cross[Status == 1,]

#' Okay, 34,600 people finished the bicycle test. Now I want to check if I can get those with minimal risk (cycle at 50% level of the absolute maximum work load according to **Protocol**)
#' 
myTab_cross[1:10,Protocol]
myTab_cross = myTab_cross[!is.na(Protocol),]
myTab_cross = myTab_cross[grepl("50%",Protocol),]

#' I see the following information
#' 
#' - sex of participant (M=male, F=female)
#' - target load (I think)
#' - protocol ID (I think)
#' - warm-up: time and load
#' - ramp-up: time and target load, and maximum load
#' - cool-down: time
#' 
myTab_cross[,sex_ECG := substr(Protocol,1,1)]
myTab_cross[,table(sex_ECG,Sex)]

myTab_cross[,targetLoad := gsub("/.*","",Protocol)]
myTab_cross[,targetLoad := gsub("[M,F]","",targetLoad)]
myTab_cross[,table(targetLoad)]
myTab_cross[,targetLoad := as.numeric(targetLoad)]

myTab_cross[,protocolID := gsub(" :.*","",Protocol)]
myTab_cross[,protocolID := gsub(".*/","",protocolID)]
myTab_cross[,protocolID := as.numeric(protocolID)]

myTab_cross[,warmUp_duration := gsub(".* : ","",Protocol)]
myTab_cross[,warmUp_duration := gsub(",.*","",warmUp_duration)]
myTab_cross[,table(warmUp_duration,sex_ECG)]
myTab_cross[,warmUp_load := 40]
myTab_cross[sex_ECG == "F",warmUp_load := 30]
myTab_cross[,warmUp_duration := 2]
myTab_cross[,warmUp_duration := as.numeric(warmUp_duration)]

myTab_cross[,rampUp := gsub(".*W, ","",Protocol)]
myTab_cross[,rampUp := gsub(",.*","",rampUp)]
myTab_cross[,table(rampUp,sex_ECG)]
myTab_cross[,rampUp_duration := substr(rampUp,1,1)]
myTab_cross[,rampUp_targetLoad := gsub("4 mins ramp to ","",rampUp)]
myTab_cross[,rampUp_targetLoad := gsub("W.*","",rampUp_targetLoad)]
myTab_cross[,rampUp_targetLoad := as.numeric(rampUp_targetLoad)]
myTab_cross[,rampUp_maxLoad := gsub(".* of ","",rampUp)]
myTab_cross[,rampUp_maxLoad := gsub("W.*","",rampUp_maxLoad)]
myTab_cross[,rampUp_maxLoad := as.numeric(rampUp_maxLoad)]
myTab_cross[,rampUp := NULL]

myTab_cross[,coolDown := gsub(".*, ","",Protocol)]
myTab_cross[,table(coolDown)]
myTab_cross[,coolDown_duration := 1]
myTab_cross[,coolDown := NULL]

#' Okay, 30,750 people finished the bicycle test with 50% ramp-up. There are still some things I am not sure of: 
#' 
#' - Duration: the duration should be approx 420 seconds ((2 mins warm-up + 4 mins ramp-up + 1 min cool-down)*60 = 420), but there are people with less then 300 minutes - how could they have successfully completed the test?
#' - NR_Trend: again, how can there be cases with just a handful of trends?
#' - maxHR: over 200 seems a bit of the top... 
#' - maxLoad: I do not trust 0 (might be typo), but I cannot trust this. 
#' 
#' So I decide to exclude people with Duration <300, maxHR>200, and maxLoad==0
#' 
myTab_cross = myTab_cross[Duration>=300 & Duration <= 500 & MaxHR<=200 & MaxLoad !=0,]

#' Okay, 30,531 people finished the bicycle test with 50% ramp-up and had summary results in normal ranges. Let's see what information is within the phases.
myTab_cross[,table(Phase1)]
myTab_cross[,table(Phase2)]
myTab_cross[,table(Phase3)]

#' It is the same information for all samples. I decide to remove those samples with phase name *MANUAL* (something must have gone wrong here). 
myTab_cross = myTab_cross[Phase2!="MANUAL" & Phase3!="MANUAL",]

#' Okay, 30,523 people finished the bicycle test with 50% ramp-up and had summary results in normal ranges and normal phase names. Check duration of phases and number of stages within this phase
myTab_cross[,table(Phase1Duration,Phase1NR_Stages)]
myTab_cross[,table(Phase3Duration,Phase3NR_Stages)]

ggplot(myTab_cross, aes(x=Phase2Duration, y=Phase2NR_Stages, col = as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="duration of exercise", y="number of stages in exercise")

#' I will restrict my sample set to people with duration of 15 seconds in phase 1 (1 stage) and duration of 60 seconds in phase 3 (1 stage). Phase 2 looks ok (the more intense the protocol, the more stages there are). 
#'    
myTab_cross = myTab_cross[Phase1Duration == 15 & Phase1NR_Stages == 1,]
myTab_cross = myTab_cross[Phase3Duration == 60 & Phase3NR_Stages == 1,]

#' # Get some plots ####
#' ***
#' Okay, now I have my final set of 30,307 people. Let's  get some plots! (all seperated by sex)
#' 
ggplot(myTab_cross, aes(x=MaxLoad, fill=as.factor(Sex))) +
  geom_histogram(position="identity", colour="grey40", alpha=0.5, bins = 10) +
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="maximum workload during fitness test (watt)", y="Counts")

ggplot(myTab_cross, aes(x=MaxHR, fill=as.factor(Sex))) +
  geom_histogram(position="identity", colour="grey40", alpha=0.5, bins = 10) +
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="maximum heart rate during fitness test (beats per minute)", y="Counts")

ggplot(myTab_cross, aes(x=Duration, fill=as.factor(Sex))) +
  geom_histogram(position="identity", colour="grey40", alpha=0.5, bins = 10) +
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="duration of fitness test (in seconds)", y="Counts")

ggplot(myTab_cross, aes(x=NR_Trend, fill=as.factor(Sex))) +
  geom_histogram(position="identity", colour="grey40", alpha=0.5, bins = 10) +
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="number of trend entries", y="Counts")

ggplot(myTab_cross, aes(x=NR_Trend, y=Duration, col= as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="number of trend entries", y="duration of fitness test (in seconds)")

ggplot(myTab_cross, aes(x=MaxHR, y=MaxLoad, col= as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="maximum heart rate (beats per minute)", y="maximum workload (watt)")

ggplot(myTab_cross, aes(x=targetLoad, y=MaxLoad, col= as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="target workload (according to protocol)", y="maximum workload (watt, as measured)")

ggplot(myTab_cross, aes(x=as.factor(protocolID), y=MaxLoad)) + 
  geom_boxplot() +
  facet_grid(. ~ Sex,scales = "free",
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="Protocol ID", y="maximum workload (watt, as measured)")

ggplot(myTab_cross, aes(x=rampUp_targetLoad, y=rampUp_maxLoad, col= as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="target workload (according to ramp-up phase)", y="maximum workload (watt, as expected)")

ggplot(myTab_cross, aes(x=rampUp_maxLoad, y=MaxLoad, col= as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="maximum workload (watt, as expected)", y="maximum workload (watt, as measured)")

ggplot(myTab_cross, aes(x=NR_Trend, y=Phase2NR_Stages, col = as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="number of trends in exercise", y="number of stages in exercise")

#' # Save data ####
#' ***
save(myTab_cross,file = paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_HR_cross.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
