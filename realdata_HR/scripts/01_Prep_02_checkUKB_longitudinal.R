#' ---
#' title: "Check UK Biobank data - part 2"
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
#' Here, I load the ECG trend-wise data. 
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load data ####
#' ***
myTab_20 = fread(UKB_phenotypes, 
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

#' # Get longitudinal columns ####
#' *** 
#' In this script, I just want to check the tree data fields with information per trend (up to 114 per person).
#' 
myAnnot[grepl("f.5983.0",colNm),]
myAnnot[grepl("f.5983.0",colNm), parameter := paste0("HR_trend",1:112)]
myAnnot[grepl("f.5983.0",colNm), comment := "This field contains heart rate (beats per minute) at the time of the trend entry."]

myAnnot[grepl("f.5984.0",colNm),]
myAnnot[grepl("f.5984.0",colNm), parameter := paste0("Load_trend",1:114)]
myAnnot[grepl("f.5984.0",colNm), comment := "This field contains load in Watts at the time of the trend entry."]

myAnnot[grepl("f.5985.0",colNm),]
myAnnot[grepl("f.5985.0",colNm), parameter := paste0("Speed_trend",1:114)]
myAnnot[grepl("f.5985.0",colNm), comment := "This field contains revolutions per minute of the treadmill at the time of the trend entry."]

myAnnot[grepl("f.5986.0",colNm),]
myAnnot[grepl("f.5986.0",colNm), parameter := paste0("PhaseTime_trend",1:114)]
myAnnot[grepl("f.5986.0",colNm), comment := "This field contains the time spent within the phase of the trend entry."]

myAnnot[grepl("f.5987.0",colNm),]
myAnnot[grepl("f.5987.0",colNm), parameter := paste0("PhaseName_trend",1:114)]
myAnnot[grepl("f.5987.0",colNm), comment := "This field contains the name of the phase when trend entry is taken."]

myAnnot[grepl("f.5988.0",colNm),]
myAnnot[grepl("f.5988.0",colNm), parameter := paste0("StageName_trend",1:114)]
myAnnot[grepl("f.5988.0",colNm), comment := "This field contains the name of the stage when trend entry is taken."]

#' # Load all samples ####
#' ***
x = myAnnot[!is.na(parameter),colNR]
myTab_long <- fread(UKB_phenotypes, 
                     header=TRUE, sep="\t",select = x)
names(myTab_long)

matched = match(names(myTab_long),myAnnot$colNm)
myAnnot[matched,parameter]
names(myTab_long) = myAnnot[matched,parameter]

#' # Filter samples ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_HR_cross.RData"))
myTab_long = myTab_long[ID %in% myTab_cross$ID,]
save(myTab_long, file = paste0(UKB_phenotypes_filtered,"/01_Prep_02_UKB_HR_wideFormat.RData"))

#' # Transform wide to long ####
#' ***
#' According to the data dictionary from the BSU application, the first 112 trends should be correct (108-112 have always the same count: 18-16-8-6-1). At the end, I just want rows with complete data (HR, load, speed, phaseTime, phaseName, and stageName). Idealy this number of trends per person would then match the number of trends in my previous table. 
#' 
#' Here, I only want the columns ID, trend number, and then the values 
#' 
#' ## Heart Rate ####
Names_HR = names(myTab_long)[grepl("HR",names(myTab_long))]
myTab_HR <- melt(myTab_long,
                 id.vars=names(myTab_long)[c(1,2)],
                 measure.vars=Names_HR,
                 variable.name="trend")

myTab_HR[,trend := gsub("HR_trend","",trend)]
myTab_HR[,trend := as.numeric(trend)]
names(myTab_HR)[4] = "HR"

#' ## Load ####
Names_load = names(myTab_long)[grepl("Load_",names(myTab_long))]
myTab_load <- melt(myTab_long,
                 id.vars=names(myTab_long)[c(1,2)],
                 measure.vars=Names_load,
                 variable.name="trend")

myTab_load[,trend := gsub("Load_trend","",trend)]
myTab_load[,trend := as.numeric(trend)]
names(myTab_load)[4] = "Load"

#' ## Speed ####
Names_speed = names(myTab_long)[grepl("Speed_",names(myTab_long))]
myTab_speed <- melt(myTab_long,
                   id.vars=names(myTab_long)[c(1,2)],
                   measure.vars=Names_speed,
                   variable.name="trend")

myTab_speed[,trend := gsub("Speed_trend","",trend)]
myTab_speed[,trend := as.numeric(trend)]
names(myTab_speed)[4] = "Speed"

#' ## PhaseTime ####
Names_time = names(myTab_long)[grepl("PhaseTime_",names(myTab_long))]
myTab_time <- melt(myTab_long,
                   id.vars=names(myTab_long)[c(1,2)],
                   measure.vars=Names_time,
                   variable.name="trend")

myTab_time[,trend := gsub("PhaseTime_trend","",trend)]
myTab_time[,trend := as.numeric(trend)]
names(myTab_time)[4] = "PhaseTime"

#' ## PhaseName ####
Names_name1 = names(myTab_long)[grepl("PhaseName_",names(myTab_long))]
myTab_name1 <- melt(myTab_long,
                   id.vars=names(myTab_long)[c(1,2)],
                   measure.vars=Names_name1,
                   variable.name="trend")

myTab_name1[,trend := gsub("PhaseName_trend","",trend)]
myTab_name1[,trend := as.numeric(trend)]
names(myTab_name1)[4] = "PhaseName"

#' ## StageName ####
Names_name2 = names(myTab_long)[grepl("StageName_",names(myTab_long))]
myTab_name2 <- melt(myTab_long,
                   id.vars=names(myTab_long)[c(1,2)],
                   measure.vars=Names_name2,
                   variable.name="trend")

myTab_name2[,trend := gsub("StageName_trend","",trend)]
myTab_name2[,trend := as.numeric(trend)]
names(myTab_name2)[4] = "StageName"

#' ## Merge data ####
myTab_HR[,dumID := paste(ID,trend,sep=":")]
myTab_load[,dumID := paste(ID,trend,sep=":")]
myTab_speed[,dumID := paste(ID,trend,sep=":")]
myTab_time[,dumID := paste(ID,trend,sep=":")]
myTab_name1[,dumID := paste(ID,trend,sep=":")]
myTab_name2[,dumID := paste(ID,trend,sep=":")]

stopifnot(myTab_load$dumID == myTab_speed$dumID)
stopifnot(myTab_load$dumID == myTab_time$dumID)
stopifnot(myTab_load$dumID == myTab_name1$dumID)
stopifnot(myTab_load$dumID == myTab_name2$dumID)

matched = match(myTab_HR$dumID,myTab_load$dumID)
myTab_HR[,load := myTab_load[matched,Load]]
myTab_HR[,speed := myTab_speed[matched,Speed]]
myTab_HR[,phaseTime := myTab_time[matched,PhaseTime]]
myTab_HR[,phaseName := myTab_name1[matched,PhaseName]]
myTab_HR[,stageName := myTab_name2[matched,StageName]]

#' # Filter trend entries ####
#' ***
myTab_long = copy(myTab_HR)
save(myTab_long, file = paste0(UKB_phenotypes_filtered,"/01_Prep_02_UKB_HR_longFormat.RData"))

#' Filter NA entries
myTab_long = myTab_long[!is.na(HR),]
myTab_long[,table(is.na(load))]
myTab_long = myTab_long[!is.na(load),]
myTab_long[,table(is.na(speed))]
myTab_long[,table(is.na(phaseTime))]
myTab_long[,table(is.na(phaseName))]
myTab_long[,table(is.na(stageName))]

#' Extract HR before test (2nd of the 2 values) and last HR value after 1 min of rest, then remove those two phases
#' 
meanRestingHR = copy(myTab_long)[phaseName == "Pretest",]
setorder(meanRestingHR,-phaseTime)
meanRestingHR = meanRestingHR[!duplicated(ID),]
matched = match(myTab_cross$ID,meanRestingHR$ID)
myTab_cross[,HR_beforeTest := meanRestingHR[matched,HR]]

afterWorkoutHR = copy(myTab_long)[phaseName == "Rest",]
setorder(afterWorkoutHR,-phaseTime)
afterWorkoutHR = afterWorkoutHR[!duplicated(ID),]
matched = match(myTab_cross$ID,afterWorkoutHR$ID)
myTab_cross[,HR_afterTest := afterWorkoutHR[matched,HR]]

#' Exclude samples which have no before / after HR (something went wrong in the measurements)
myTab_cross = myTab_cross[!is.na(HR_beforeTest) & !is.na(HR_afterTest),]
myTab_long = myTab_long[ID %in% myTab_cross$ID,]
myTab_long = myTab_long[phaseName == "Exercise",]

ggplot(myTab_cross, aes(x=HR_beforeTest, y=HR_afterTest, col = as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="heart rate before exercise", y="heart rate after exercise")

ggplot(myTab_cross, aes(x=HR_beforeTest, y=MaxHR, col = as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="heart rate before exercise", y="maximum heart rate during exercise")

ggplot(myTab_cross, aes(x=HR_afterTest, y=MaxHR, col = as.factor(protocolID))) +
  geom_point()+
  facet_grid(. ~ Sex,
             labeller = as_labeller(c("0" = "Women", "1" = "Men"))) + 
  theme(legend.position = "none")+
  labs(x="heart rate after exercise", y="maximum heart rate during exercise")

#' Try and get some trajectories
plotData = copy(myTab_long)
plotData = plotData[ID %in% myTab_cross[protocolID %in% c(66,88),ID],]
plotData[,myShape := 1]
plotData[stageName=="Constant",myShape := myShape+1]
plotData[,myShape := as.factor(myShape)]

ggplot(plotData, aes(x=phaseTime, y=HR, group=ID,shape=myShape)) +
  geom_line(aes(alpha=0.01,col=as.factor(Sex))) + 
  geom_point(aes(colour = as.factor(Sex)))+
  labs(x="Time (in seconds)",y="heart rate (in bpm)", color="Sex",shape="Stage",
       title = paste0("UKB Field 5983 filtered for protocol ID 66 (women; 30 -> 90) and 88 (men, 40 -> 120)")) +
  scale_shape_manual(values=c(21,24),
                     labels = c("Ramp-up","Constant"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("women","men"))+
  theme_classic() 

ggplot(plotData, aes(x=phaseTime, y=load, group=ID,shape=myShape)) +
  geom_line(aes(alpha=0.01,col=as.factor(Sex))) + 
  geom_point(aes(colour = as.factor(Sex)))+
  labs(x="Time (in seconds)",y="load (in watt)", color="Sex",shape="Stage",
       title = paste0("UKB Field 5984 filtered for protocol ID 66 (women; 30 -> 90) and 88 (men, 40 -> 120)")) +
  scale_shape_manual(values=c(21,24),
                     labels = c("Ramp-up","Constant"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("women","men"))+
  theme_classic() 

ggplot(plotData, aes(x=phaseTime, y=speed, group=ID,shape=myShape)) +
  geom_line(aes(alpha=0.01,col=as.factor(Sex))) + 
  geom_point(aes(colour = as.factor(Sex)))+
  labs(x="Time (in seconds)",y="revolutions (per minute)", color="Sex",shape="Stage",
       title = paste0("UKB Field 5985 filtered for protocol ID 66 (women; 30 -> 90) and 88 (men, 40 -> 120)")) +
  scale_shape_manual(values=c(21,24),
                     labels = c("Ramp-up","Constant"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("women","men"))+
  theme_classic()

#' # Save data ####
#' ***
save(myTab_long, myTab_cross, file = paste0(UKB_phenotypes_filtered,"/01_Prep_02_UKB_HR.RData"))
write.table(myTab_cross$ID,file = paste0(UKB_phenotypes_filtered,"/01_Prep_02_SampleList_HR.txt"), 
            col.names = F, row.names = F, quote = F)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
