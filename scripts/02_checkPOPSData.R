#' ---
#' title: "Check POPS data"
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
#' I receive the POPS data from Ulla (see email from 26/09/2023). This individual level data cannot be shared via github or any other repository. Hence, I stored it under /data/IndividualLevelData/, which is ignored by git.
#' 
#' Here, I want to load the data and check out the time-(in)dependent trajectories for the growth parameters: 
#' 
#' - mothers weight
#' - estimated fetal weight (EFW)
#' - head circumference (HC)
#' - femur length (FL)
#' - biparietal diameter (BPD)
#' - abdominal circumference (AC)
#' 
#' 
#' # Initialize ####
#' ***
#' 
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2023-","23-",tag)
tag = gsub("-","",tag)

#' # Load POS data ####
#' ***
#' 
data = data.table(read_xlsx("../data/IndividualLevelData/POPS_data_extract_JP_26sep2023.xlsx"))
names(data)
names(data)[grepl("an_scan2_",names(data))]

#' The gestational date is not seperated by "_", so I will change that column name

names(data)[c(15,40,65)]
names(data)[c(15,40,65)] = c("an_scan2_ga","an_scan3_ga","an_scan4_ga")


#' # Re-organize data ####
#' ***
#' The data was already filtered by Ulla Sovio for. Exclusion criteria were: 
#' 
#' - No biometry available at any time point,
#' - No information on the mode of delivery,
#' - Preterm birth,
#' - Non-cephalic presentation at delivery,
#' - Prelabor CD,
#' - Antepartum stillbirth, or
#' - Preexisting diabetes.
#' 
#' Gestational diabetes might affect fetal growth and will be used as covariable only. 
#' 
#' Right now, I have all data in one table. I want to split the data a bit: 
#' 
#' - longitudinal exposure data: data per scan (26 parameters)
#' - outcome data: data at birth (sex, GA, birth weight, caesarean section)
#' - covariable data: more information about mother and partner (height, weight, presentation at birth, gestational diabetes)
#' 
#' Big goal: Understand the data I got (each parameter!)
#' 
#' ## Longitudinal exposure data ####
#' ***
myVariables = names(data)[grepl("an_scan2",names(data))]

dumTab = foreach(i = 1:length(myVariables))%do%{
  #i=12
  myVar = myVariables[i]
  myVars = c(myVar,gsub("scan2","scan3",myVar),gsub("scan2","scan4",myVar))
  dummy = myVars
  dummy
}
names(dumTab) = gsub("an_scan2_","",myVariables)

myTab_X <- melt(data,
              id.vars=names(data)[c(1,87,90,96)],
              measure.vars=dumTab,
              variable.name="scan")
myTab_X = myTab_X[,c(1:5,15,6:14,16:31)]

for(i in 1:(length(myVariables)-1)){
  #i=2
  message("Working on variable: ",names(myTab_X)[i+6])
  myVar = names(myTab_X)[i+6]
  myTab_X[,value := get(myVar)]
  
  ggp1  = ggplot(myTab_X, aes(x=ga, y=value, col=pn_sex) ) +
    geom_point()+
    labs(x="Gestational Week",y=myVar, color="Babys sex") +
    theme(legend.position = "none") + theme_classic() +
    geom_smooth(method = "loess",
                formula = y ~ x)
  print(ggp1)
  
  ggp2 = ggplot(myTab_X, aes(x=ga, y=value, col=pn_sex, group=POPSID))+
    geom_line(show.legend = TRUE) + 
    labs(x="Gestational Week",y=myVar, color="Babys sex") +
    theme(legend.position = "none") + theme_classic()
  print(ggp2)
  
  ggp3  = ggplot(myTab_X, aes(x=ga, y=value, col=as.factor(pn_emcsall)) ) +
    geom_point()+
    labs(x="Gestational Week",y=myVar, color="Emergency CS") +
    theme(legend.position = "none") + theme_classic() +
    geom_smooth(method = "loess",
                formula = y ~ x)
  print(ggp3)
  
  ggp4 = ggplot(myTab_X, aes(x=ga, y=value, col=as.factor(pn_emcsall), group=POPSID))+
    geom_line(show.legend = TRUE) + 
    labs(x="Gestational Week",y=myVar, color="Emergency CS") +
    theme(legend.position = "none") + theme_classic()
  print(ggp4)
  
  myTab_X[!is.na(value),interaction.plot(x.factor = scan,
                                         trace.factor = pn_sex,
                                         response = value,
                                         xlab = "#Scan",
                                         ylab = paste0("Mean of ",myVar),
                                         col = c("blue","red"),
                                         trace.label = "Babys sex")]
  myTab_X[!is.na(value),interaction.plot(x.factor = scan,
                                         trace.factor = as.factor(pn_emcsall),
                                         response = value,
                                         xlab = "#Scan",
                                         ylab = paste0("Mean of ",myVar),
                                         col = c("darkgreen","purple"),
                                         trace.label = "Emergency CS")]
  
}

save(myTab_X,file = "../data/IndividualLevelData/02_LongitudinalExposure.RData")

#' ## Outcome data ####
#' ***
myTab_Y = copy(data)
myTab_Y = myTab_Y[,!grepl("an_scan",names(data)),with=F]

myTab_Y = myTab_Y[,c(1:3,20:22,4:19,23:24)]

#' Check if there are any pre-term CS
myTab_Y[,table(pn_elcs)]

#' Check if there are non-cephalic presentations at birth 
myTab_Y[,table(pt_presdel)]

#' Okay, a lot of entries are NA. Those without NA look ok. 
#' 
#' Check gestational diabetes 
myTab_Y[,table(pn_gdm_diet_medication,is.na(pn_gdm_ga))]

#' Okay, there are 2 women with GDM with no GA of diagnosis. Maybe not that critical
#' 
myTab_Y[,boxplot(pn_gdm_ga)]

#' Check height of parents
ggp5  = ggplot(myTab_Y, aes(x=an_height, y=an_height_partner, col=pn_sex) ) +
  geom_point()+
  labs(x="Height Mother",y="Height Partner", color="Babys sex") +
  theme(legend.position = "none") + theme_classic() +
  geom_smooth(method = "loess",
              formula = y ~ x)
print(ggp5)

#' Check height of mother vs birth weight
ggp6  = ggplot(myTab_Y, aes(x=an_height, y=pn_bw, col=pn_sex) ) +
  geom_point()+
  labs(x="Height Mother",y="Birth Weight", color="Babys sex") +
  theme(legend.position = "none") + theme_classic() +
  geom_smooth(method = "loess",
              formula = y ~ x)
print(ggp6)

ggp7  = ggplot(myTab_Y, aes(x=an_height, y=pn_bw, col=as.factor(pn_emcsall)) ) +
  geom_point()+
  labs(x="Height Mother",y="Birth Weight", color="Emergency CS") +
  theme(legend.position = "none") + theme_classic() +
  geom_smooth(method = "loess",
              formula = y ~ x)
print(ggp7)

ggp8  = ggplot(myTab_Y, aes(x=an_height_partner, y=pn_bw, col=as.factor(pn_emcsall)) ) +
  geom_point()+
  labs(x="Height Partner",y="Birth Weight", color="Emergency CS") +
  theme(legend.position = "none") + theme_classic() +
  geom_smooth(method = "loess",
              formula = y ~ x)
print(ggp8)

save(myTab_Y,file = "../data/IndividualLevelData/02_Outcome.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
