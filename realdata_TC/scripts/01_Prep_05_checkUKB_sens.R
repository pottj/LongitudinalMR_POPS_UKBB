#' ---
#' title: "Get nice subset of TC data"
#' subtitle: "Longitudinal MVMR"
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
#' Here, I want to filter the data to a more structured data set. 
#' 
#' - no statin treatment (just to simplefy my life)
#' - ceiling of age --> per year just one TC measurement (first each)
#' - only age of 50-70 (statin treatment in older people could be higher ...)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load UKB data ####
#' ***
load(paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_GP_TC.RData"))

#' filter for no statin treatment
myTab6 = myTab6[lipLowMed==0,]

#' ceiling of age and filter for 1 value per year
myTab6[,age2:= ceiling(exposure_age)]
myTab6[,dumID2 := paste0(BSU_ID,"_",age2)]
myTab6 = myTab6[!duplicated(dumID2),]

#' filter for age group
myTab6[,min(age2)]
myTab6[,max(age2)]
myTab6 = myTab6[age2>=50 & age2<=70,]

#' check number of samples and number of time points
test = myTab6[,.N,by=BSU_ID]
hist(test$N)
test[,table(N>2)]

myTab6 = myTab6[BSU_ID %in% test[N>2,BSU_ID]]
myTab7 = myTab7[ID %in% test[N>2,BSU_ID]]

#' Try and get some trajectories
matched = match(myTab6$BSU_ID,myTab7$ID)
myTab_long = cbind(myTab6,myTab7[matched,c(2,7:16)])
myTab_long[sex==0,sex:=2]
myTab_long[,exposure_type := NULL]

plotData = copy(myTab_long)
plotData[,myShape := as.factor(lipLowMed)]

# test plotting with subset of samples
ggplot1 = ggplot(plotData, aes(x=age2, y=exposure_value, group=BSU_ID)) +
  geom_line(aes(col=as.factor(sex))) + 
  geom_point(aes(colour = as.factor(sex)))+
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  labs(x="Age (in years, ceiled)",
       y="Total cholesterol (in mmol/l)", 
       color="Sex",
       title = paste0("Trajectories of TC levels (sensitivity subset)")) +
  scale_shape_manual(values=c(21,24),
                     labels = c("no","yes"))+
  scale_colour_manual(values = c("steelblue","darkred"),
                      labels = c("men","women"))+
  theme_classic() 
plot(ggplot1)

ggplot2 = ggplot(plotData, aes(x=exposure_age, y=exposure_value, group=BSU_ID)) +
  geom_line(aes(col=as.factor(sex))) + 
  geom_point(aes(colour = as.factor(sex)))+
  geom_hline(yintercept = 5.17,linetype="dashed", linewidth=0.5)+
  geom_hline(yintercept = 6.18,linetype="dotted", linewidth=0.5)+
  labs(x="Age (in years, not rounded)",
       y="Total cholesterol (in mmol/l)", 
       color="Sex",
       title = paste0("Trajectories of TC levels (sensitivity subset)")) +
  scale_shape_manual(values=c(21,24),
                     labels = c("no","yes"))+
  scale_colour_manual(values = c("steelblue","darkred"),
                      labels = c("men","women"))+
  theme_classic() 
plot(ggplot2)

filename = paste0("../results/_figures/01_Trajectory_TC_sens_scan.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(ggplot1)
dev.off()

filename = paste0("../results/_figures/01_Trajectory_TC_sens_age.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(ggplot2)
dev.off()

#' # Save data ####
#' ***
save(myTab6, myTab7, file = paste0(UKB_phenotypes_filtered,"/01_Prep_05_UKB_GP_TC_sens.RData"))

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

