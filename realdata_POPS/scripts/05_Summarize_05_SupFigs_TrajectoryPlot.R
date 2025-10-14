#' ---
#' title: "Get Supplemental Figures for real data (trajectory plot)"
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
#' **Supplemental Figures: Trajectory Plot**
#' 
#' - MAIN (all samples): colored by sex
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile.R")

#' # Load data ####
#' ***
myGD_files = list.files(path = POPS_phenotypes,pattern = "01_Prep_04")
loaded1 = load(paste0(POPS_phenotypes,myGD_files[1]))
loaded1

myTab_X[,logefwcomb := log(efwcomb)]
myTab_X = myTab_X[!is.na(logefwcomb)]
test1 = myTab_X[,.N,by=POPSID]
myTab_X = myTab_X[POPSID %in% test1[N>1,POPSID]]

#' # Get plot data ####
#' ***
plotData = copy(myTab_X)
names(plotData)

#' # Get Plot ####
#' ***
ggp2  = ggplot(plotData, aes(x=ga, y=logefwcomb, group=POPSID,
                             col=pn_sex,fill=pn_sex)) +
  geom_line(aes(alpha=0.01)) + 
  geom_point(shape = 21,color="black")+
  labs(x="Gestational Age (weeks)", y="Estimated Foetal Weight (log-transformed)") +
  scale_fill_manual(values = c("darkred","steelblue"),
                      labels = c("FEMALE", "MALE"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("FEMALE", "MALE"))+
  theme_classic() +
  theme(legend.position = "none")
ggp2

filename = paste0("../results/_figures/SupFigs/Trajectory_EFWL_MAIN.png")
png(filename = filename,width = 1800, height = 1000, res=200)
plot(ggp2)
dev.off()

#' # Get plot data ####
#' ***
myTab_X = myTab_X[POPSID %in% test1[N==3,POPSID]]
plotData = copy(myTab_X)

#' # Get Plot ####
#' ***
ggp2  = ggplot(plotData, aes(x=ga, y=logefwcomb, group=POPSID,
                             col=pn_sex,fill=pn_sex)) +
  geom_line(aes(alpha=0.01)) + 
  geom_point(shape = 21,color="black")+
  labs(x="Gestational Age (weeks)", y="Estimated Foetal Weight (log-transformed)") +
  scale_fill_manual(values = c("darkred","steelblue"),
                    labels = c("FEMALE", "MALE"))+
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("FEMALE", "MALE"))+
  theme_classic() +
  theme(legend.position = "none")
ggp2

filename = paste0("../results/_figures/SupFigs/Trajectory_EFWL_SENS_sampleSet.png")
png(filename = filename,width=1500,height=1200,res = 200)
plot(ggp2)
dev.off()

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME of script (in minutes): " ,round(difftime(Sys.time(), time0,units = "mins"),2))


