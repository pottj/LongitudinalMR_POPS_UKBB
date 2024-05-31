#' ---
#' title: "Get Trajectories"
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
#' I want the trajectories for EFW (in kg) including the final birth weight, colored by eCS. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2024-","24-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
myGD_files = list.files(path = POPS_phenotypes,pattern = "01_Prep_04")
myGD_files = myGD_files[grepl("240517",myGD_files)]
loaded1 = load(paste0(POPS_phenotypes,myGD_files[1]))
loaded1
loaded2 = load(paste0(POPS_phenotypes,myGD_files[2]))
loaded2

#' # Get plot data ####
#' ***
names(myTab_X)
myTab_X =myTab_X[,c(1,2,4,5,6,18,23,24)]
dummy = myTab_X[!is.na(efwcomb),.N,by=POPSID]
myTab_X = myTab_X[POPSID %in% dummy[N>=2,POPSID]]

names(myTab_Y)
myTab_Y[,scan := 4]
myTab_Y =myTab_Y[,c(1,9,18,55,10,12,16,15)]
names(myTab_Y) = names(myTab_X)
myTab_Y = myTab_Y[POPSID %in% myTab_X$POPSID]

plotData= rbind(myTab_X,myTab_Y)

plotData[,myShape := 1]
plotData[scan==4,myShape := myShape+1]
plotData[,myShape := as.factor(myShape)]

#' # Do plotting ####
#' ***
#' ## Absolute values ####
plotData[,sd(efwcomb/1000,na.rm=T),by=scan]

ggp1  = ggplot(plotData[scan!=4,], aes(x=ga, y=efwcomb/1000,  
                             group=POPSID,shape=myShape)) +
  #facet_wrap(~pn_emcsall,scales = "free_y")+
  geom_line(aes(alpha=0.01,col=as.factor(pn_emcsall))) + 
  geom_point(aes(colour = "black",fill = as.factor(pn_emcsall)))+
  labs(x="Gestational Week", y="Weight (in kg)", 
       fill="eCS?",shape="",color="eCS?") +
  scale_shape_manual(values=c(21,24),
                     labels = c("estimated \nfetal weight","birth weight"))+
  scale_fill_manual(values = c("steelblue","darkred"),
                    labels = c("no", "yes"))+
  scale_colour_manual(values = c("steelblue","darkred","black"),
                      labels = c("no", "yes","idk"))+
  theme(legend.position = "none") + theme_classic() +
  theme(legend.position = "none")

ggp1

#' I really dont know how to get rid of the black color.
#' 
png(file=paste0("../results/_figures/05_Trajectory_EFWR.png"),
    width=1500,height=1200,res = 200)
print(ggp1)
dev.off()

#' ## Log-transformed values
plotData[,sd(log(efwcomb/1000),na.rm=T),by=scan]

ggp2  = ggplot(plotData[scan!=4,], aes(x=ga, y=log(efwcomb/1000),  
                             group=POPSID,shape=myShape)) +
  #facet_wrap(~pn_emcsall,scales = "free_y")+
  geom_line(aes(alpha=0.01,col=as.factor(pn_emcsall))) + 
  geom_point(aes(colour = "black",fill = as.factor(pn_emcsall)))+
  labs(x="Gestational Week", y="Weight (log-transformed)", 
       fill="eCS?",shape="",color="eCS?") +
  scale_shape_manual(values=c(21,24),
                     labels = c("estimated \nfetal weight","birth weight"))+
  scale_fill_manual(values = c("steelblue","darkred"),
                    labels = c("no", "yes"))+
  scale_colour_manual(values = c("steelblue","darkred","black"),
                      labels = c("no", "yes","idk"))+
  theme(legend.position = "none") + theme_classic() +
  theme(legend.position = "none")
ggp2

#' I really dont know how to get rid of the black color.
#' 
png(file=paste0("../results/_figures/05_Trajectory_EFWL.png"),
    width=1500,height=1200,res = 200)
print(ggp2)
dev.off()

#' ## Z-scores ####
plotData[,sd(efwcombZv2,na.rm=T),by=scan]

ggp3  = ggplot(plotData[scan!=4,], aes(x=ga, y=efwcombZv2,  
                             group=POPSID,shape=myShape)) +
  #facet_wrap(~pn_emcsall,scales = "free_y")+
  geom_line(aes(alpha=0.01,col=as.factor(pn_emcsall))) + 
  geom_point(aes(colour = "black",fill = as.factor(pn_emcsall)))+
  labs(x="Gestational Week", y="Weight (GA-adjusted Z-scores)", 
       fill="eCS?",shape="",color="eCS?") +
  scale_shape_manual(values=c(21,24),
                     labels = c("estimated \nfetal weight","birth weight"))+
  scale_fill_manual(values = c("steelblue","darkred"),
                    labels = c("no", "yes"))+
  scale_colour_manual(values = c("steelblue","darkred","black"),
                      labels = c("no", "yes","idk"))+
  theme(legend.position = "none") + theme_classic() +
  theme(legend.position = "none")
ggp3

#' I really dont know how to get rid of the black color.
#' 
png(file=paste0("../results/_figures/05_Trajectory_EFWZ.png"),
    width=1500,height=1200,res = 200)
print(ggp3)
dev.off()

#' ## Centiles ####
plotData[,sd(efwcombv2_cent,na.rm=T),by=scan]

ggp4  = ggplot(plotData[scan!=4,], aes(x=ga, y=efwcombv2_cent,  
                             group=POPSID,shape=myShape)) +
  #facet_wrap(~pn_emcsall,scales = "free_y")+
  geom_line(aes(alpha=0.01,col=as.factor(pn_emcsall))) + 
  geom_point(aes(colour = "black",fill = as.factor(pn_emcsall)))+
  labs(x="Gestational Week", y="Weight (GA-adjusted Centiles)", 
       fill="eCS?",shape="",color="eCS?") +
  scale_shape_manual(values=c(21,24),
                     labels = c("estimated \nfetal weight","birth weight"))+
  scale_fill_manual(values = c("steelblue","darkred"),
                    labels = c("no", "yes"))+
  scale_colour_manual(values = c("steelblue","darkred","black"),
                      labels = c("no", "yes","idk"))+
  theme(legend.position = "none") + theme_classic() +
  theme(legend.position = "none")
ggp4

#' I really dont know how to get rid of the black color.
#' 
png(file=paste0("../results/_figures/05_Trajectory_EFWC.png"),
    width=1500,height=1200,res = 200)
print(ggp4)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
