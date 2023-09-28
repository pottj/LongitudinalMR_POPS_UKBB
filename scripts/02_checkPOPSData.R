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
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_laptop.R")
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


#' # Reorganize data ####
#' ***
#' I want a long format, with each row being the data from 1 individual at 1 specific time point. 
#' 
#' - time 1: pre-pregnancy (mothers weight only)
#' - time 2: scan 2
#' - time 3: scan 3
#' - time 4: scan 4
#' - time 5: post-birth (birth weight)
#' 
#' To make the melting easier, I rename the columns for pre-pregnancy weight, and birth weight
#' 
setnames(data,"an_weight","an_scan1_weight")
setnames(data,"pn_ga_wk","an_scan5_ga")
setnames(data,"pn_bw","an_scan5_efw")
setnames(data,"pn_bwpct","an_scan5_efw_centileHadlock")

myVariables = names(data)[grepl("an_scan2",names(data))]

dumTab = foreach(i = 1:length(myVariables))%do%{
  #i=1
  myVar = myVariables[i]
  myVars = c(myVar,gsub("2","3",myVar),gsub("2","4",myVar))
  if(i==1) myVars = c("an_scan1_weight",myVars)
  if(i==10) myVars = c(myVars,"an_scan5_ga")
  if(i==11) myVars = c(myVars,"an_scan5_efw")
  if(i==15) myVars = c(myVars,"an_scan5_efw_centileHadlock")
  dummy = myVars
  dummy
}
names(dumTab) = gsub("an_scan2_","",myVariables)

data2 <- melt(data,
              # ID variables - all the variables to keep but not split apart on
              id.vars=names(data)[c(1:3,5,87)],
              # The source columns
              measure.vars=dumTab[c(1:11,13:15,22:23,26)],
              # Name of the destination column that will identify the original
              # column that the measurement came from
              variable.name="scan"
)

#' # Plotting ####
#' ***
#' 
#' ## EFW ####
#' 
ggp<-ggplot(data2, aes(ga, efw) ) +
  geom_segment(aes(x = 20.00000, y = 405.1557, xend = 22.57143, yend = 1192.0366, colour = "red"))+
  geom_segment(aes(x = 22.57143, y = 1192.0366, xend = 36.14286, yend = 3828.4207, colour = "red"))+
  geom_segment(aes(x = 36.14286, y = 3828.4207, xend = 41.29000, yend = 4580.0000, colour = "red"))+
  geom_point()+
  geom_point(data=data2[POPSID==211], 
             aes(x=ga,y=efw), 
             color='red',
             size=3)+
  scale_colour_manual(values=c("red"),label="211")+
  labs(x="Gestational Week",y="Estimated fetal weight", color="POPSID") +
  geom_smooth(method = "loess",
              formula = y ~ x)
ggp

plot1 = ggplot()+
  geom_line(aes(y=efw, 
                x=ga, 
                group=POPSID, 
                colour=pn_sex), 
            data=data2, 
            show.legend = TRUE) + 
  labs(color="Baby sex") +
  theme(legend.position = "none") + theme_classic()
plot1

#' ## Femur length ####
#' 
ggp<-ggplot(data2, aes(ga, fl) ) +
  geom_segment(aes(x = 20.00000, y = 35, xend = 22.57143, yend = 51, colour = "red"))+
  geom_segment(aes(x = 22.57143, y = 51, xend = 36.14286, yend = 73, colour = "red"))+
  geom_point()+
  geom_point(data=data2[POPSID==211], 
             aes(x=ga,y=fl), 
             color='red',
             size=3)+
  scale_colour_manual(values=c("red"),label="211")+
  labs(x="Gestational Week",y="Femur length", color="POPSID") +
  geom_smooth(method = "loess",
              formula = y ~ x)
ggp

plot1 = ggplot()+
  geom_line(aes(y=fl, 
                x=ga, 
                group=POPSID, 
                colour=pn_sex), 
            data=data2, 
            show.legend = TRUE) + 
  labs(color="Baby sex") +
  theme(legend.position = "none") + theme_classic()
plot1

#' ## Head circumference ####
#' 
ggp<-ggplot(data2, aes(ga, hc) ) +
  geom_segment(aes(x = 20.00000, y = 185, xend = 22.57143, yend = 267, colour = "red"))+
  geom_segment(aes(x = 22.57143, y = 267, xend = 36.14286, yend = 344, colour = "red"))+
  geom_point()+
  geom_point(data=data2[POPSID==211], 
             aes(x=ga,y=hc), 
             color='red',
             size=3)+
  scale_colour_manual(values=c("red"),label="211")+
  labs(x="Gestational Week",y="Head circumference", color="POPSID") +
  geom_smooth(method = "loess",
              formula = y ~ x)
ggp

plot1 = ggplot()+
  geom_line(aes(y=hc, 
                x=ga, 
                group=POPSID, 
                colour=pn_sex), 
            data=data2, 
            show.legend = TRUE) + 
  labs(color="Baby sex") +
  theme(legend.position = "none") + theme_classic()
plot1

#' ## Abdominal circumference ####
#' 
ggp<-ggplot(data2, aes(ga, ac) ) +
  geom_segment(aes(x = 20.00000, y = 163, xend = 22.57143, yend = 244, colour = "red"))+
  geom_segment(aes(x = 22.57143, y = 244, xend = 36.14286, yend = 364, colour = "red"))+
  geom_point()+
  geom_point(data=data2[POPSID==211], 
             aes(x=ga,y=ac), 
             color='red',
             size=3)+
  scale_colour_manual(values=c("red"),label="211")+
  labs(x="Gestational Week",y="Abdominal circumference", color="POPSID") +
  geom_smooth(method = "loess",
              formula = y ~ x)
ggp

plot1 = ggplot()+
  geom_line(aes(y=ac, 
                x=ga, 
                group=POPSID, 
                colour=pn_sex), 
            data=data2, 
            show.legend = TRUE) + 
  labs(color="Baby sex") +
  theme(legend.position = "none") + theme_classic()
plot1

#' ## Biparietal diameter ####
#' 
ggp<-ggplot(data2, aes(ga, bpd) ) +
  geom_segment(aes(x = 20.00000, y = 50, xend = 22.57143, yend = 70, colour = "red"))+
  geom_segment(aes(x = 22.57143, y = 70, xend = 36.14286, yend = 97, colour = "red"))+
  geom_point()+
  geom_point(data=data2[POPSID==211], 
             aes(x=ga,y=bpd), 
             color='red',
             size=3)+
  scale_colour_manual(values=c("red"),label="211")+
  labs(x="Gestational Week",y="Biparietal diameter", color="POPSID") +
  geom_smooth(method = "loess",
              formula = y ~ x)
ggp

plot1 = ggplot()+
  geom_line(aes(y=bpd, 
                x=ga, 
                group=POPSID, 
                colour=pn_sex), 
            data=data2, 
            show.legend = TRUE) + 
  labs(color="Baby sex") +
  theme(legend.position = "none") + theme_classic()
plot1

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
