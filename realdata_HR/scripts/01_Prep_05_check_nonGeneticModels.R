#' ---
#' title: "Check UKB non-genetic models"
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
#' Here, I want to estimate the none-genetic effects in a **gamlssIA** model for HR. My plan is to use the most simple model as main model (age, sex, and time), but include others in a sensitivity model (age, sex, time, load, speed, )
#' 
#' To be completed over time. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load and prep UKB data ####
#' ***
#' Load data (genetic data & phenotype data)
load(paste0(UKB_phenotypes_filtered,"/01_Prep_02_UKB_HR.RData"))

#' # Prepare data ####
#' ***
#' Check protocols again - I would like to have a sample size of roughly 10,000 (similar to simulation)
myTab_cross[,table(protocolID,targetLoad)]
myTab_cross[,factor := targetLoad/warmUp_load]
myTab_cross[,table(protocolID,factor>2 & factor<=3)]
myTab_cross[,table(factor>2 & factor<=3,sex_ECG)]
myTab_cross = myTab_cross[factor>2 & factor<=3,]

#' ## Get set for constant load
#' 
#' - factor >2 & factor<=3
#' - stage name == "constant"
#' - load > 30 (exclude trends at the very beginning)
#' - number of trends >8 & <=13
#' 
myTab_cross2 = copy(myTab_cross)
myTab_long2 = copy(myTab_long)
myTab_long2 = myTab_long2[ID %in% myTab_cross2$ID,]
myTab_long2 = myTab_long2[stageName == "Constant",]
myTab_long2 = myTab_long2[load>=30,]
test = myTab_long2[,.N,by=ID]
hist(test$N)
myTab_long2 = myTab_long2[ID %in% test[N>8 & N<14,ID]]
myTab_cross2 = myTab_cross2[ID %in% myTab_long2$ID,]

#' ## Get set for ramp-up
#' 
#' - factor >2 & factor<=3
#' - stage name != "constant"
#' - speed within 40-60 at all times
#' - minimal heart rate >=50

myTab_cross3 = copy(myTab_cross)
myTab_long3 = copy(myTab_long)
myTab_long3 = myTab_long3[ID %in% myTab_cross3$ID,]
myTab_long3 = myTab_long3[stageName != "Constant",]

test = myTab_long3[,max(speed),by=ID]
hist(test$V1)
myTab_long3 = myTab_long3[ID %in% test[V1<=80,ID]]
myTab_cross3 = myTab_cross3[ID %in% myTab_long3$ID,]

test = myTab_long3[,min(speed),by=ID]
hist(test$V1)
myTab_long3 = myTab_long3[ID %in% test[V1>=40,ID]]
myTab_cross3 = myTab_cross3[ID %in% myTab_long3$ID,]

test = myTab_long3[,min(HR),by=ID]
hist(test$V1)
myTab_long3 = myTab_long3[ID %in% test[V1>=50,ID]]
myTab_cross3 = myTab_cross3[ID %in% myTab_long3$ID,]

#' ## Get overlap 
#' 
#' I want only samples which are ok for both constant and ramp-up phase
#' 
myTab_long4 = copy(myTab_long)
myTab_long4 = myTab_long4[ID %in% myTab_long3$ID & ID %in% myTab_long2$ID,]
myTab_long4[,dumID2 := paste(dumID,stageName,sep=":")]
myTab_long3[,dumID2 := paste(dumID,stageName,sep=":")]
myTab_long2[,dumID2 := paste(dumID,stageName,sep=":")]
myTab_long4 = myTab_long4[dumID2 %in% myTab_long3$dumID2 | dumID2 %in% myTab_long2$dumID2,]
myTab_cross4 = copy(myTab_cross)
myTab_cross4 = myTab_cross4[ID %in% myTab_long4$ID]
myTab_cross4[,table(sex_ECG)]

ggplot(myTab_long4[stageName=="Constant"], aes(x=phaseTime, y=HR, group=ID)) +
  geom_line(aes(alpha=0.01,col=as.factor(Sex))) + 
  geom_point(aes(colour = as.factor(Sex)))+
  labs(x="Time (in seconds)",y="Heart rate (in bpm)", color="Sex",
       title = paste0("HR (stage = constant)")) +
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("women","men"))+
  theme_classic() 

ggplot(myTab_long4[stageName=="Constant"], aes(x=phaseTime, y=speed, group=ID)) +
  geom_line(aes(alpha=0.01,col=as.factor(Sex))) + 
  geom_point(aes(colour = as.factor(Sex)))+
  labs(x="Time (in seconds)",y="Revolutions (per minute)", color="Sex",
       title = paste0("Speed (stage = constant)")) +
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("women","men"))+
  theme_classic() 

ggplot(myTab_long4[stageName!="Constant"], aes(x=phaseTime, y=HR, group=ID)) +
  geom_line(aes(alpha=0.01,col=as.factor(Sex))) + 
  geom_point(aes(colour = as.factor(Sex)))+
  labs(x="Time (in seconds)",y="Heart rate (in bpm)", color="Sex",
       title = paste0("HR (stage = ramp-up)")) +
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("women","men"))+
  theme_classic() 

ggplot(myTab_long4[stageName!="Constant"], aes(x=phaseTime, y=speed, group=ID)) +
  geom_line(aes(alpha=0.01,col=as.factor(Sex))) + 
  geom_point(aes(colour = as.factor(Sex)))+
  labs(x="Time (in seconds)",y="Revolutions (per minute)", color="Sex",
       title = paste0("Speed (stage = ramp-up)")) +
  scale_colour_manual(values = c("darkred","steelblue"),
                      labels = c("women","men"))+
  theme_classic() 

save(myTab_long4, myTab_cross4, file = paste0(UKB_phenotypes_filtered,"/01_Prep_05_UKB_HR_filtered.RData"))

#' # GAMLSS non-genetics ####
#' ***
matched = match(myTab_long4$ID,myTab_cross4$ID)
myTab_long4 = cbind(myTab_long4,myTab_cross4[matched,c(25:35,39)])
myTab_long4[Sex==0,Sex:=2]

#' ## Constant stage
#' 
#' - Minimal model: just Age and Sex and phaseTime (RI + RS)
#' - Full model: Age, Sex, speed, protocolID, and phaseTime
#' - sex-combined, and sex-stratified
#' 
#' ### Sex combined ####
data0 = copy(myTab_long4)
data0 = data0[stageName == "Constant",]

time1<-Sys.time()
mod0_RS = gamlss(HR ~  Sex + Age + phaseTime + re(random=~1+phaseTime|as.factor(ID)),   
                  sigma.formula = ~ Sex + Age + phaseTime + re(random=~1|as.factor(ID)), 
                  data = data0, family = "NO")
message("\nTIME FOR CONSTANT LOAD IN SEX COMBINED SETTING & RI + RS: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

time1<-Sys.time()
mod0_RI = gamlss(HR ~  Sex + Age + phaseTime + re(random=~1|as.factor(ID)), 
                  sigma.formula = ~ Sex + Age + phaseTime, 
                  data = data0, family = "NO")
message("\nTIME FOR CONSTANT LOAD IN SEX COMBINED SETTING & RI only: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

summary(mod0_RS)
summary(mod0_RI)

pred_RS <- predict(mod0_RS, newdata=data0, type = "response")
pred_RI <- predict(mod0_RI, newdata=data0, type = "response")

plot(data0$HR,pred_RS)
plot(data0$HR,pred_RI)
plot(pred_RS,pred_RI)
plot(data0$phaseTime,pred_RS)
plot(data0$phaseTime,pred_RI)
plot(data0$phaseTime,data0$HR)

#' ### Men only ####
data1 = copy(data0)
data1 = data1[Sex == 1,]

time1<-Sys.time()
mod1_RS = gamlss(HR ~  Age + phaseTime + re(random=~1+phaseTime|as.factor(ID)),   
                  sigma.formula = ~ Age + phaseTime + re(random=~1|as.factor(ID)), 
                  data = data1, family = "NO")
message("\nTIME FOR CONSTANT LOAD IN MALE SETTING & RI + RS: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

time1<-Sys.time()
mod1_RI = gamlss(HR ~  Age + phaseTime + re(random=~1|as.factor(ID)),   
                  sigma.formula = ~ Age + phaseTime, 
                  data = data1, family = "NO")
message("\nTIME FOR CONSTANT LOAD IN MALE SETTING & RI only: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

summary(mod1_RS)
summary(mod1_RI)

pred_RS <- predict(mod1_RS, newdata=data1, type = "response")
pred_RI <- predict(mod1_RI, newdata=data1, type = "response")

plot(data1$HR,pred_RS)
plot(data1$HR,pred_RI)
plot(pred_RS,pred_RI)
plot(data1$phaseTime,pred_RS)
plot(data1$phaseTime,pred_RI)
plot(data1$phaseTime,data1$HR)

#' ### Women only ####
data2 = copy(data0)
data2 = data2[Sex == 2,]

time1<-Sys.time()
mod2_RS = gamlss(HR ~  Age + phaseTime + re(random=~1+phaseTime|as.factor(ID)),   
                  sigma.formula = ~ Age + phaseTime + re(random=~1|as.factor(ID)), 
                  data = data2, family = "NO")
message("\nTIME FOR CONSTANT LOAD IN FEMALE SETTING & RI + RS: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

time1<-Sys.time()
mod2_RI = gamlss(HR ~  Age + phaseTime + re(random=~1|as.factor(ID)), 
                  sigma.formula = ~ Age + phaseTime, 
                  data = data2, family = "NO")
message("\nTIME FOR CONSTANT LOAD IN FEMALE SETTING & RI only: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

summary(mod2_RS)
summary(mod2_RI)

pred_RS <- predict(mod2_RS, newdata=data2, type = "response")
pred_RI <- predict(mod2_RI, newdata=data2, type = "response")

plot(data2$HR,pred_RS)
plot(data2$HR,pred_RI)
plot(pred_RS,pred_RI)
plot(data2$phaseTime,pred_RS)
plot(data2$phaseTime,pred_RI)
plot(data2$phaseTime,data2$HR)

#' ## Ramp-up stages
#' 
#' - Minimal model: just Age and Sex and trend and load
#' - Full model: Age, Sex, speed, protocolID, and trend and load
#' - sex-combined, and sex-stratified
#' 
#' ### Sex combined ####
data0 = copy(myTab_long4)
data0 = data0[stageName != "Constant",]

time1<-Sys.time()
mod0_RS = gamlss(HR ~  Sex + Age + load + phaseTime + re(random=~1+phaseTime|as.factor(ID)),   
                  sigma.formula = ~ Sex + Age + load + phaseTime + re(random=~1|as.factor(ID)), 
                  data = data0, family = "NO")
message("\nTIME FOR RAMP-UP IN SEX COMBINED SETTING & RI + RS: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

time1<-Sys.time()
mod0_RI = gamlss(HR ~  Sex + Age + load + phaseTime + re(random=~1|as.factor(ID)), 
                  sigma.formula = ~ Sex + Age + load + phaseTime, 
                  data = data0, family = "NO")
message("\nTIME FOR RAMP-UP IN SEX COMBINED SETTING & RI only: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

summary(mod0_RS)
summary(mod0_RI)

pred_RS <- predict(mod0_RS, newdata=data0, type = "response")
pred_RI <- predict(mod0_RI, newdata=data0, type = "response")

plot(data0$HR,pred_RS)
plot(data0$HR,pred_RI)
plot(pred_RS,pred_RI)
plot(data0$phaseTime,pred_RS)
plot(data0$phaseTime,pred_RI)
plot(data0$phaseTime,data0$HR)

#' ### Men only ####
data1 = copy(data0)
data1 = data1[Sex == 1,]

time1<-Sys.time()
mod1_RS = gamlss(HR ~  Age + load + phaseTime + re(random=~1+phaseTime|as.factor(ID)),   
                  sigma.formula = ~ Age + load + phaseTime + re(random=~1|as.factor(ID)), 
                  data = data1, family = "NO")
message("\nTIME FOR RAMP-UP IN MALE SETTING & RI + RS: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

time1<-Sys.time()
mod1_RI = gamlss(HR ~  Age + load + phaseTime + re(random=~1|as.factor(ID)), 
                  sigma.formula = ~ Age + load + phaseTime, 
                  data = data1, family = "NO")
message("\nTIME FOR RAMP-UP IN MALE SETTING & RI only: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

summary(mod1_RS)
summary(mod1_RI)

pred_RS <- predict(mod1_RS, newdata=data1, type = "response")
pred_RI <- predict(mod1_RI, newdata=data1, type = "response")

plot(data1$HR,pred_RS)
plot(data1$HR,pred_RI)
plot(pred_RS,pred_RI)
plot(data1$phaseTime,pred_RS)
plot(data1$phaseTime,pred_RI)
plot(data1$phaseTime,data1$HR)

#' ### Women only ####
data2 = copy(data0)
data2 = data2[Sex == 2,]

time1<-Sys.time()
mod2_RS = gamlss(HR ~  Age + load + phaseTime + re(random=~1+phaseTime|as.factor(ID)),   
                  sigma.formula = ~ Age + load + phaseTime + re(random=~1|as.factor(ID)), 
                  data = data2, family = "NO")
message("\nTIME FOR RAMP-UP IN FEMALE SETTING & RI + RS: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

time1<-Sys.time()
mod2_RI = gamlss(HR ~  Age + load + phaseTime + re(random=~1|as.factor(ID)), 
                  sigma.formula = ~ Age + load + phaseTime, 
                  data = data2, family = "NO")
message("\nTIME FOR RAMP-UP IN FEMALE SETTING & RI only: \n" ,round(difftime(Sys.time(),time1,units = "mins"),3)," minutes")

summary(mod2_RS)
summary(mod2_RI)

pred_RS <- predict(mod2_RS, newdata=data2, type = "response")
pred_RI <- predict(mod2_RI, newdata=data2, type = "response")

plot(data2$HR,pred_RS)
plot(data2$HR,pred_RI)
plot(pred_RS,pred_RI)
plot(data2$phaseTime,pred_RS)
plot(data2$phaseTime,pred_RI)
plot(data2$phaseTime,data2$HR)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
