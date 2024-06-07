#' ---
#' title: "Evaluation POPS simulation"
#' subtitle: "Summary of all simulations"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: html_document
#' ---
#'
#' # Introduction ####
#' 
#' # Initialize ####
#' ***

rm(list = ls())
time0<-Sys.time()

source("../../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Get to-do list ####
#' ***
#' get the file names (same names for all results directories)
mySims = list.files(path = "../results/",pattern = "_SimAll.RData")

ToDoFile = data.table(NR = 1:length(mySims),
                      file = mySims)
dummy = unlist(strsplit(mySims,"_"))
ToDoFile[,SimNR := dummy[seq(2,length(dummy),5)]]
ToDoFile[,SimX := dummy[seq(3,length(dummy),5)]]
ToDoFile[,SimY := dummy[seq(4,length(dummy),5)]]
ToDoFile[SimY == "CM1",theta1 := 0.3]
ToDoFile[SimY != "CM1",theta1 := 1.2]
ToDoFile[SimY == "CM4",theta2 := 0]
ToDoFile[SimY != "CM4",theta2 := 0.3]
ToDoFile[SimY %in% c("CM1","CM2"),theta3 := 0]
ToDoFile[SimY %in% c("CM3","CM4"),theta3 := 1]
ToDoFile

#' # Load files ####
#' ***
#' 
dumTab = foreach(i = 1:dim(ToDoFile)[1])%do%{
  #i=1
  myRow = ToDoFile[i,]
  message("Working on ",myRow$file)
  load(paste0("../results/",myRow$file))
  
  SimTab[exposure %in% c("mean"),exposure := "X1"]
  SimTab[exposure %in% c("slope"),exposure := "X2"]
  SimTab[exposure %in% c("var"),exposure := "X3"]

  # correct estimates for exposure age
  # SimTab[exposure == "X2",beta_IVW2 := beta_IVW/60]
  # SimTab[exposure == "X2",SE_IVW2 := SE_IVW/60]
  
  # add true value
  SimTab[,theta1 := 0]
  SimTab[outcome %in% c("Y2","Y5","Y6","Y8"),theta1 := myRow$theta1]
  SimTab[,theta2 := 0]
  SimTab[outcome %in% c("Y3","Y5","Y7","Y8"),theta2 := myRow$theta2]
  SimTab[,theta3 := 0]
  SimTab[outcome %in% c("Y4","Y6","Y7","Y8"),theta3 := myRow$theta3]
  
  SimTab[exposure == "X1", theta := theta1]
  SimTab[exposure == "X2", theta := theta2]
  SimTab[exposure == "X3", theta := theta3]
  SimTab[,theta1 := NULL]
  SimTab[,theta2 := NULL]
  SimTab[,theta3 := NULL]
  
  # get significant simulations
  SimTab[,dumID := paste(exposure,outcome,sep="_")]
  tab1 = SimTab[,.N,by = dumID]
  tab2 = SimTab[pval_IVW<0.05,.N,by = dumID]
  matched = match(tab1$dumID,tab2$dumID)
  tab1[,N_sig := tab2[matched,N]]
  tab1[is.na(N_sig),N_sig := 0]
  tab1[,N_sig_proc := N_sig/N ]
    
  # get bias 
  SimTab[,bias := beta_IVW2 - theta]
  tab3 = SimTab[,mean(bias),by = dumID]
  tab4 = SimTab[,sd(beta_IVW2),by = dumID]
  tab6 = SimTab[,mean(beta_IVW2),by = dumID]
  matched = match(tab1$dumID,tab3$dumID)
  tab1[,bias := tab3[matched,V1]]
  matched = match(tab1$dumID,tab4$dumID)
  tab1[,bias_SE := tab4[matched,V1]]
  tab1[,bias_SE := bias_SE * sqrt(1/N)]
    
  tab1[,mean_betaIVW2 := tab6[matched,V1]]
  tab1[,sd_betaIVW2 := tab4[matched,V1]]
  
  tab7 = SimTab[,sd(beta_IVW),by = dumID]
  tab8 = SimTab[,mean(beta_IVW),by = dumID]
  matched = match(tab1$dumID,tab8$dumID)
  tab1[,mean_betaIVW := tab8[matched,V1]]
  matched = match(tab1$dumID,tab7$dumID)
  tab1[,sd_betaIVW := tab7[matched,V1]]
    
  # get median numbers of SNPs
  stopifnot(sum(is.na(SimTab$NR_SNPs_type))==0)
  tab5 = SimTab[,summary(NR_SNPs_type),by=dumID]
  tab5[,V2 := rep(c("min","1stQ","median","mean","3rdQ","max"),length(unique(dumID)))]
  tab5_median = copy(tab5)[V2=="median"]
  tab5_1stQ = copy(tab5)[V2=="1stQ"]
  tab5_3rdQ = copy(tab5)[V2=="3rdQ"]
  
  matched = match(tab1$dumID,tab5_median$dumID)
  tab1[,NR_SNPs_median := as.numeric(tab5_median[matched,V1])]
  tab1[,NR_SNPs_1stQ := as.numeric(tab5_1stQ[matched,V1])]
  tab1[,NR_SNPs_3rdQ := as.numeric(tab5_3rdQ[matched,V1])]
    
  # get median cond F-Stat
  #stopifnot(sum(is.na(SimTab$condFStats))==0)
  tab5 = SimTab[!is.na(condFStats),summary(condFStats),by=dumID]
  tab6 = SimTab[,sum(is.na(condFStats)),by=dumID]
  
  tab5[,V2 := rep(c("min","1stQ","median","mean","3rdQ","max"),length(unique(dumID)))]
  tab5_median = copy(tab5)[V2=="median"]
  tab5_1stQ = copy(tab5)[V2=="1stQ"]
  tab5_3rdQ = copy(tab5)[V2=="3rdQ"]
  
  matched = match(tab1$dumID,tab5_median$dumID)
  tab1[,condFStats_median := as.numeric(tab5_median[matched,V1])]
  tab1[,condFStats_1stQ := as.numeric(tab5_1stQ[matched,V1])]
  tab1[,condFStats_3rdQ := as.numeric(tab5_3rdQ[matched,V1])]
  tab1[,condFstats_NA := as.numeric(tab6$V1)]
  
  # add result to myRow
  dummy = unlist(strsplit(tab1$dumID,"_"))
  tab1[,exposure := dummy[seq(1,length(dummy),2)]]
  tab1[,outcome := dummy[seq(2,length(dummy),2)]]
    
  tab1[,Sim_NR := myRow$SimNR]
  tab1[,Sim_X := myRow$SimX]
  tab1[,Sim_Y := myRow$SimY]
    
  tab1_wide = dcast(tab1, Sim_NR + Sim_X + Sim_Y + outcome ~ exposure, 
                      value.var = names(tab1)[2:17])
  tab1_wide
  
}

myTab = rbindlist(dumTab, fill = T)
names(myTab)
myTab = myTab[,c(1:4,
                 7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,
                 8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)]

#' # Check 1: Detection rates ####
#' ***

dumTab = copy(myTab)
dumTab[,dumID1 := paste(Sim_X,Sim_Y,sep="_")]

dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X1")
dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X2")

dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1])
x = dim(dumTab2)[2]
dumMat = as.matrix(dumTab2[,-1])
colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X2"),each=12),sep=" - ")

filt1 = grepl("X12_",colnames(dumMat))
corrplot(dumMat[,filt1], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filt2 = grepl("X123_",colnames(dumMat))
corrplot(dumMat[,filt2], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filt3 = grepl("X13_",colnames(dumMat))
corrplot(dumMat[,filt3], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)

filename = paste0("../results/_figures/DetectionRates_X12.png")
png(filename = filename,width = 2600, height = 1400, res=200)
corrplot(dumMat[,filt1], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

filename = paste0("../results/_figures/DetectionRates_X123.png")
png(filename = filename,width = 2600, height = 1400, res=200)
corrplot(dumMat[,filt2], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

filename = paste0("../results/_figures/DetectionRates_X133.png")
png(filename = filename,width = 2600, height = 1400, res=200)
corrplot(dumMat[,filt3], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()


#' # Check 2: Bias ####
#' ***
#' I want to plot the bias - again per outcome but only for continuous ones (does not really make sense for binary outcomes)
#' 
dumTab = copy(myTab)
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,bias_X1 := bias_X2]
dumTab2[,bias_SE_X1 := bias_SE_X2]
dumTab2[,type := "slope"]

dumTab4 = rbind(dumTab,dumTab2)
dumTab4 = dumTab4[!is.na(bias_X1),]

plot5 = ggplot(dumTab4[Sim_X=="X12"], aes(x=Sim_Y, y=bias_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=bias_X1-1.96*bias_SE_X1, ymax=bias_X1+1.96*bias_SE_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Bias") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/Bias_X12.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

plot5 = ggplot(dumTab4[Sim_X=="X123"], aes(x=Sim_Y, y=bias_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=bias_X1-1.96*bias_SE_X1, ymax=bias_X1+1.96*bias_SE_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Bias") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/Bias_X123.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

plot5 = ggplot(dumTab4[Sim_X=="X13"], aes(x=Sim_Y, y=bias_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=bias_X1-1.96*bias_SE_X1, ymax=bias_X1+1.96*bias_SE_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Bias") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/Bias_X13.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

#' # Check 3: Estimate ####
#' ***
#' I want to plot the corrected causal effect estimates - again per outcome!
#' 
dumTab = copy(myTab)
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,mean_betaIVW2_X1 := mean_betaIVW2_X2]
dumTab2[,sd_betaIVW2_X1 := sd_betaIVW2_X2]
dumTab2[,type := "slope"]

dumTab4 = rbind(dumTab,dumTab2)
dumTab4 = dumTab4[!is.na(mean_betaIVW2_X1),]

data_hlines = data.frame(type = c(rep("mean",4),"slope"),
                         mylines = c(0.3,1.2,-1.2,-0.3,0.3))

plot5 = ggplot(dumTab4[Sim_X == "X12",], aes(x=Sim_Y, y=mean_betaIVW2_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW2_X1-1.96*sd_betaIVW2_X1, ymax=mean_betaIVW2_X1+1.96*sd_betaIVW2_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Corrected estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/CorrectedEstimates_X12.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

plot5 = ggplot(dumTab4[Sim_X == "X123",], aes(x=Sim_Y, y=mean_betaIVW2_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW2_X1-1.96*sd_betaIVW2_X1, ymax=mean_betaIVW2_X1+1.96*sd_betaIVW2_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Corrected estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/CorrectedEstimates_X123.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

plot5 = ggplot(dumTab4[Sim_X == "X13",], aes(x=Sim_Y, y=mean_betaIVW2_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW2_X1-1.96*sd_betaIVW2_X1, ymax=mean_betaIVW2_X1+1.96*sd_betaIVW2_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Corrected estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/CorrectedEstimates_X13.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

#' # Check 4: raw estimate ####
#' ***
#' I want to plot the uncorrected causal effect estimates - again per outcome!
dumTab = copy(myTab)
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,mean_betaIVW_X1 := mean_betaIVW_X2]
dumTab2[,sd_betaIVW_X1 := sd_betaIVW_X2]
dumTab2[,type := "slope"]

dumTab4 = rbind(dumTab,dumTab2)
dumTab4 = dumTab4[!is.na(mean_betaIVW_X1),]

data_hlines = data.frame(type = c(rep("mean",4),"slope"),
                         mylines = c(0.3,1.2,-1.2,-0.3,0.3*4*17.5))

plot5 = ggplot(dumTab4[Sim_X=="X12"], aes(x=Sim_Y, y=mean_betaIVW_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*sd_betaIVW_X1, ymax=mean_betaIVW_X1+1.96*sd_betaIVW_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab("Raw estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/RawEstimates_X12.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

plot5 = ggplot(dumTab4[Sim_X=="X123"], aes(x=Sim_Y, y=mean_betaIVW_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*sd_betaIVW_X1, ymax=mean_betaIVW_X1+1.96*sd_betaIVW_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab("Raw estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/RawEstimates_X123.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

plot5 = ggplot(dumTab4[Sim_X=="X13"], aes(x=Sim_Y, y=mean_betaIVW_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*sd_betaIVW_X1, ymax=mean_betaIVW_X1+1.96*sd_betaIVW_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab("Raw estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/RawEstimates_X13.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

#' # Check 5: conditional F-statistics ####
#' ***
dumTab = copy(myTab)
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,condFStats_median_X1 := condFStats_median_X2]
dumTab2[,condFStats_1stQ_X1 := condFStats_1stQ_X2]
dumTab2[,condFStats_3rdQ_X1 := condFStats_3rdQ_X2]
dumTab2[,type := "slope"]

dumTab4 = rbind(dumTab,dumTab2)
dumTab4 = dumTab4[!is.na(condFStats_median_X1),]

dumTab4[,dumID := paste(Sim_X,Sim_Y, sep=" - ")]

plot5 = ggplot(dumTab4[outcome=="Y2"], aes(x=dumID, y=condFStats_median_X1)) +
  facet_wrap(~ type,scales = "free") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(yintercept = 10,color="red",linetype = "dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=condFStats_1stQ_X1, 
                    ymax=condFStats_3rdQ_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Conditional F-Statistics")
plot5

filename = paste0("../results/_figures/CondFStats.png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot5)
dev.off()

#' # Save ####
#' ***
#' I want to save the big table as for the supplemental data!
#' 
save(myTab,file="../results/_tables/Simulation_complete.RData")

# excel_fn = paste0("../results/_tables/Simulation_complete.xlsx")
# WriteXLS("myTab", 
#          ExcelFileName=excel_fn, 
#          SheetNames="Sim_res", 
#          AutoFilter=T, 
#          BoldHeaderRow=T,
#          FreezeRow=1)

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
