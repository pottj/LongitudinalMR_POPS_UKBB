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

source("../../SourceFile_HPC.R")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Get to-do list ####
#' ***
#' get the file names (same names for all results directories)
mySims = list.files(path = "../results/",pattern = "_SimAll.RData")
mySims

ToDoFile = data.table(NR = 1:length(mySims),
                      file = mySims)
dummy = unlist(strsplit(mySims,"_"))
ToDoFile[,tag := dummy[seq(2,length(dummy),4)]]
ToDoFile[,SNP := dummy[seq(3,length(dummy),4)]]

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

  # filter bad MR runs (cond. F-Stat <0 for mean and slope)
  # SimTab[,dumID0 := paste(n_sim,type,sep="__")]
  # dum1 = SimTab[exposure=="X1",mean(condFStats),by=dumID0]
  # dum2 = SimTab[exposure=="X2",mean(condFStats),by=dumID0]
  # dum3 = SimTab[exposure=="X3",mean(condFStats),by=dumID0]
  # dum1[,dumID1 := dumID0]
  # dum1[grepl("Only",dumID0),dumID1 := gsub("mean","",dumID0)]
  # dum2[,dumID1 := dumID0]
  # dum2[grepl("Only",dumID0),dumID1 := gsub("Slope","",dumID0)]
  # dum3[,dumID1 := dumID0]
  # dum3[grepl("Only",dumID0),dumID1 := gsub("Var","",dumID0)]
  # dum1 = dum1[!is.na(V1) & V1>0,]
  # dum2 = dum2[!is.na(V1) & V1>0,]
  # dum3 = dum3[!is.na(V1) & V1>0,]
  # mySims = c(dum1[dumID1 %in% dum2$dumID1,dumID0],dum2[dumID1 %in% dum1$dumID1,dumID0])
  # SimTab = SimTab[dumID0 %in% mySims,]
  # SimTab[,table(type)]
  
  x=length(unique(SimTab$n_sim))
  if(x<5){
    tab1_wide = data.table(Sim_growth = myRow$growth,
                           Sim_age = myRow$age,
                           Sim_reg = myRow$regression,
                           Sim_SNPset = myRow$SNP,
                           Sim_comment = "not enough simulations with cond FStat>10")
    
  }else{
    # add true value
    SimTab[outcome %in% c("Y1","Y3","Y1_bin","Y3_bin"),theta1 := 0]
    SimTab[outcome %in% c("Y2","Y4","Y2_bin","Y4_bin"),theta1 := 0.3]
    SimTab[outcome %in% c("Y1","Y2","Y1_bin","Y2_bin"),theta2 := 0]
    SimTab[outcome %in% c("Y3","Y4","Y3_bin","Y4_bin"),theta2 := 0.3]
    SimTab[exposure == "X1", theta := theta1]
    SimTab[exposure == "X2", theta := theta2]
    SimTab[exposure == "X3", theta := 0]
    SimTab[,theta1 := NULL]
    SimTab[,theta2 := NULL]
    
    # get significant simulations
    SimTab[,dumID := paste(type, exposure,outcome,sep="_")]
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
    tab1[,dumID := gsub("_bin","bin",dumID)]
    dummy = unlist(strsplit(tab1$dumID,"_"))
    tab1[,exposure := dummy[seq(3,length(dummy),4)]]
    tab1[,outcome := dummy[seq(4,length(dummy),4)]]
    tab1[,Sim_SNPset := dummy[seq(2,length(dummy),4)]]
    tab1[,Sim_type := dummy[seq(1,length(dummy),4)]]
    
    if(length(unique(SimTab$comment))==1){  
      tab1[,Sim_comment := unique(SimTab$comment)]
    }else{
      tab9 = table(SimTab$dumID,SimTab$comment)
      rownames(tab9) = gsub("_bin","bin",rownames(tab9))
      tab9_text = paste0(colnames(tab9)[1],": ",tab9[,1]," | ",colnames(tab9)[2],": ",tab9[,2])
      matched =  match(tab1$dumID,rownames(tab9))
      stopifnot(rownames(tab9)[matched]==tab1$dumID)
      tab1[,Sim_comment := tab9_text[matched]]
    }
    
    tab1_wide = dcast(tab1, Sim_type + Sim_SNPset + Sim_comment + outcome ~ exposure, 
                      value.var = names(tab1)[2:17])
  }
  
  tab1_wide
  
}

myTab = rbindlist(dumTab, fill = T)
myTab[,table(Sim_comment)]
myTab[,outcome_type := "continuous"]
myTab[grepl("bin",outcome),outcome_type := "binary"]

names(myTab)
myTab = myTab[,c(1:5,53,
                 11,14,17,20,23,26,29,32,35,38,41,44,47,50,
                 12,15,18,21,24,27,30,33,36,39,42,45,48,51,
                 13,16,19,22,25,28,31,34,37,40,43,46,49,52)]
myTab
myTab[grepl("p-value threshold = 0.05: 0",Sim_comment),Sim_comment := "p-value threshold = 1e-06"]

myTab = myTab[Sim_type != "sens6",]
myTab[Sim_type == "main",Sim_type := "Main"]
myTab[Sim_type == "sens1",Sim_type := "S1 - wrong time"]
myTab[Sim_type == "sens2",Sim_type := "S2 - wrong growth"]
myTab[Sim_type == "sens3",Sim_type := "S3 - no slope"]
myTab[Sim_type == "sens4",Sim_type := "S4 - var time indep."]
myTab[Sim_type == "sens5",Sim_type := "S5 - no var"]

#' # Check 1: Detection rates ####
#' ***
#' 
myOutcome_types = c("continuous","binary")

for(i in 1:length(myOutcome_types)){
  #i=1
  dumTab = copy(myTab)
  dumTab = dumTab[outcome_type==myOutcome_types[i],]
  
  dumTab[,dumID1 := paste(Sim_type,Sim_SNPset,sep="_")]
  
  dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X1")
  dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X2")
  dumTab_X3 <- dcast(dumTab, outcome ~ dumID1, value.var="N_sig_proc_X3")
  
  dumTab2 = rbind(dumTab_X1,dumTab_X2,dumTab_X3)
  setorder(dumTab2,outcome)
  dumMat = as.matrix(dumTab2[,2:13])
  rownames(dumMat) = paste(dumTab2$outcome,rep(c("X1","X2","X3"),3),sep=" - ")
  dumMat = dumMat[4:12,]
  colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45,na.label = " ")
  
  filename = paste0("../results/_figures/DetectionRates_condFStat0_",myOutcome_types[i],".png")
  png(filename = filename,width = 1600, height = 1600, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45,na.label = " ")
  dev.off()
}

#' # Check 2: Bias ####
#' ***
#' I want to plot the bias - again per outcome but only for continuous ones (does not really make sense for binary outcomes)
#' 
dumTab = copy(myTab)
dumTab = dumTab[outcome_type==myOutcome_types[1],]
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,bias_X1 := bias_X2]
dumTab2[,bias_SE_X1 := bias_SE_X2]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab)
dumTab3[,bias_X1 := bias_X3]
dumTab3[,bias_SE_X1 := bias_SE_X3]
dumTab3[,type := "var"]

dumTab4 = rbind(dumTab,dumTab2,dumTab3)
dumTab4 = dumTab4[!is.na(bias_X1),]

dumTab4[,dumID := paste(Sim_type,Sim_SNPset, sep=" - ")]

plot5 = ggplot(dumTab4, 
               aes(x=dumID, y=bias_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free") +
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

filename = paste0("../results/_figures/Bias_condFStat0_",myOutcome_types[1],".png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot5)
dev.off()

#' # Check 3: Estimate ####
#' ***
#' I want to plot the corrected causal effect estimates - again per outcome!
#' 
dumTab = copy(myTab)
dumTab = dumTab[outcome_type==myOutcome_types[1],]
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,mean_betaIVW2_X1 := mean_betaIVW2_X2]
dumTab2[,sd_betaIVW2_X1 := sd_betaIVW2_X2]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab)
dumTab3[,mean_betaIVW2_X1 := mean_betaIVW2_X3]
dumTab3[,sd_betaIVW2_X1 := sd_betaIVW2_X3]
dumTab3[,type := "var"]

dumTab4 = rbind(dumTab,dumTab2,dumTab3)
dumTab4 = dumTab4[!is.na(mean_betaIVW2_X1),]

dumTab4[,dumID := paste(Sim_type,Sim_SNPset, sep=" - ")]

plot5 = ggplot(dumTab4, 
               aes(x=dumID, y=mean_betaIVW2_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free") +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW2_X1-1.96*sd_betaIVW2_X1, ymax=mean_betaIVW2_X1+1.96*sd_betaIVW2_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Corrected estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/CorrectedEstimates_condFStat0_",myOutcome_types[1],".png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot5)
dev.off()

#' # Check 4: raw estimate ####
#' ***
#' I want to plot the uncorrected causal effect estimates - again per outcome!
dumTab = copy(myTab)
dumTab = dumTab[outcome_type==myOutcome_types[1],]
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,mean_betaIVW_X1 := mean_betaIVW_X2]
dumTab2[,sd_betaIVW_X1 := sd_betaIVW_X2]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab)
dumTab3[,mean_betaIVW_X1 := mean_betaIVW_X3]
dumTab3[,sd_betaIVW_X1 := sd_betaIVW_X3]
dumTab3[,type := "var"]

dumTab4 = rbind(dumTab,dumTab2,dumTab3)
dumTab4 = dumTab4[!is.na(mean_betaIVW_X1),]

dumTab4[,dumID := paste(Sim_type,Sim_SNPset, sep=" - ")]

plot5 = ggplot(dumTab4, 
               aes(x=dumID, y=mean_betaIVW_X1, color = outcome)) +
  facet_wrap(~ type,scales = "free") +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*sd_betaIVW_X1, ymax=mean_betaIVW_X1+1.96*sd_betaIVW_X1), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  #theme(axis.text.x = element_text(angle = 45)) +
  xlab("Scenario") + ylab("Raw estimate") +
  labs(color = "Outcome")
plot5

filename = paste0("../results/_figures/RawEstimates_condFStat0_",myOutcome_types[1],".png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot5)
dev.off()

#' # Check 5: conditional F-statistics ####
#' ***
dumTab = copy(myTab)
dumTab = dumTab[outcome_type==myOutcome_types[1],]
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,condFStats_median_X1 := condFStats_median_X2]
dumTab2[,condFStats_1stQ_X1 := condFStats_1stQ_X2]
dumTab2[,condFStats_3rdQ_X1 := condFStats_3rdQ_X2]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab)
dumTab3[,condFStats_median_X1 := condFStats_median_X3]
dumTab3[,condFStats_1stQ_X1 := condFStats_1stQ_X3]
dumTab3[,condFStats_3rdQ_X1 := condFStats_3rdQ_X3]
dumTab3[,type := "var"]

dumTab4 = rbind(dumTab,dumTab2,dumTab3)
dumTab4 = dumTab4[!is.na(condFStats_median_X1),]

dumTab4[,dumID := paste(Sim_type,Sim_SNPset, sep=" - ")]

plot5 = ggplot(dumTab4[outcome=="Y2"], 
               aes(x=dumID, y=condFStats_median_X1)) +
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

filename = paste0("../results/_figures/CondFStats_condFStat0_",myOutcome_types[1],".png")
png(filename = filename,width = 3200, height = 1600, res=200)
print(plot5)
dev.off()

#' # Save ####
#' ***
#' I want to save the big table as for the supplemental data!
#' 
save(myTab,file="../results/_tables/Simulation_complete_F0.RData")

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
