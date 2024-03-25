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
ToDoFile[,age := dummy[seq(2,length(dummy),6)]]
ToDoFile[,SNP := dummy[seq(3,length(dummy),6)]]
ToDoFile[,growth := dummy[seq(4,length(dummy),6)]]
ToDoFile[,regression := dummy[seq(5,length(dummy),6)]]

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
  tab5 = SimTab[,median(NR_SNPs_type),by=dumID]
  matched = match(tab1$dumID,tab5$dumID)
  tab1[,NR_SNPs := tab5[matched,V1]]
  
  # add result to myRow
  tab1[,dumID := gsub("_bin","bin",dumID)]
  dummy = unlist(strsplit(tab1$dumID,"_"))
  tab1[,exposure := dummy[seq(1,length(dummy),3)]]
  tab1[,exposure_type := dummy[seq(2,length(dummy),3)]]
  tab1[,outcome := dummy[seq(3,length(dummy),3)]]
  
  tab1[,Sim_growth := myRow$growth]
  tab1[,Sim_age := myRow$age]
  tab1[,Sim_reg := myRow$regression]
  tab1[,Sim_SNPset := myRow$SNP]
  
  tab1_wide = dcast(tab1, Sim_growth + Sim_age + Sim_reg + Sim_SNPset + exposure + outcome ~ exposure_type, 
                    value.var = names(tab1)[2:11])
  tab1_wide
  
}

myTab = rbindlist(dumTab, fill = T)
myTab[,exposure := "continuous"]
myTab[grepl("bin",outcome),exposure := "binary"]
setnames(myTab,"exposure","outcome_type")

names(myTab)
myTab = myTab[,c(1:6,
                 10,13,16,19,22,25,28,31,34,
                 11,14,17,20,23,26,29,32,35,
                 12,15,18,21,24,27,30,33,36)]
myTab

#' # Check 1: Detection rates ####
#' ***
#' 
ToDoList = data.table(age = c("scan","scan","age","age"),
                      outcome_type = c("continuous","binary","continuous","binary"))

for(i in 1:dim(ToDoList)[1]){
  #i=1
  dumTab = copy(myTab)
  dumTab = dumTab[Sim_age==ToDoList[i,age],]
  dumTab = dumTab[outcome_type==ToDoList[i,outcome_type],]
  
  dumTab[,dumID1 := paste(Sim_SNPset,Sim_growth,sep="_")]
  dumTab[,dumID2 := paste(outcome,Sim_reg,sep="_")]
  
  dumTab_X1 <- dcast(dumTab, dumID2 ~ dumID1, value.var="N_sig_proc_X1")
  dumTab_X2 <- dcast(dumTab, dumID2 ~ dumID1, value.var="N_sig_proc_X2")
  
  dumTab2 = cbind(dumTab_X1,dumTab_X2[,c(2:5)])
  names(dumTab2)[2:9] = paste(names(dumTab2)[2:9],rep(c("X1","X2"),each=4),sep="_") 
  dumMat = as.matrix(dumTab2[,2:9])
  
  rownames(dumMat) = gsub("_"," - ",dumTab2$dumID2)
  colnames(dumMat) = gsub("_"," - ",colnames(dumMat))
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  
  filename = paste0("../results/_figures/DetectionRates_",ToDoList[i,age],"_",ToDoList[i,outcome_type],".png")
  png(filename = filename,width = 1600, height = 1600, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),col = COL1('Reds'), addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
}

#' **Summary**
#' 
#' - lin vs. quad: no difference in detection rates
#' - scan vs age: using real age increases the power for slope effect, but lowers it for mean effect
#' - dist vs. shared: works ok for Y2 (only effect of the mean), but cannot detect Y3 (only slope) or Y4 (mean and slope) correctly!
#' - linMixed vs gamlssIA: small differences, but linMixed seems to perform similar even without adjusting for heteroscedacicity
#' 
#' 
#' # Check 2: Bias ####
#' ***
#' I want to plot the bias - again per outcome!
#' 
for(i in 1:dim(ToDoList)[1]){
  #i=1
  dumTab = copy(myTab)
  dumTab = dumTab[Sim_age==ToDoList[i,age],]
  dumTab = dumTab[outcome_type==ToDoList[i,outcome_type],]
  
  dumTab[,type := "mean"]
  
  dumTab2 = copy(dumTab)
  dumTab2[,bias_X1 := bias_X2]
  dumTab2[,bias_SE_X1 := bias_SE_X2]
  dumTab2[,type := "slope"]
  
  dumTab4 = rbind(dumTab,dumTab2)
  dumTab4 = dumTab4[!is.na(bias_X1),]
  
  dumTab4[,dumID := paste(Sim_SNPset,Sim_growth, Sim_reg, sep=" - ")]
  
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
  
  filename = paste0("../results/_figures/Bias_Scan_",ToDoList[i,age],"_",ToDoList[i,outcome_type],".png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
}

#' # Check 3: Estimate ####
#' ***
#' I want to plot the corrected causal effect estimates - again per outcome!
#' 
for(i in 1:dim(ToDoList)[1]){
  #i=1
  dumTab = copy(myTab)
  dumTab = dumTab[Sim_age==ToDoList[i,age],]
  dumTab = dumTab[outcome_type==ToDoList[i,outcome_type],]
  
  dumTab[,type := "mean"]
  
  dumTab2 = copy(dumTab)
  dumTab2[,mean_betaIVW2_X1 := mean_betaIVW2_X2]
  dumTab2[,sd_betaIVW2_X1 := sd_betaIVW2_X2]
  dumTab2[,type := "slope"]
  
  dumTab4 = rbind(dumTab,dumTab2)
  dumTab4 = dumTab4[!is.na(mean_betaIVW2_X1),]
  
  dumTab4[,dumID := paste(Sim_SNPset,Sim_growth, Sim_reg, sep=" - ")]
  
  plot5 = ggplot(dumTab4, 
                 aes(x=dumID, y=mean_betaIVW2_X1,color = outcome)) +
    facet_wrap(~ type,scales = "free") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_betaIVW2_X1-1.96*sd_betaIVW2_X1, 
                      ymax=mean_betaIVW2_X1+1.96*sd_betaIVW2_X1), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) + 
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #theme(axis.text.x = element_text(angle = 45)) +
    xlab("Scenario") + ylab("Corrected estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0("../results/_figures/CorrectedEstimates_Scan_",ToDoList[i,age],"_",ToDoList[i,outcome_type],".png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
}

#' # Check 4: raw estimate ####
#' ***
#' I want to plot the uncorrected causal effect estimates - again per outcome!

for(i in 1:dim(ToDoList)[1]){
  #i=1
  dumTab = copy(myTab)
  dumTab = dumTab[Sim_age==ToDoList[i,age],]
  dumTab = dumTab[outcome_type==ToDoList[i,outcome_type],]
  
  dumTab[,type := "mean"]
  
  dumTab2 = copy(dumTab)
  dumTab2[,mean_betaIVW_X1 := mean_betaIVW_X2]
  dumTab2[,sd_betaIVW_X1 := sd_betaIVW_X2]
  dumTab2[,type := "slope"]

  dumTab4 = rbind(dumTab,dumTab2)
  dumTab4 = dumTab4[!is.na(mean_betaIVW_X1),]
  
  dumTab4[,dumID := paste(Sim_SNPset,Sim_growth, Sim_reg, sep=" - ")]
  
  plot5 = ggplot(dumTab4, 
                 aes(x=dumID, y=mean_betaIVW_X1,color = outcome)) +
    facet_wrap(~ type,scales = "free") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*sd_betaIVW_X1, 
                      ymax=mean_betaIVW_X1+1.96*sd_betaIVW_X1), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) + 
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    #theme(axis.text.x = element_text(angle = 45)) +
    xlab("Scenario") + ylab("Raw estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0("../results/_figures/RawEstimates_Scan_",ToDoList[i,age],"_",ToDoList[i,outcome_type],".png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
}

#' # Save ####
#' ***
#' I want to save the big table as for the supplemental data!
#' 
save(myTab,file="../results/_tables/Simulation_complete.RData")

excel_fn = paste0("../results/_tables/Simulation_complete.xlsx")
WriteXLS("myTab", 
         ExcelFileName=excel_fn, 
         SheetNames="Sim_res", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
