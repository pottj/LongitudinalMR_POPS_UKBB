#' ---
#' title: "Evaluation Sensitivity Runs"
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

#' # Get data ####
#' ***
#' Index j: all 10 sensitivity runs
#' Index i: all 12 settings per run
#' 
mySensitivityRuns = list.dirs(path = "../..",recursive = F)
mySensitivityRuns = mySensitivityRuns[-1]

myNames = gsub("../../","",mySensitivityRuns)

dumTab0 = foreach(j = 1:length(myNames))%do%{
  #j=7
  message("Working on Sensitivity run ",myNames[j]," (",j," of ",length(myNames),")")
  outdir_results = paste0("../results/_figures_",myNames[j])
  if(dir.exists(outdir_results)==F){
    dir.create(outdir_results)
    message("Created results folder ",outdir_results, " for simulation ",myNames[j])
  }else{
    message("Using pre-existing results folder ",outdir_results, " for simulation ",myNames[j])
  }
  
  # Find the 12 simulation data sets
  mySims = list.files(path = paste0(mySensitivityRuns[j],"/results/"),pattern = "_SimAll.RData")
  ToDoFile = data.table(NR = 1:length(mySims),
                        file = mySims)
  dummy = unlist(strsplit(mySims,"_"))
  ToDoFile[,SimNR := dummy[seq(2,length(dummy),5)]]
  ToDoFile[,SimX := dummy[seq(3,length(dummy),5)]]
  ToDoFile[,SimY := dummy[seq(4,length(dummy),5)]]
  
  ToDoFile[SimY == "CM1",theta1 := 0.3]
  ToDoFile[SimY == "CM2",theta1 := 1.2]
  ToDoFile[SimY == "CM3",theta1 := -1.2]
  ToDoFile[SimY == "CM4",theta1 := -0.3]
  ToDoFile[,theta2 := 0.3]
  ToDoFile[,theta3 := 1]
  
  # Loop to load all the data sets
  dumTab = foreach(i = 1:dim(ToDoFile)[1])%do%{
    #i=10
    myRow = ToDoFile[i,]
    message("Working on ",myRow$file)
    load(paste0(mySensitivityRuns[j],"/results/",myRow$file))
    
    # restrict to the first 100 simulations 
    SimTab = SimTab[n_sim <=100,]
    
    SimTab[exposure %in% c("mean"),exposure := "X1"]
    SimTab[exposure %in% c("slope"),exposure := "X2"]
    SimTab[exposure %in% c("var"),exposure := "X3"]
    
    # correct estimates for outcome age
    if(myNames[j] %in% c("samplesAsPOPS","weakInstrumentBias")){
      SimTab[exposure == "X2",beta_IVW := beta_IVW/30]
      SimTab[exposure == "X2",SE_IVW := SE_IVW/30]
    }else{
      SimTab[exposure == "X2",beta_IVW := beta_IVW/70]
      SimTab[exposure == "X2",SE_IVW := SE_IVW/70]
    }
    
    
    # correct estimates for Allele score factor
    SimTab[exposure == "X3",beta_IVW := beta_IVW*0.5]
    SimTab[exposure == "X3",SE_IVW := SE_IVW*0.5]
    
    # add true value
    SimTab[,theta1 := 0]
    SimTab[outcome %in% c("Y2","Y5","Y6","Y8"),theta1 := myRow$theta1]
    SimTab[,theta2 := 0]
    SimTab[outcome %in% c("Y3","Y5","Y7","Y8"),theta2 := myRow$theta2]
    SimTab[,theta3 := 0]
    SimTab[outcome %in% c("Y4","Y6","Y7","Y8"),theta3 := myRow$theta3]
    
    # finalize theta
    SimTab[exposure == "X1", theta := theta1]
    SimTab[exposure == "X2", theta := theta2]
    SimTab[exposure == "X3", theta := theta3]
    SimTab[,theta1 := NULL]
    SimTab[,theta2 := NULL]
    SimTab[,theta3 := NULL]
    
    # get power (theta != 0) / type I error (theta == 0) 
    SimTab[,dumID := paste(exposure,outcome,sep="_")]
    tab1 = SimTab[,.N,by = dumID]
    tab2 = SimTab[pval_IVW<0.05,.N,by = dumID]
    matched = match(tab1$dumID,tab2$dumID)
    tab1[,N_sig := tab2[matched,N]]
    tab1[is.na(N_sig),N_sig := 0]
    tab1[,power := N_sig/N ]
    tab1[,power_SE := sqrt((power * (1-power))/N)]
    
    # get coverage
    SimTab[,lower := beta_IVW - 1.96*SE_IVW]
    SimTab[,upper := beta_IVW + 1.96*SE_IVW]
    SimTab[,covered := 0]
    SimTab[lower<=theta & theta<= upper,covered := 1]
    tab2 = SimTab[,sum(covered),by = dumID]
    tab1[,coverage := tab2$V1/N]
    tab1[,coverage_SE := sqrt((coverage * (1-coverage)) / N)]
    
    # get bias 
    SimTab[,dif := beta_IVW - theta]
    tab3 = SimTab[,mean(dif),by = dumID]
    tab1[,bias := tab3[,V1]]
    tab3 = SimTab[,sd(beta_IVW),by = dumID]
    tab1[,bias_SE := tab3[,V1]]
    tab1[,bias_SE := bias_SE * sqrt(1/N)]
    
    # get empirical SE
    tab1[,empSE := tab3[,V1]]
    tab1[,empSE_SE := empSE / sqrt(2*(N-1))]
    
    # get mean of estimates
    tab4 = SimTab[,mean(beta_IVW),by = dumID]
    tab1[,mean_betaIVW := tab4[,V1]]
    
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
                      value.var = names(tab1)[2:19])
    tab1_wide
    
  }
  
  myTab = rbindlist(dumTab, fill = T)
  names(myTab)
  
  # in case of "noSlope" and "noVar", there is a different number of columns
  if(myNames[j] %in% c("noSlope","noVariability","noSlopeCorrelatedShared")){
    myTab = myTab[,c(1:4,
                     7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,
                     8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]
    
  }else{
    myTab = myTab[,c(1:4,
                     8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,
                     9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,
                     10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58)]
  }
 
  # Check 1: Detection rates / power / type 1 error
  {
    dumTab = copy(myTab)
    dumTab[,dumID1 := paste(Sim_X,Sim_Y,sep="_")]
    
    if(myNames[j] == "noSlope" | myNames[j]=="noSlopeCorrelatedShared"){
      dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="power_X1")
      dumTab_X3 <- dcast(dumTab, outcome ~ dumID1, value.var="power_X3")
      dumTab2 = cbind(dumTab_X1,dumTab_X3[,-1])
      dumMat = as.matrix(dumTab2[,-1])
      colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X3"),each=12),sep=" - ")
    }else if(myNames[j] == "noVariability"){
      dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="power_X1")
      dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="power_X2")
      dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1])
      dumMat = as.matrix(dumTab2[,-1])
      colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X2"),each=12),sep=" - ")
    }else{
      dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="power_X1")
      dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="power_X2")
      dumTab_X3 <- dcast(dumTab, outcome ~ dumID1, value.var="power_X3")
      dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1],dumTab_X3[,-1])
      dumMat = as.matrix(dumTab2[,-1])
      colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X2","X3"),each=12),sep=" - ")
    }
    
    filt1 = grepl("X12_",colnames(dumMat))
    filt2 = grepl("X123_",colnames(dumMat))
    filt3 = grepl("X13_",colnames(dumMat))
    
    filename = paste0(outdir_results,"/Power_X12.png")
    png(filename = filename,width = 2600, height = 1400, res=200)
    corrplot(dumMat[,filt1], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
             col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
             addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
    dev.off()
    
    filename = paste0(outdir_results,"/Power_X123.png")
    png(filename = filename,width = 2600, height = 1400, res=200)
    corrplot(dumMat[,filt2], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
             col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
             addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
    dev.off()
    
    filename = paste0(outdir_results,"/Power_X13.png")
    png(filename = filename,width = 2600, height = 1400, res=200)
    corrplot(dumMat[,filt3], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
             col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
             addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
    dev.off()
    
  }
  
  # Check 2: Coverage
  {
    dumTab = copy(myTab)
    dumTab[,dumID1 := paste(Sim_X,Sim_Y,sep="_")]
    
    if(myNames[j] == "noSlope" | myNames[j] == "noSlopeCorrelatedShared"){
      dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_X1")
      dumTab_X3 <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_X3")
      dumTab2 = cbind(dumTab_X1,dumTab_X3[,-1])
      dumMat = as.matrix(dumTab2[,-1])
      colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X3"),each=12),sep=" - ")
    }else if(myNames[j] == "noVariability"){
      dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_X1")
      dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_X2")
      dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1])
      dumMat = as.matrix(dumTab2[,-1])
      colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X2"),each=12),sep=" - ")
    }else{
      dumTab_X1 <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_X1")
      dumTab_X2 <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_X2")
      dumTab_X3 <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_X3")
      dumTab2 = cbind(dumTab_X1,dumTab_X2[,-1],dumTab_X3[,-1])
      dumMat = as.matrix(dumTab2[,-1])
      colnames(dumMat) = paste(colnames(dumMat),rep(c("X1","X2","X3"),each=12),sep=" - ")
    }
    
    filt1 = grepl("X12_",colnames(dumMat))
    filt2 = grepl("X123_",colnames(dumMat))
    filt3 = grepl("X13_",colnames(dumMat))
    
    filename = paste0(outdir_results,"/Coverage_X12.png")
    png(filename = filename,width = 2600, height = 1400, res=200)
    corrplot(dumMat[,filt1], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
             col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
             addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
    dev.off()
    
    filename = paste0(outdir_results,"/Coverage_X123.png")
    png(filename = filename,width = 2600, height = 1400, res=200)
    corrplot(dumMat[,filt2], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
             col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
             addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
    dev.off()
    
    filename = paste0(outdir_results,"/Coverage_X13.png")
    png(filename = filename,width = 2600, height = 1400, res=200)
    corrplot(dumMat[,filt3], is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
             col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
             addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
    dev.off()
    
  }
  
  # Check 3: Bias 
  {
    dumTab = copy(myTab)
    dumTab[,type := "mean"]
    
    if(!grepl("noSlope",myNames[j])){
      dumTab2 = copy(dumTab)
      dumTab2[,bias_X1 := bias_X2]
      dumTab2[,bias_SE_X1 := bias_SE_X2]
      dumTab2[,type := "slope"]
    }
    
    if(myNames[j] != "noVariability"){
      dumTab3 = copy(dumTab)
      dumTab3[,bias_X1 := bias_X3]
      dumTab3[,bias_SE_X1 := bias_SE_X3]
      dumTab3[,type := "var"]
      
    }
    
    if(grepl("noSlope",myNames[j])){
      dumTab4 = rbind(dumTab,dumTab3)
    }else if(myNames[j] == "noVariability"){ 
      dumTab4 = rbind(dumTab,dumTab2)
    }else{
      dumTab4 = rbind(dumTab,dumTab2,dumTab3)  
    }
    dumTab4 = dumTab4[!is.na(bias_X1),]
    
    plot5 = ggplot(dumTab4[Sim_X=="X12"], aes(x=Sim_Y, y=bias_X1, color = outcome)) +
      facet_wrap(~ type,scales = "free_y") +
      geom_hline(yintercept = 0,color="grey") +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=bias_X1-1.96*bias_SE_X1, ymax=bias_X1+1.96*bias_SE_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      xlab("Scenario") + ylab("Bias") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/Bias_X12.png")
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
      xlab("Scenario") + ylab("Bias") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/Bias_X123.png")
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
      xlab("Scenario") + ylab("Bias") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/Bias_X13.png")
    png(filename = filename,width = 2800, height = 1600, res=200)
    print(plot5)
    dev.off()
    
  }
  
  # Check 4: empSE
  {
    dumTab = copy(myTab)
    dumTab[,type := "mean"]
    
    if(!grepl("noSlope",myNames[j])){
      dumTab2 = copy(dumTab)
      dumTab2[,empSE_X1 := empSE_X2]
      dumTab2[,empSE_SE_X1 := empSE_SE_X2]
      dumTab2[,type := "slope"]
    }
    
    if(myNames[j] != "noVariability"){
      dumTab3 = copy(dumTab)
      dumTab3[,empSE_X1 := empSE_X3]
      dumTab3[,empSE_SE_X1 := empSE_SE_X3]
      dumTab3[,type := "var"]
    }
    
    if(grepl("noSlope",myNames[j])){
      dumTab4 = rbind(dumTab,dumTab3)
    }else if(myNames[j] == "noVariability"){ 
      dumTab4 = rbind(dumTab,dumTab2)
    }else{
      dumTab4 = rbind(dumTab,dumTab2,dumTab3)  
    }
    dumTab4 = dumTab4[!is.na(empSE_X1),]
    
    plot5 = ggplot(dumTab4[Sim_X=="X12"], aes(x=Sim_Y, y=empSE_X1, color = outcome)) +
      facet_wrap(~ type,scales = "free_y") +
      geom_hline(yintercept = 0,color="grey") +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=empSE_X1-1.96*empSE_SE_X1, ymax=empSE_X1+1.96*empSE_SE_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      #theme(axis.text.x = element_text(angle = 45)) +
      xlab("Scenario") + ylab("empirical SE") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/empSE_X12.png")
    png(filename = filename,width = 2800, height = 1600, res=200)
    print(plot5)
    dev.off()
    
    plot5 = ggplot(dumTab4[Sim_X=="X123"], aes(x=Sim_Y, y=empSE_X1, color = outcome)) +
      facet_wrap(~ type,scales = "free_y") +
      geom_hline(yintercept = 0,color="grey") +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=empSE_X1-1.96*empSE_SE_X1, ymax=empSE_X1+1.96*empSE_SE_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      #theme(axis.text.x = element_text(angle = 45)) +
      xlab("Scenario") + ylab("empirical SE") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/empSE_X123.png")
    png(filename = filename,width = 2800, height = 1600, res=200)
    print(plot5)
    dev.off()
    
    plot5 = ggplot(dumTab4[Sim_X=="X13"], aes(x=Sim_Y, y=empSE_X1, color = outcome)) +
      facet_wrap(~ type,scales = "free_y") +
      geom_hline(yintercept = 0,color="grey") +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=empSE_X1-1.96*empSE_SE_X1, ymax=empSE_X1+1.96*empSE_SE_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      #theme(axis.text.x = element_text(angle = 45)) +
      xlab("Scenario") + ylab("empirical SE") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/empSE_X13.png")
    png(filename = filename,width = 2800, height = 1600, res=200)
    print(plot5)
    dev.off()
    
  }
  
  # Check 5: Estimates
  {
    dumTab = copy(myTab)
    dumTab[,type := "mean"]
    
    if(!grepl("noSlope",myNames[j])){
      dumTab2 = copy(dumTab)
      dumTab2[,mean_betaIVW_X1 := mean_betaIVW_X2]
      dumTab2[,empSE_X1 := empSE_X2]
      dumTab2[,type := "slope"]
    }
    
    if(myNames[j] != "noVariability"){
      dumTab3 = copy(dumTab)
      dumTab3[,mean_betaIVW_X1 := mean_betaIVW_X3]
      dumTab3[,empSE_X1 := empSE_X3]
      dumTab3[,type := "var"]
    } 
    
    if(grepl("noSlope",myNames[j])){
      dumTab4 = rbind(dumTab,dumTab3)
    }else if(myNames[j] == "noVariability"){ 
      dumTab4 = rbind(dumTab,dumTab2)
    }else{
      dumTab4 = rbind(dumTab,dumTab2,dumTab3)  
    }
    dumTab4 = dumTab4[!is.na(mean_betaIVW_X1),]
    
    data_hlines = data.frame(type = c(rep("mean",4),"slope","var"),
                             mylines = c(0.3,1.2,-1.2,-0.3,0.3,1))
    
    plot5 = ggplot(dumTab4[Sim_X == "X12",], aes(x=Sim_Y, y=mean_betaIVW_X1, color = outcome)) +
      facet_wrap(~ type,scales = "free_y") +
      geom_hline(yintercept = 0,color="grey") +
      geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*empSE_X1, ymax=mean_betaIVW_X1+1.96*empSE_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      xlab("Scenario") + ylab("Corrected estimate") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/CorrectedEstimates_X12.png")
    png(filename = filename,width = 2800, height = 1600, res=200)
    print(plot5)
    dev.off()
    
    plot5 = ggplot(dumTab4[Sim_X == "X123",], aes(x=Sim_Y, y=mean_betaIVW_X1, color = outcome)) +
      facet_wrap(~ type,scales = "free_y") +
      geom_hline(yintercept = 0,color="grey") +
      geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*empSE_X1, ymax=mean_betaIVW_X1+1.96*empSE_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      xlab("Scenario") + ylab("Corrected estimate") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/CorrectedEstimates_X123.png")
    png(filename = filename,width = 2800, height = 1600, res=200)
    print(plot5)
    dev.off()
    
    plot5 = ggplot(dumTab4[Sim_X == "X13",], aes(x=Sim_Y, y=mean_betaIVW_X1, color = outcome)) +
      facet_wrap(~ type,scales = "free_y") +
      geom_hline(yintercept = 0,color="grey") +
      geom_hline(data = data_hlines, col="black", linetype="dotted", aes(yintercept = mylines)) +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=mean_betaIVW_X1-1.96*empSE_X1, ymax=mean_betaIVW_X1+1.96*empSE_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      xlab("Scenario") + ylab("Corrected estimate") +
      labs(color = "Outcome")
    plot5
    
    filename = paste0(outdir_results,"/CorrectedEstimates_X13.png")
    png(filename = filename,width = 2800, height = 1600, res=200)
    print(plot5)
    dev.off()
    
  }
  
  # Check 6: Conditional F-statistics
  {
    dumTab = copy(myTab)
    dumTab[,type := "mean"]
    
    if(!grepl("noSlope",myNames[j])){
      dumTab2 = copy(dumTab)
      dumTab2[,condFStats_median_X1 := condFStats_median_X2]
      dumTab2[,condFStats_1stQ_X1 := condFStats_1stQ_X2]
      dumTab2[,condFStats_3rdQ_X1 := condFStats_3rdQ_X2]
      dumTab2[,type := "slope"]
    }
    
    if(myNames[j] != "noVariability"){
      dumTab3 = copy(dumTab)
      dumTab3[,condFStats_median_X1 := condFStats_median_X3]
      dumTab3[,condFStats_1stQ_X1 := condFStats_1stQ_X3]
      dumTab3[,condFStats_3rdQ_X1 := condFStats_3rdQ_X3]
      dumTab3[,type := "var"]
    }

    if(grepl("noSlope",myNames[j])){
      dumTab4 = rbind(dumTab,dumTab3)
    }else if(myNames[j] == "noVariability"){ 
      dumTab4 = rbind(dumTab,dumTab2)
    }else{
      dumTab4 = rbind(dumTab,dumTab2,dumTab3)  
    }
    dumTab4 = dumTab4[!is.na(condFStats_median_X1),]
    
    dumTab4[,dumID := paste(Sim_X,Sim_Y, sep=" - ")]
    
    plot5 = ggplot(dumTab4[outcome=="Y2"], aes(x=Sim_Y, y=condFStats_median_X1, col = Sim_X)) +
      facet_wrap(~ type) +
      geom_hline(yintercept = 0,color="grey") +
      geom_hline(yintercept = 10,color="black",linetype = "dashed") +
      geom_point(position=position_dodge(0.5),size=3) +
      geom_errorbar(aes(ymin=condFStats_1stQ_X1, 
                        ymax=condFStats_3rdQ_X1), width=.2,
                    position=position_dodge(0.5)) +
      theme_bw(base_size = 15) + 
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      xlab("Scenario") + ylab("Conditional F-Statistics")
    plot5
    
    filename = paste0(outdir_results,"/CondFStats.png")
    png(filename = filename,width = 3200, height = 1600, res=200)
    print(plot5)
    dev.off()
    
  }
  
  # Save data
  save(myTab,file=paste0("../results/_tables/Simulation_complete_",myNames[j],".RData"))
  
} 
#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
