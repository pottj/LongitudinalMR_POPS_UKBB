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

suppressPackageStartupMessages(library(data.table))
setDTthreads(1)
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggplot2))

#' # Get data ####
#' ***
variable_parameters = fread("../temp/parameters_variable.txt")
variable_parameters[,NR := 1:dim(variable_parameters)[1]]
myNames = variable_parameters$scenario

fixed_parameters = fread("../temp/parameters_fixed.txt")
Y_theta = fixed_parameters[parameter %in% paste0("Y_theta_",c("M","S","V")),value]

mySims = list.files(path = "../result/",pattern = "_SimAll.RData")
mySims_NR = gsub("_.*","",mySims)
mySims_NR = as.numeric(mySims_NR)
ordering = order(mySims_NR)
mySims = mySims[ordering]

dumTab0 = foreach(j = 1:length(mySims))%do%{
  #j=1
  message("Working on scenario ",myNames[j]," (",j," of ",length(myNames),")")
  load(paste0("../result/",mySims[j]))
  
  # rename exposure columns
  setnames(SimTab,"exposure","exposure_type")
  setnames(SimTab,"comment2","exposure")
  
  # add true theta value
  SimTab[,theta1 := 0]
  SimTab[outcome %in% c("Y2","Y5","Y6","Y8"),theta1 := Y_theta[1]]
  SimTab[,theta2 := 0]
  SimTab[outcome %in% c("Y3","Y5","Y7","Y8"),theta2 := Y_theta[2]]
  SimTab[,theta3 := 0]
  SimTab[outcome %in% c("Y4","Y6","Y7","Y8"),theta3 := Y_theta[3]]
  
  # finalize theta
  SimTab[exposure_type == "mean", theta := theta1]
  SimTab[exposure_type == "slope", theta := theta2]
  SimTab[exposure_type == "var", theta := theta3]
  SimTab[,theta1 := NULL]
  SimTab[,theta2 := NULL]
  SimTab[,theta3 := NULL]
  
  # get power (theta != 0) / type I error (theta == 0) 
  SimTab[,dumID := paste(exposure,exposure_type,outcome,sep="_")]
  tab1 = SimTab[,.N,by = dumID]
  tab2 = SimTab[pval<0.05,.N,by = dumID]
  matched = match(tab1$dumID,tab2$dumID)
  tab1[,N_sig := tab2[matched,N]]
  tab1[is.na(N_sig),N_sig := 0]
  tab1[,power := N_sig/N ]
  tab1[,power_SE := sqrt((power * (1-power))/N)]
  
  # get coverage
  SimTab[,lower := beta - 1.96*SE]
  SimTab[,upper := beta + 1.96*SE]
  SimTab[,covered := 0]
  SimTab[lower<=theta & theta<= upper,covered := 1]
  tab2 = SimTab[,sum(covered),by = dumID]
  tab1[,coverage := tab2$V1/N]
  tab1[,coverage_SE := sqrt((coverage * (1-coverage)) / N)]
  
  # get bias 
  SimTab[,dif := beta - theta]
  tab3 = SimTab[,mean(dif),by = dumID]
  tab1[,bias := tab3[,V1]]
  tab3 = SimTab[,sd(beta),by = dumID]
  tab1[,bias_SE := tab3[,V1]]
  tab1[,bias_SE := bias_SE * sqrt(1/N)]
  
  # get empirical SE
  tab1[,empSE := tab3[,V1]]
  tab1[,empSE_SE := empSE / sqrt(2*(N-1))]
  
  # get mean of estimates
  tab4 = SimTab[,mean(beta),by = dumID]
  tab1[,mean_beta := tab4[,V1]]
  
  # get median numbers of SNPs
  stopifnot(sum(is.na(SimTab$NR_SNPs_type))==0)
  tab5 = SimTab[,summary(NR_SNPs_type),by=dumID]
  tab5[,V2 := rep(c("min","1stQ","median","mean","3rdQ","max"),length(unique(dumID)))]
  tab5_median = copy(tab5)[V2=="median"]
  tab5_1stQ = copy(tab5)[V2=="1stQ"]
  tab5_3rdQ = copy(tab5)[V2=="3rdQ"]
  
  matched = match(tab1$dumID,tab5_median$dumID)
  tab1[,NR_SNPs_median := as.numeric(tab5_median[matched,V1])]
  # tab1[,NR_SNPs_1stQ := as.numeric(tab5_1stQ[matched,V1])]
  # tab1[,NR_SNPs_3rdQ := as.numeric(tab5_3rdQ[matched,V1])]
  
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
  # tab1[,condFstats_NA := as.numeric(tab6$V1)]
  
  # add result to myRow
  dummy = unlist(strsplit(tab1$dumID,"_"))
  tab1[,exposure := dummy[seq(1,length(dummy),3)]]
  tab1[,exposure_type := dummy[seq(2,length(dummy),3)]]
  tab1[,outcome := dummy[seq(3,length(dummy),3)]]
  
  tab1[,Sim_NR := variable_parameters$NR[j]]
  tab1[,Sim_name := variable_parameters$scenario[j]]
  
  tab1_wide = dcast(tab1, Sim_NR + Sim_name + exposure + outcome ~ exposure_type, 
                    value.var = names(tab1)[2:16])
  tab1_wide
}

myTab = rbindlist(dumTab0, fill=T)
names(myTab)
myTab = myTab[grepl("Y",outcome)]

matchingTab = data.table(number1 = c(1:11),
                         number2 = c("0","2B","2A","1A","1B","3A","3B","4A","4B","5A","5B"))
matched = match(myTab$Sim_NR,matchingTab$number1)
table(is.na(matched))
myTab[,Sim_NR2 := matchingTab[matched,number2]]

#' # Check 1: Detection rates ####
#' ***
#' Also known as power (if there should be an effect) or type 1 error (if there should not be an effect)
{
  outdir_results = "../result/_figures_power/"
  if(dir.exists(outdir_results)==F){
    dir.create(outdir_results)
    message("Created figure folder ",outdir_results, " for power heat maps")
  }else{
    message("Using pre-existing figure folder ",outdir_results, " for power heat maps ")
  }
  
  dumTab = copy(myTab)
  dumTab[,dumID1 := paste(exposure,Sim_NR2,sep="_")]
  
  # Power to detect mean effect over all scenarios
  dumTab_mean <- dcast(dumTab, outcome ~ dumID1, value.var="power_mean")
  dumMat = as.matrix(dumTab_mean[,-1])
  filename = paste0(outdir_results,"/Mean_allScenarios.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # Power to detect slope effect over all scenarios
  dumTab_slope <- dcast(dumTab, outcome ~ dumID1, value.var="power_slope")
  filt = grepl("_3A",names(dumTab_slope))
  dumTab_slope = dumTab_slope[,!filt,with = F]
  dumMat = as.matrix(dumTab_slope[,-1])
  filename = paste0(outdir_results,"/Slope_allScenarios.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # Power to detect variability effect over all scenarios
  dumTab_var <- dcast(dumTab, outcome ~ dumID1, value.var="power_var")
  filt = grepl("_3B",names(dumTab_var))
  dumTab_var = dumTab_var[,!filt,with = F]
  dumMat = as.matrix(dumTab_var[,-1])
  filename = paste0(outdir_results,"/Variability_allScenarios.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # Power to detect effects of MSV exposure over all scenarios and all exposure types
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("MSV",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:11] = gsub("MSV","mean",colnames(dumMat)[1:11])
  colnames(dumMat)[12:21] = gsub("MSV","slope",colnames(dumMat)[12:21])
  colnames(dumMat)[22:31] = gsub("MSV","var",colnames(dumMat)[22:31])
  filename = paste0(outdir_results,"/MSV_allScenarios_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # Power to detect effects of MS exposure over all scenarios and all exposure types
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("MS_",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:11] = gsub("MS","mean",colnames(dumMat)[1:11])
  colnames(dumMat)[12:21] = gsub("MS","slope",colnames(dumMat)[12:21])
  colnames(dumMat)[22:31] = gsub("MS","var",colnames(dumMat)[22:31])
  filename = paste0(outdir_results,"/MS_allScenarios_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # Power to detect effects of MV exposure over all scenarios and all exposure types
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("MV_",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:11] = gsub("MV","mean",colnames(dumMat)[1:11])
  colnames(dumMat)[12:21] = gsub("MV","slope",colnames(dumMat)[12:21])
  colnames(dumMat)[22:31] = gsub("MV","var",colnames(dumMat)[22:31])
  filename = paste0(outdir_results,"/MV_allScenarios_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
 
  # Power to detect effects in main scenario 
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("_0",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:3] = paste0(colnames(dumMat)[1:3],"_mean")
  colnames(dumMat)[4:6] = paste0(colnames(dumMat)[4:6],"_slope")
  colnames(dumMat)[7:9] = paste0(colnames(dumMat)[7:9],"_var")
  filename = paste0(outdir_results,"/Main_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off() 
}

#' # Check 2: Coverage ####
#' ***
#' Is the true effect within the confidence intervals?
{
  outdir_results = "../result/_figures_coverage/"
  if(dir.exists(outdir_results)==F){
    dir.create(outdir_results)
    message("Created figure folder ",outdir_results, " for Coverage heat maps")
  }else{
    message("Using pre-existing figure folder ",outdir_results, " for Coverage heat maps ")
  }
  
  dumTab = copy(myTab)
  dumTab[,dumID1 := paste(exposure,Sim_NR2,sep="_")]
  
  # mean effect over all scenarios
  dumTab_mean <- dcast(dumTab, outcome ~ dumID1, value.var="coverage_mean")
  dumMat = as.matrix(dumTab_mean[,-1])
  filename = paste0(outdir_results,"/Mean_allScenarios.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
           col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # slope effect over all scenarios
  dumTab_slope <- dcast(dumTab, outcome ~ dumID1, value.var="power_slope")
  filt = grepl("_6",names(dumTab_slope))
  dumTab_slope = dumTab_slope[,!filt,with = F]
  dumMat = as.matrix(dumTab_slope[,-1])
  filename = paste0(outdir_results,"/Slope_allScenarios.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
           col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # variability effect over all scenarios
  dumTab_var <- dcast(dumTab, outcome ~ dumID1, value.var="power_var")
  filt = grepl("_7",names(dumTab_var))
  dumTab_var = dumTab_var[,!filt,with = F]
  dumMat = as.matrix(dumTab_var[,-1])
  filename = paste0(outdir_results,"/Variability_allScenarios.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
           col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # effects of MSV exposure over all scenarios and all exposure types
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("MSV",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:11] = gsub("MSV","mean",colnames(dumMat)[1:11])
  colnames(dumMat)[12:21] = gsub("MSV","slope",colnames(dumMat)[12:21])
  colnames(dumMat)[22:31] = gsub("MSV","var",colnames(dumMat)[22:31])
  filename = paste0(outdir_results,"/MSV_allScenarios_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
           col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # effects of MS exposure over all scenarios and all exposure types
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("MS_",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:11] = gsub("MS","mean",colnames(dumMat)[1:11])
  colnames(dumMat)[12:21] = gsub("MS","slope",colnames(dumMat)[12:21])
  colnames(dumMat)[22:31] = gsub("MS","var",colnames(dumMat)[22:31])
  filename = paste0(outdir_results,"/MS_allScenarios_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
           col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # effects of MV exposure over all scenarios and all exposure types
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("MV_",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:11] = gsub("MV","mean",colnames(dumMat)[1:11])
  colnames(dumMat)[12:21] = gsub("MV","slope",colnames(dumMat)[12:21])
  colnames(dumMat)[22:31] = gsub("MV","var",colnames(dumMat)[22:31])
  filename = paste0(outdir_results,"/MV_allScenarios_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
           col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # effects in main scenario 
  dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
  filt = grepl("_0",names(dumTab2))
  dumTab2 = dumTab2[,filt,with = F]
  dumMat = as.matrix(dumTab2)
  colnames(dumMat)[1:3] = paste0(colnames(dumMat)[1:3],"_mean")
  colnames(dumMat)[4:6] = paste0(colnames(dumMat)[4:6],"_slope")
  colnames(dumMat)[7:9] = paste0(colnames(dumMat)[7:9],"_var")
  filename = paste0(outdir_results,"/Main_allTypes.png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'),
           col= colorRampPalette(c("darkblue","#FFF5F0","#67000D"))(5),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off() 
}

#' # Check 3: Bias ####
#' ***
#' Is the observed difference random or directed?
{
  outdir_results = "../result/_figures_bias/"
  if(dir.exists(outdir_results)==F){
    dir.create(outdir_results)
    message("Created figure folder ",outdir_results, " for bias scatter plot")
  }else{
    message("Using pre-existing figure folder ",outdir_results, " for bias scatter plot")
  }
  
  # create dummy table
  dumTab1 = copy(myTab)
  dumTab1[,type := "mean"]
  dumTab1[,bias := bias_mean]
  dumTab1[,bias_SE := bias_SE_mean]
  dumTab2 = copy(myTab)
  dumTab2[,type := "slope"]
  dumTab2[,bias := bias_slope]
  dumTab2[,bias_SE := bias_SE_slope]
  dumTab3 = copy(myTab)
  dumTab3[,type := "var"]
  dumTab3[,bias := bias_var]
  dumTab3[,bias_SE := bias_SE_var]
  dumTab4 = rbind(dumTab1,dumTab2,dumTab3)  
  dumTab4 = dumTab4[!is.na(bias),]
  dumTab4 = dumTab4[,c(1:4,50:53)]
  
  # mean effect over all scenarios
  plot5 = ggplot(dumTab4[type=="mean"], aes(x=Sim_NR2, y=bias, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Mean_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # slope effect over all scenarios
  plot5 = ggplot(dumTab4[type=="slope"], aes(x=Sim_NR2, y=bias, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Slope_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # variability effect over all scenarios
  plot5 = ggplot(dumTab4[type=="var"], aes(x=Sim_NR2, y=bias, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Var_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MSV exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MSV"], aes(x=Sim_NR2, y=bias, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MSV_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MS exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MS"], aes(x=Sim_NR2, y=bias, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MS_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MV exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MV"], aes(x=Sim_NR2, y=bias, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MV_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects in main scenario 
  plot5 = ggplot(dumTab4[Sim_NR2=="0"], aes(x=exposure, y=bias, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Main_allTypes.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
}

#' # Check 4: empSE ####
#' ***
#' Not sure what it tells me - look up again!
{
  outdir_results = "../result/_figures_empSE/"
  if(dir.exists(outdir_results)==F){
    dir.create(outdir_results)
    message("Created figure folder ",outdir_results, " for empSE scatter plot")
  }else{
    message("Using pre-existing figure folder ",outdir_results, " for empSE scatter plot")
  }
  
  # create dummy table
  dumTab1 = copy(myTab)
  dumTab1[,type := "mean"]
  dumTab1[,empSE := empSE_mean]
  dumTab1[,empSE_SE := empSE_SE_mean]
  dumTab2 = copy(myTab)
  dumTab2[,type := "slope"]
  dumTab2[,empSE := empSE_slope]
  dumTab2[,empSE_SE := empSE_SE_slope]
  dumTab3 = copy(myTab)
  dumTab3[,type := "var"]
  dumTab3[,empSE := empSE_var]
  dumTab3[,empSE_SE := empSE_SE_var]
  dumTab4 = rbind(dumTab1,dumTab2,dumTab3)  
  dumTab4 = dumTab4[!is.na(empSE),]
  dumTab4 = dumTab4[,c(1:4,50:53)]

  # mean effect over all scenarios
  plot5 = ggplot(dumTab4[type=="mean"], aes(x=Sim_NR2, y=empSE, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=empSE-1.96*empSE_SE, ymax=empSE+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("empSE") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Mean_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # slope effect over all scenarios
  plot5 = ggplot(dumTab4[type=="slope"], aes(x=Sim_NR2, y=empSE, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=empSE-1.96*empSE_SE, ymax=empSE+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("empSE") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Slope_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # variability effect over all scenarios
  plot5 = ggplot(dumTab4[type=="var"], aes(x=Sim_NR2, y=empSE, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=empSE-1.96*empSE_SE, ymax=empSE+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("empSE") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Var_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MSV exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MSV"], aes(x=Sim_NR2, y=empSE, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=empSE-1.96*empSE_SE, ymax=empSE+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("empSE") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MSV_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MS exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MS"], aes(x=Sim_NR2, y=empSE, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=empSE-1.96*empSE_SE, ymax=empSE+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("empSE") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MS_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MV exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MV"], aes(x=Sim_NR2, y=empSE, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=empSE-1.96*empSE_SE, ymax=empSE+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("empSE") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MV_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects in main scenario 
  plot5 = ggplot(dumTab4[Sim_NR2=="0"], aes(x=exposure, y=empSE, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=empSE-1.96*empSE_SE, ymax=empSE+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure") + ylab("empSE") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Main_allTypes.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
}

#' # Check 5: Estimates ####
#' ***
#' I just want to visualize the estimates too, but this is not a performance metric
{
  outdir_results = "../result/_figures_estimates/"
  if(dir.exists(outdir_results)==F){
    dir.create(outdir_results)
    message("Created figure folder ",outdir_results, " for estimates scatter plot")
  }else{
    message("Using pre-existing figure folder ",outdir_results, " for estimates scatter plot")
  }
  
  # create dummy table
  dumTab1 = copy(myTab)
  dumTab1[,type := "mean"]
  dumTab1[,mean_beta := mean_beta_mean]
  dumTab1[,empSE_SE := empSE_SE_mean]
  dumTab2 = copy(myTab)
  dumTab2[,type := "slope"]
  dumTab2[,mean_beta := mean_beta_slope]
  dumTab2[,empSE_SE := empSE_SE_slope]
  dumTab3 = copy(myTab)
  dumTab3[,type := "var"]
  dumTab3[,mean_beta := mean_beta_var]
  dumTab3[,empSE_SE := empSE_SE_var]
  dumTab4 = rbind(dumTab1,dumTab2,dumTab3)  
  dumTab4 = dumTab4[!is.na(mean_beta),]
  dumTab4 = dumTab4[,c(1:4,50:53)]
  
  # mean effect over all scenarios
  plot5 = ggplot(dumTab4[type=="mean"], aes(x=Sim_NR2, y=mean_beta, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = Y_theta[1],color="black",linetype = "dotted")+
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_beta-1.96*empSE_SE, ymax=mean_beta+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("MVMR estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Mean_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # slope effect over all scenarios
  plot5 = ggplot(dumTab4[type=="slope"], aes(x=Sim_NR2, y=mean_beta, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = Y_theta[2],color="black",linetype = "dotted")+
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_beta-1.96*empSE_SE, ymax=mean_beta+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("MVMR estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Slope_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # variability effect over all scenarios
  plot5 = ggplot(dumTab4[type=="var"], aes(x=Sim_NR2, y=mean_beta, color = outcome)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = Y_theta[3],color="black",linetype = "dotted")+
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_beta-1.96*empSE_SE, ymax=mean_beta+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("MVMR estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Var_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MSV exposure over all scenarios and all exposure types
  data_hlines = data.frame(type = unique(dumTab4$type),
                           mylines = Y_theta)
  plot5 = ggplot(dumTab4[exposure=="MSV"], aes(x=Sim_NR2, y=mean_beta, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(data = data_hlines, aes(yintercept = mylines),linetype="dashed", show.legend = FALSE) +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_beta-1.96*empSE_SE, ymax=mean_beta+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("MVMR estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MSV_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MS exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MS"], aes(x=Sim_NR2, y=mean_beta, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(data = data_hlines, aes(yintercept = mylines),linetype="dashed", show.legend = FALSE) +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_beta-1.96*empSE_SE, ymax=mean_beta+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("MVMR estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MS_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects of MV exposure over all scenarios and all exposure types
  plot5 = ggplot(dumTab4[exposure=="MV"], aes(x=Sim_NR2, y=mean_beta, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(data = data_hlines, aes(yintercept = mylines),linetype="dashed", show.legend = FALSE) +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_beta-1.96*empSE_SE, ymax=mean_beta+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("MVMR estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/MV_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # effects in main scenario 
  plot5 = ggplot(dumTab4[Sim_NR2=="0"], aes(x=exposure, y=mean_beta, color = outcome)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(data = data_hlines, aes(yintercept = mylines),linetype="dashed", show.legend = FALSE) +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=mean_beta-1.96*empSE_SE, ymax=mean_beta+1.96*empSE_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure") + ylab("MVMR estimate") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Main_allTypes.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
}

#' # Check 6: Conditional F-statistics ####
#' ***
#' Are my MVMRs well powered?
{
  outdir_results = "../result/_figures_condF/"
  if(dir.exists(outdir_results)==F){
    dir.create(outdir_results)
    message("Created figure folder ",outdir_results, " for estimates scatter plot")
  }else{
    message("Using pre-existing figure folder ",outdir_results, " for estimates scatter plot")
  }
  
  # create dummy table
  dumTab1 = copy(myTab)
  dumTab1[,type := "mean"]
  dumTab1[,condFStats_median := condFStats_median_mean]
  dumTab1[,condFStats_1stQ := condFStats_1stQ_mean]
  dumTab1[,condFStats_3rdQ := condFStats_3rdQ_mean]
  dumTab2 = copy(myTab)
  dumTab2[,type := "slope"]
  dumTab2[,condFStats_median := condFStats_median_slope]
  dumTab2[,condFStats_1stQ := condFStats_1stQ_slope]
  dumTab2[,condFStats_3rdQ := condFStats_3rdQ_slope]
  dumTab3 = copy(myTab)
  dumTab3[,type := "var"]
  dumTab3[,condFStats_median := condFStats_median_var]
  dumTab3[,condFStats_1stQ := condFStats_1stQ_var]
  dumTab3[,condFStats_3rdQ := condFStats_3rdQ_var]
  dumTab4 = rbind(dumTab1,dumTab2,dumTab3)  
  dumTab4 = dumTab4[!is.na(condFStats_median),]
  dumTab4 = dumTab4[,c(1:4,50:54)]
  dumTab4 = dumTab4[outcome=="Y8",]
  
  # mean effect over all scenarios
  plot5 = ggplot(dumTab4[type=="mean"], aes(x=exposure, y=condFStats_median,color=exposure)) +
    facet_wrap(~ Sim_NR2,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype = "dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, 
                      ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Conditional F-Statistics") +
    labs(color = "Exposure")
  plot5
  
  filename = paste0(outdir_results,"/Mean_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # slope effect over all scenarios
  plot5 = ggplot(dumTab4[type=="slope"], aes(x=exposure, y=condFStats_median,color=exposure)) +
    facet_wrap(~ Sim_NR2,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype = "dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, 
                      ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Conditional F-Statistics") +
    labs(color = "Exposure")
  plot5
  
  filename = paste0(outdir_results,"/Slope_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # variability effect over all scenarios
  plot5 = ggplot(dumTab4[type=="var"], aes(x=exposure, y=condFStats_median,color=exposure)) +
    facet_wrap(~ Sim_NR2,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype = "dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, 
                      ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Conditional F-Statistics") +
    labs(color = "Exposure")
  plot5
  
  filename = paste0(outdir_results,"/Var_allScenarios.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # MSV effect over all scenarios
  plot5 = ggplot(dumTab4[exposure=="MSV"], aes(x=Sim_NR2, y=condFStats_median,color=exposure)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype = "dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, 
                      ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Conditional F-Statistics") +
    labs(color = "Exposure")
  plot5
  
  filename = paste0(outdir_results,"/MSV_allScenarios_allTypes.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # MS effect over all scenarios
  plot5 = ggplot(dumTab4[exposure=="MS"], aes(x=Sim_NR2, y=condFStats_median,color=exposure)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype = "dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, 
                      ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Conditional F-Statistics") +
    labs(color = "Exposure")
  plot5
  
  filename = paste0(outdir_results,"/MS_allScenarios_allTypes.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # MV effect over all scenarios
  plot5 = ggplot(dumTab4[exposure=="MV"], aes(x=Sim_NR2, y=condFStats_median,color=exposure)) +
    facet_wrap(~ type,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype = "dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, 
                      ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Conditional F-Statistics") +
    labs(color = "Exposure")
  plot5
  
  filename = paste0(outdir_results,"/MV_allScenarios_allTypes.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
  # MV effect over all scenarios
  plot5 = ggplot(dumTab4[Sim_NR2=="0"], aes(x=type, y=condFStats_median,color=exposure)) +
    facet_wrap(~ exposure,scales = "free_y") +
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype = "dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, 
                      ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 15) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario") + ylab("Conditional F-Statistics") +
    labs(color = "Exposure")
  plot5
  
  filename = paste0(outdir_results,"/Main_allTypes.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot5)
  dev.off()
  
}

#' # Save data ####
#' ***
outdir_results = "../result/_tables/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
  message("Created table folder ",outdir_results, " for simulation summary")
}else{
  message("Using pre-existing table folder ",outdir_results, " for simulation summary")
}
save(myTab,file=paste0("../result/_tables/Simulation_complete.RData"))
  
#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
