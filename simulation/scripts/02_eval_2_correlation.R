#' ---
#' title: "Check correlation in simulations"
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
suppressPackageStartupMessages(library(ggplot2))

#' # Get data ####
#' ***
variable_parameters = fread("../temp/parameters_variable.txt")
variable_parameters[,NR := 1:dim(variable_parameters)[1]]
myNames = paste(variable_parameters$NR,variable_parameters$scenario,sep="_")

fixed_parameters = fread("../temp/parameters_fixed.txt")
Y_theta = fixed_parameters[parameter %in% paste0("Y_theta_",c("M","S","V")),value]

dumTab0 = foreach(j = 1:length(myNames))%do%{
  #j=1
  message("Working on scenario ",myNames[j]," (",j," of ",length(myNames),")")
  
  #load correlation info
  myCors = list.files(path = paste0("../result/",myNames[j]),
                      pattern = "04_CorrelationSNPEffects",recursive = T)
  dumTab1 = foreach(i = 1:length(myCors))%do%{
    #i=1
    load(paste0("../result/",myNames[j],"/",myCors[i]))
    myCorTable[,dumID := paste(exposure,data1,data2,sep="_")]
    replicate_NR = gsub("Simulation_","",myCors[i])
    replicate_NR = gsub("/.*","",replicate_NR)
    myCorTable[,replicate_NR := as.numeric(replicate_NR)]
    myCorTable
  }
  myCorTab = rbindlist(dumTab1)
  
  # get power / type I error  
  tab1 = myCorTab[,.N,by = dumID]
  tab2 = myCorTab[pval<0.05,.N,by = dumID]
  matched = match(tab1$dumID,tab2$dumID)
  tab1[,N_sig := tab2[matched,N]]
  tab1[is.na(N_sig),N_sig := 0]
  tab1[,power := N_sig/N ]
  
  # get mean of estimates
  tab3 = myCorTab[,sd(cor),by = dumID]
  tab4 = myCorTab[,mean(cor),by = dumID]
  tab1[,mean_correlation := tab4[,V1]]
  tab1[,empSE := tab3[,V1]]
  
  # add result to myRow
  dummy = unlist(strsplit(tab1$dumID,"_"))
  tab1[,exposure := dummy[seq(1,length(dummy),5)]]
  tab1[,exposure_type1 := dummy[seq(3,length(dummy),5)]]
  tab1[,exposure_type2 := dummy[seq(5,length(dummy),5)]]
  tab1[,exposure_types := paste(exposure_type1,exposure_type2,sep=" vs. ")]
  
  tab1[,Sim_NR := variable_parameters$NR[j]]
  tab1[,Sim_name := variable_parameters$scenario[j]]
  
  tab1
}

myTab = rbindlist(dumTab0, fill=T)
names(myTab)

matchingTab = data.table(number1 = c(1:11),
                         number2 = c("0","2B","2A","1A","1B","3A","3B","4A","4B","5A","5B"))
matched = match(myTab$Sim_NR,matchingTab$number1)
table(is.na(matched))
myTab[,Sim_NR2 := matchingTab[matched,number2]]

#' # Check correlation per combination ####
#' ***

outdir_results = "../result/_figures_genCor/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
  message("Created figure folder ",outdir_results, " for genetic correlation")
}else{
  message("Using pre-existing figure folder ",outdir_results, " for genetic correlation ")
}
myTypes = unique(myTab$exposure_types)

# mean vs slope over all scenarios
plot5 = ggplot(myTab[exposure_types==myTypes[1]], 
               aes(x=Sim_NR2, y=mean_correlation, color = exposure)) +
  facet_wrap(~ exposure,scales = "free_y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(yintercept = -0.9,color="black",linetype="dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_correlation-1.96*empSE, 
                    ymax=mean_correlation+1.96*empSE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab(paste0("Genetic correlation, \n",myTypes[1])) +
  labs(color = "Exposure")
plot5

filename = paste0(outdir_results,"/Mean_Slope.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

# mean vs var over all scenarios
plot5 = ggplot(myTab[exposure_types==myTypes[2]], 
               aes(x=Sim_NR2, y=mean_correlation, color = exposure)) +
  facet_wrap(~ exposure,scales = "free_y") +
  geom_hline(yintercept = 0.5,color="grey",linetype="dotted") +
  geom_hline(yintercept = 0,color="black",linetype="dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_correlation-1.96*empSE, 
                    ymax=mean_correlation+1.96*empSE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab(paste0("Genetic correlation, \n",myTypes[2])) +
  labs(color = "Exposure")
plot5

filename = paste0(outdir_results,"/Mean_Var.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

# slope vs var over all scenarios
plot5 = ggplot(myTab[exposure_types==myTypes[3]], 
               aes(x=Sim_NR2, y=mean_correlation, color = exposure)) +
  facet_wrap(~ exposure,scales = "free_y") +
  geom_hline(yintercept = -0.3,color="grey",linetype="dotted") +
  geom_hline(yintercept = 0,color="black",linetype="dashed") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_correlation-1.96*empSE, 
                    ymax=mean_correlation+1.96*empSE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab(paste0("Genetic correlation, \n",myTypes[3])) +
  labs(color = "Exposure")
plot5

filename = paste0(outdir_results,"/Slope_Var.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

# all in one plot
data_hlines = data.frame(exposure_types = unique(myTab$exposure_types),
                         mylines1 = c(-0.9,0,0),
                         mylines2 = c(-0.9,0.5,-0.3))
plot5 = ggplot(myTab, 
               aes(x=Sim_NR2, y=mean_correlation, color = exposure)) +
  facet_wrap(~ exposure_types,scales = "free_y") +
  geom_hline(data = data_hlines, aes(yintercept = mylines1),
             linetype="dashed", show.legend = FALSE) +
  geom_hline(data = data_hlines, aes(yintercept = mylines2),
             linetype="dotted", show.legend = FALSE) +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=mean_correlation-1.96*empSE, 
                    ymax=mean_correlation+1.96*empSE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab("Genetic correlation") +
  labs(color = "Exposure")
plot5

filename = paste0(outdir_results,"/GenCor_all.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot5)
dev.off()

#' # Save data ####
#' ***
outdir_results = "../result/_tables/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
  message("Created table folder ",outdir_results, " for simulation summary")
}else{
  message("Using pre-existing table folder ",outdir_results, " for simulation summary")
}
save(myTab,file=paste0("../result/_tables/Simulation_GeneticCorrelation.RData"))

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
