#' ---
#' title: "Supplemental simulation figure: estimates"
#' subtitle: "Main Figures"
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

outdir_results = "../results/_figures/SupFigures/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
  message("Created results folder ",outdir_results, " for supplemental figures ")
}else{
  message("Using pre-existing results folder ",outdir_results, " for supplemental figures ")
}

#' # Get data ####
#' ***
myMainScenarios = list.files(path = "../results/_tables/",pattern = "main")
mynoSlopeScenarios = list.files(path = "../results/_tables/",pattern = "noSlope")
myScenarios = c(myMainScenarios,mynoSlopeScenarios)

dumTab = foreach(i=1:length(myScenarios))%do%{
  #i=1
  load(paste0("../results/_tables/",myScenarios[i]))
  myFlag = gsub("Simulation_complete_","",myScenarios[i])
  myFlag = gsub(".RData","",myFlag)
  myTab[,flag:=myFlag]
  myTab
}
myTab = rbindlist(dumTab, fill=T)
myTab = myTab[Sim_Y %in% c("CM2")]
myTab[,table(flag)]

#' # Get data ####
#' ***
dumTab = copy(myTab)
dumTab[,type := "mean"]

dumTab2 = copy(dumTab)
dumTab2[,mean_betaIVW_X1 := mean_betaIVW_X2]
dumTab2[,empSE_X1 := empSE_X2]
dumTab2[,type := "slope"]

dumTab3 = copy(dumTab)
dumTab3[,mean_betaIVW_X1 := mean_betaIVW_X3]
dumTab3[,empSE_X1 := empSE_X3]
dumTab3[,type := "variability"]

dumTab4 = rbind(dumTab,dumTab2,dumTab3)
dumTab4 = dumTab4[!is.na(mean_betaIVW_X1),]

# reduce to relevant columns
myNames = c("Sim_X", "flag", "outcome", "type", "mean_betaIVW_X1", "empSE_X1")
colsOut = setdiff(names(dumTab4),myNames)
dumTab4[,get("colsOut"):=NULL]
setcolorder(dumTab4,myNames)
setnames(dumTab4,"mean_betaIVW_X1","value")
setnames(dumTab4,"empSE_X1","SE")
#dumTab4[,type := paste0("estimate - ",type)]

estimateData = copy(dumTab4)

#' # Get Plots ####
#' ***
plotData = copy(estimateData)

data_hlines = data.frame(type = c("mean","slope","variability"),
                         mylines = c(1.2,0.3*70,1))


#' I want to have 6 facets - exposure type by GWAS model

plotData[,type2:= "main"]
plotData[grepl("noSlope",flag),type2:= "no slope"]

plotData[,flag2:= "A"]
plotData[grepl("B",flag),flag2:= "B"]

#' I want the raw estimates to better compare slope and variability
#' 
plotData[type == "slope",value := value * 70]
plotData[type == "slope",SE := SE * 70]

# plotData[flag2=="A" & type =="variability", value := value * 0.5]
# plotData[flag2=="A" & type =="variability", SE := SE * 0.5]
# plotData[flag2=="B" & type =="variability", value := value * 0.3]
# plotData[flag2=="B" & type =="variability", SE := SE * 0.3]

plotData[, type3 := paste(type2, type, sep = " - ")]

plot12 = ggplot(plotData[Sim_X == "X12",], aes(x=flag2, y=value, color = outcome)) +
  facet_grid2(type2~ type,scales = "free",independent = "y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="darkred", linetype="dotted", 
             linewidth=1.3, aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=value-1.96*SE, ymax=value+1.96*SE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab("Causal estimate") +
  labs(color = "Outcome")
plot12

plot123 = ggplot(plotData[Sim_X == "X123",], aes(x=flag2, y=value, color = outcome)) +
  facet_grid2(type2~ type,scales = "free",independent = "y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="darkred", linetype="dotted", 
             linewidth=1.3, aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=value-1.96*SE, ymax=value+1.96*SE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab("Causal estimate") +
  labs(color = "Outcome")
plot123

plot13 = ggplot(plotData[Sim_X == "X13",], aes(x=flag2, y=value, color = outcome)) +
  facet_grid2(type2~ type,scales = "free",independent = "y") +
  geom_hline(yintercept = 0,color="grey") +
  geom_hline(data = data_hlines, col="darkred", linetype="dotted", 
             linewidth=1.3, aes(yintercept = mylines)) +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=value-1.96*SE, ymax=value+1.96*SE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) + 
  #scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Scenario") + ylab("Causal estimate") +
  labs(color = "Outcome")
plot13

filename = paste0(outdir_results,"/RawEstimate_X12.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot12)
dev.off()

filename = paste0(outdir_results,"/RawEstimate_X13.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot13)
dev.off()

filename = paste0(outdir_results,"/RawEstimate_X123.png")
png(filename = filename,width = 2800, height = 1600, res=200)
print(plot123)
dev.off()

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
