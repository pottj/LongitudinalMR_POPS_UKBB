#' ---
#' title: "Supplemental simulation figure: Power heatmaps"
#' subtitle: "Supplemental Figures"
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
  message("Created results folder ",outdir_results, " for main figures ")
}else{
  message("Using pre-existing results folder ",outdir_results, " for main figures ")
}

#' # Get data ####
#' ***
load("../results/SupTabs.RData")

myTab = copy(tab2)
myTab = myTab[Sim_Y %in% c("CM2")]

todoList = data.table(nr = 1:9,
                      exposure = rep(c("X123","X12","X13"),each=3),
                      exposure_type = rep(c("mean","slope","var"),3))
todoList

#' # Get plots ####
#' ***
colNames_ordered = c("A","B","A_sampleSize","B_sampleSize","A_noSlope","B_noSlope","A_noVar","A_distSNPs","A_noCor",
                     "A_posCor","A_GxE","B_GxE","A_RIsigma","A_binary20","A_binary50")

plotting = foreach(i=1:9)%do%{
  #i=2
  myRow = todoList[i,]
  plotData = copy(myTab)
  plotData = plotData[Sim_X == myRow$exposure]
  plotData = plotData[!is.na(get(paste0("power_",myRow$exposure_type))),]
  dumTab = dcast(plotData, outcome ~ check, value.var=paste0("power_",myRow$exposure_type))
  dumTab[,outcome := NULL]
  names(dumTab) = gsub("main_","",names(dumTab))
  names(dumTab) = gsub("sens_","",names(dumTab))
  colNames_ordered2 = colNames_ordered[colNames_ordered %in% names(dumTab)]
  setcolorder(dumTab,colNames_ordered2)
  
  dumMat = as.matrix(dumTab)
  rownames(dumMat) = paste0("Y",1:8)
  
  filename = paste0(outdir_results,"/Power_",myRow$exposure,"_",myRow$exposure_type,".png")
  png(filename = filename,width = 2600, height = 1400, res=200)
  corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
           col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
           addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
  dev.off()
  
  # I want to export the matrix to later modify in power point or spread sheet
  dumMat
}

plot1 = as.data.table(plotting[[1]])
plot2 = as.data.table(plotting[[2]])
plot3 = as.data.table(plotting[[3]])
plot4 = as.data.table(plotting[[4]])
plot5 = as.data.table(plotting[[5]])
plot6 = as.data.table(plotting[[6]])
plot7 = as.data.table(plotting[[7]])
plot8 = as.data.table(plotting[[8]])
plot9 = as.data.table(plotting[[9]])

plot_X123 = rbind(plot1,plot2,plot3,fill=T)
plot_X12 = rbind(plot4,plot5,plot6,fill=T)
plot_X13 = rbind(plot7,plot8,plot9,fill=T)

WriteXLS(c("plot_X123","plot_X12","plot_X13","plot1","plot2","plot3","plot4","plot5","plot6","plot7","plot8","plot9"),
         ExcelFileName=paste0("../results/SupFigs_Power.xlsx"),
         SheetNames=c("X123","X12","X13",paste(todoList$exposure,todoList$exposure_type, sep="_")),
         AutoFilter=F,
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
