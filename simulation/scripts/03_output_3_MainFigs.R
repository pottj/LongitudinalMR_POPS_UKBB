#' ---
#' title: "Main Figures of Simulation"
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
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggplot2))

#' # Load data ####
#' ***
load(file="../result/_tables/SupTabs.RData")

#' # Power plots ####
#' ***
#' 
outdir_results = "../result/_MainFigures/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
}
dumTab = copy(tab2)
dumTab = dumTab[Scenario_NR == "0",]

dumTab_mean <- dcast(dumTab, outcome ~ exposure, value.var="power_mean")
dumTab_slope <- dcast(dumTab, outcome ~ exposure, value.var="power_slope")
dumTab_var <- dcast(dumTab, outcome ~ exposure, value.var="power_var")

dumTab2 = cbind(dumTab_mean,dumTab_slope[,-1],dumTab_var[,-1])
dumMat = as.matrix(dumTab2[,-1])
dumMat = dumMat[,order(colnames(dumMat))]
dumMat = dumMat[,c(4:6,1:3,7:9)]
colnames(dumMat)[c(1,4,7)] = paste0("mean_",colnames(dumMat)[c(1,4,7)])
colnames(dumMat)[c(2,5,8)] = paste0("slope_",colnames(dumMat)[c(2,5,8)])
colnames(dumMat)[c(3,6,9)] = paste0("var_",colnames(dumMat)[c(3,6,9)])
rownames(dumMat) = paste("Y_",c(1:8))
filename = paste0(outdir_results,"/Power.png")
png(filename = filename,width = 1800, height = 1400, res=200)
corrplot(dumMat, is.corr = FALSE,col.lim = c(0, 1),#col = COL1('Reds'), 
         col= colorRampPalette(c("#FFF5F0","#FB6A4A","#67000D"))(10),
         addCoef.col = 'grey50',method = 'color',tl.col = "black",tl.srt = 45)
dev.off()

#' # Bias plots ####
#' ***
#' 
dumTab1 = copy(tab2)
dumTab1[,type := "mean"]
dumTab1[,bias := bias_mean]
dumTab1[,bias_SE := bias_SE_mean]

dumTab2 = copy(tab2)
dumTab2[,type := "slope"]
dumTab2[,bias := bias_slope]
dumTab2[,bias_SE := bias_SE_slope]

dumTab3 = copy(tab2)
dumTab3[,type := "var"]
dumTab3[,bias := bias_var]
dumTab3[,bias_SE := bias_SE_var]

dumTab4 = rbind(dumTab1,dumTab2,dumTab3)  
dumTab4 = dumTab4[!is.na(bias),]
dumTab4 = dumTab4[,c(1:4,41:43)]
dumTab4 = dumTab4[Scenario_NR=="0"]

dumTab4[,outcome := gsub("Y","Y_",outcome)]
dumTab4[exposure == "MSV", exposure := "X1_MSV"]
dumTab4[exposure == "MS", exposure := "X2_MS"]
dumTab4[exposure == "MV", exposure := "X3_MV"]

exposure_name <- c(
  X2_MS = "X^(MS)",
  X1_MSV = "X^(MSV)",
  X3_MV = "X^(MV)"
)

plot5 = ggplot(dumTab4, aes(x=type, y=bias, color = outcome)) +
  facet_wrap(~ exposure,scales = "free_y",labeller = labeller(exposure = exposure_name)) +
  geom_hline(yintercept = 0,color="grey") +
  geom_point(position=position_dodge(0.5),size=3) +
  geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                position=position_dodge(0.5)) +
  theme_bw(base_size = 15) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  xlab("Exposure type") + ylab("Bias") +
  labs(color = "Outcome")
plot5

filename = paste0(outdir_results,"/Bias.png")
png(filename = filename,width = 2000, height = 1200, res=200)
print(plot5)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
