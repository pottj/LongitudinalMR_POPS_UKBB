#' ---
#' title: "Supplemental Figures of Simulation"
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

outdir_results = "../result/_SupplementalFigures/"
if(dir.exists(outdir_results)==F){
  dir.create(outdir_results)
}

#' # Prep data ####
#' ***
#' get power in long format
myTab = copy(tab2)
myTab[exposure == "MSV", exposure := "X1_MSV"]
myTab[exposure == "MS", exposure := "X2_MS"]
myTab[exposure == "MV", exposure := "X3_MV"]

#' # Power heat map ####
#' ***
{
  dumTab1 = copy(myTab)
  dumTab1[,type := "mean"]
  dumTab1[,power := power_mean]
  dumTab2 = copy(myTab)
  dumTab2[,type := "slope"]
  dumTab2[,power := power_slope]
  dumTab3 = copy(myTab)
  dumTab3[,type := "var"]
  dumTab3[,power := power_var]
  
  dumTab4_power = rbind(dumTab1,dumTab2,dumTab3)  
  dumTab4_power = dumTab4_power[!is.na(power),]
  dumTab4_power = dumTab4_power[,c(1:4,41:42)]
  dumTab4_power[,outcome := gsub("Y","Y_",outcome)]
  
  exposure_name <- c(
    X2_MS = "X^(MS)",
    X1_MSV = "X^(MSV)",
    X3_MV = "X^(MV)"
  )
  
  plot0 = ggplot(dumTab4_power, aes(Scenario_NR, outcome)) + 
    facet_grid(type~exposure, labeller = labeller(
      exposure = exposure_name))+ 
    geom_tile(aes(fill = power)) + 
    scale_fill_gradient(low = "#FFF5F0", high = "#67000D") + 
    theme_dark(base_size = 15) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    xlab("Scenario number") + ylab("Outcome") +
    labs(fill = "Power")
  
  filename = paste0(outdir_results,"/Power.png")
  png(filename = filename,width = 2800, height = 1600, res=200)
  print(plot0)
  dev.off()
}

#' # Bias plot ####
#' ***
{
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
  
  dumTab4_bias = rbind(dumTab1,dumTab2,dumTab3)  
  dumTab4_bias = dumTab4_bias[!is.na(bias),]
  dumTab4_bias = dumTab4_bias[,c(1:4,40:43)]
  dumTab4_bias[,outcome := gsub("Y","Y_",outcome)]
  
  plot1 = ggplot(dumTab4_bias[Scenario_NR %in% c("0","1A","1B")], 
                  aes(x=type, y=bias, color = outcome)) +
    facet_grid(Scenario_NR~exposure, labeller = labeller(
      exposure = exposure_name), scales = "free")+ 
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 20) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure type") + ylab("Bias") +
    labs(color = "Outcome")
  plot1
  
  plot2 = ggplot(dumTab4_bias[Scenario_NR %in% c("0","2A","2B")], 
                  aes(x=type, y=bias, color = outcome)) +
    facet_grid(Scenario_NR~exposure, labeller = labeller(
      exposure = exposure_name), scales = "free")+ 
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 20) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure type") + ylab("Bias") +
    labs(color = "Outcome")
  plot2
  
  plot3 = ggplot(dumTab4_bias[Scenario_NR %in% c("0","3A","3B")], 
                  aes(x=type, y=bias, color = outcome)) +
    facet_grid(Scenario_NR~exposure, labeller = labeller(
      exposure = exposure_name), scales = "free")+ 
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 20) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure type") + ylab("Bias") +
    labs(color = "Outcome")
  plot3

  plot4 = ggplot(dumTab4_bias[Scenario_NR %in% c("0","4A","4B")], 
                 aes(x=type, y=bias, color = outcome)) +
    facet_grid(Scenario_NR~exposure, labeller = labeller(
      exposure = exposure_name), scales = "free")+ 
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 20) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure type") + ylab("Bias") +
    labs(color = "Outcome")
  plot4
  
  plot5 = ggplot(dumTab4_bias[Scenario_NR %in% c("0","5A","5B")], 
                  aes(x=type, y=bias, color = outcome)) +
    facet_grid(Scenario_NR~exposure, labeller = labeller(
      exposure = exposure_name), scales = "free")+ 
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=bias-1.96*bias_SE, ymax=bias+1.96*bias_SE), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 20) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Exposure type") + ylab("Bias") +
    labs(color = "Outcome")
  plot5
  
  filename = paste0(outdir_results,"/Bias_Main_vs_SampleSize.png")
  png(filename = filename,width = 2400, height = 1600, res=200)
  print(plot1)
  dev.off()
  
  filename = paste0(outdir_results,"/Bias_Main_vs_GenCor.png")
  png(filename = filename,width = 2400, height = 1600, res=200)
  print(plot2)
  dev.off()
  
  filename = paste0(outdir_results,"/Bias_Main_vs_GAMLSS.png")
  png(filename = filename,width = 2400, height = 1600, res=200)
  print(plot3)
  dev.off()
  
  filename = paste0(outdir_results,"/Bias_Main_vs_MVMR.png")
  png(filename = filename,width = 2400, height = 1600, res=200)
  print(plot4)
  dev.off()
  
  filename = paste0(outdir_results,"/Bias_Main_vs_AgeMisSpec.png")
  png(filename = filename,width = 2400, height = 1600, res=200)
  print(plot5)
  dev.off()
  
}

#' # Conditional F-statistics plot ####
#' ***
{
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
  
  dumTab4_condF = rbind(dumTab1,dumTab2,dumTab3)
  dumTab4_condF = dumTab4_condF[!is.na(condFStats_median),]
  dumTab4_condF = dumTab4_condF[,c(1:4,40:44)]
  dumTab4_condF = dumTab4_condF[outcome=="Y8",]
  
  plot6 = ggplot(dumTab4_condF, 
                  aes(x=Scenario_NR, y=condFStats_median, color = exposure)) +
    facet_grid(type~exposure, labeller = labeller(
      exposure = exposure_name), scales = "free")+ 
    geom_hline(yintercept = 0,color="grey") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 20) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario number") + ylab("Conditional F-Statistic") +
    theme(legend.position = "none")
  plot6
  
  plot7 = ggplot(dumTab4_condF[condFStats_median<250,], 
                  aes(x=Scenario_NR, y=condFStats_median, color = exposure)) +
    facet_grid(type~exposure, labeller = labeller(
      exposure = exposure_name), scales = "free")+ 
    geom_hline(yintercept = 0,color="grey") +
    geom_hline(yintercept = 10,color="black",linetype="dashed") +
    geom_point(position=position_dodge(0.5),size=3) +
    geom_errorbar(aes(ymin=condFStats_1stQ, ymax=condFStats_3rdQ), width=.2,
                  position=position_dodge(0.5)) +
    theme_bw(base_size = 20) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    xlab("Scenario number") + ylab("Conditional F-Statistic") +
    theme(legend.position = "none")
  plot7
  
  filename = paste0(outdir_results,"/CondF_complete.png")
  png(filename = filename,width = 2400, height = 1600, res=200)
  print(plot6)
  dev.off()
  
  filename = paste0(outdir_results,"/CondF_lessThan250.png")
  png(filename = filename,width = 2400, height = 1600, res=200)
  print(plot7)
  dev.off()
  
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
