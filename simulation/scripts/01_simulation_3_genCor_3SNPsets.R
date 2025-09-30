#' ---
#' title: "Simulation Study"
#' subtitle: "tba"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: github_document
#' ---
#'
#' # Introduction ####
#' 
#' # Initialize ####
#' ***
#' ## Load R packages
{
  rm(list = ls())
  time0<-Sys.time()
  set.seed(2023)
  
  suppressPackageStartupMessages(library(data.table))
  setDTthreads(1)
  suppressPackageStartupMessages(library(foreach))
  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(doRNG))
  suppressPackageStartupMessages(library(MASS))
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(gamlss))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(MendelianRandomization))
  suppressPackageStartupMessages(library(readxl))
  
}

#' ## Helper function ####
{
  source("../helperfunctions/SimulateSNP.R")
  source("../helperfunctions/HWETest.R")
  source("../helperfunctions/GetAssociation_simulation.R")
  source("../helperfunctions/MVMR_jp_simulation.R")
}

#' # Parameter Settings ####
#' ***
#' ## Fixed parameters ####
fixed_parameters = fread("../temp/parameters_fixed.txt")

{
  n_times_default = fixed_parameters[parameter == "n_times_default",value]
  SNPs_NR = fixed_parameters[parameter == "SNPs_NR",value]
  SNPs_EAF = fixed_parameters[parameter == "SNPs_EAF",value]
  SNPs_mean = fixed_parameters[parameter %in% paste0("SNPs_mean_",c("M","S","V")),value]
  SNPs_mean = SNPs_mean*10
  SNPs_mean[2] = SNPs_mean[2]*10
  SNPs_var = fixed_parameters[parameter %in% paste0("SNPs_var_",c("M","S","V")),value]
  SNPs_var = SNPs_var*10
  SNPs_var[2] = SNPs_var[2]*10
  SNPs_centering = as.logical(fixed_parameters[parameter == "SNPs_centering",value])
  SNPs_Correct_ASs = as.logical(fixed_parameters[parameter == "SNPs_correct_AS",value])
  
  X_mean_R = fixed_parameters[parameter %in% paste0("X_mean_R",c("I","S")),value]
  X_var_R = fixed_parameters[parameter %in% paste0("X_var_R",c("I","S")),value]
  X_covar_R = fixed_parameters[parameter == "X_covar_R",value]
  X_mean_error = fixed_parameters[parameter == "X_mean_error",value]
  X_var_error = fixed_parameters[parameter == "X_var_error",value]
  X_mean_age = fixed_parameters[parameter == "X_mean_age",value]
  X_mean_age2 = fixed_parameters[parameter == "X_mean_age2",value]
  X_mean_int = fixed_parameters[parameter == "X_mean_int",value]
  
  Y_mean = fixed_parameters[parameter == "Y_mean",value]
  Y_var = fixed_parameters[parameter == "Y_var",value]
  Y_theta = fixed_parameters[parameter %in% paste0("Y_theta_",c("M","S","V")),value]
  
  MR_filterBadSNPs = as.logical(fixed_parameters[parameter == "MR_filt",value])
  MR_filterBadSNPs_treshold = fixed_parameters[parameter == "MR_pval_IVs",value]

} 

#' ## Variable parameters ####
variable_parameters = fread("../temp/parameters_variable.txt")
variable_parameters[,NR := 1:dim(variable_parameters)[1]]
variable_parameters = variable_parameters[3]

{
  n_samples = variable_parameters$n_samples
  n_times_setting = variable_parameters$n_times
  n_sim = variable_parameters$n_sim
  gc = c(variable_parameters$gc_MS,variable_parameters$gc_MV,variable_parameters$gc_SV)
  SNPs_factor = c(variable_parameters$SNPs_factor_M, variable_parameters$SNPs_factor_S, variable_parameters$SNPs_factor_V)
  X_age_bl = variable_parameters$X_meanAge_bl
  X_scanSteps = variable_parameters$X_scanSteps
  
  GAMLSS_model = variable_parameters$AssocModelX
  MVMR_model = variable_parameters$MVMRmethod 
  scenario_name = variable_parameters$scenario
  outfiles_prefix = paste(variable_parameters$NR,variable_parameters$scenario,sep="_")
}

#' ## Additional parameters (for debugging)
save_data_perSim = F
do_plotting = F

#' ## Create result directory ####
{
  outfiles_dir = "../result/"
  if(dir.exists(outfiles_dir)==F){
    dir.create(outfiles_dir)
    message("Created results folder: ",outfiles_dir)
  }else{
    message("Using pre-existing results folder: ",outfiles_dir)
  }
  if(dir.exists(paste0(outfiles_dir,outfiles_prefix))==F){
    dir.create(paste0(outfiles_dir,outfiles_prefix))
    message("Created results folder for this scenario: ",paste0(outfiles_dir,outfiles_prefix))
  }else{
    message("Using pre-existing scenario folder: ",paste0(outfiles_dir,outfiles_prefix))
  }
  
}

#' # Simulation ####
#' ***
#' I use a foreach loop to allow for parallelization - the number of cores to be used will be given in the slurm file!
#' 
#' I include a counter to get some output files - not for all simulation but just 10. This should be enough to give me a feeling of want goes wrong (if anything). 
#' 
counter = seq(1,n_sim,n_sim/10)
registerDoParallel(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

SimTab = foreach(s = 1:n_sim)%dorng%{
  #SimTab = foreach(s = 1:n_sim)%do%{
  #s=1
  #message("Working on simulation ",s)
  
  #' ## Step 0: create directory to store data
  #' ***
  {
    outdir_sim = paste0(outfiles_dir,outfiles_prefix,"/Simulation_",s,"/")
    if(save_data_perSim == T | s %in% counter){
      if(dir.exists(outdir_sim)==F){
        dir.create(outdir_sim)
        message("Created results folder ",outdir_sim, " for simulation ",s)
      }else{
        message("Using pre-existing results folder ",outdir_sim, " for simulation ",s)
      }
    }
  }
  
  #' ## Step 1: Simulate genotypes G ####
  #' ***
  {
    #' Get genotypes 
    G = matrix(data = NA, nrow = n_samples, ncol=SNPs_NR)
    
    HWETab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      
      # simulate SNP
      SNP = SimulateSNP(n_samples = n_samples,
                        p = SNPs_EAF,
                        SNP_centered = SNPs_centering)
      G[,j] = SNP
      
      table1 = table(SNP)
      names(table1)[names(table1)==0] = "AA"
      names(table1)[names(table1)==1] = "AB"
      names(table1)[names(table1)==2] = "BB"
      
      AA = table1[names(table1)=="AA"]
      AB = table1[names(table1)=="AB"]
      BB = table1[names(table1)=="BB"]
      
      if(length(AA)==0) AA = 0
      if(length(AB)==0) AB = 0
      if(length(BB)==0) BB = 0
      
      # test HWE
      test = HWETest(AA = AA, 
                     AB = AB,
                     BB = BB)
      test[,NR := j]
      test
    }
    HWETab = rbindlist(HWETab)
    
    #' Check for HWE violation
    HWETab[Pvalue <= 0.05/SNPs_NR,]
    HWETab[Pvalue <= 0.05,]
    if(save_data_perSim == T | s %in% counter) save(HWETab, file = paste0(outdir_sim, "/01_HWE.RData"))
    
    #' Check for LD violation
    CorTab = cor(G)^2
    filt = CorTab>= 0.1
    stopifnot(sum(filt) == SNPs_NR)
    stopifnot(sum(diag(CorTab)) == SNPs_NR)
    rownames(G) = paste0("ID",c(1:n_samples))
    if(save_data_perSim == T | s %in% counter) save(G, file = paste0(outdir_sim, "/01_GenotypeMatrix.RData"))
    
  }
  
  #' ## Step 2: Get allele score AS ####
  #' ***
  {
    #' Get SNP effect 
    CoVarMatrix_SNPs = diag(x=SNPs_var,nrow = 3)
    CoVarMatrix_SNPs[1,2] = gc[1]
    CoVarMatrix_SNPs[2,1] = gc[1]
    CoVarMatrix_SNPs[1,3] = gc[2]
    CoVarMatrix_SNPs[3,1] = gc[2]
    CoVarMatrix_SNPs[2,3] = gc[3]
    CoVarMatrix_SNPs[3,2] = gc[3]
    
    REff_SNPs = MASS::mvrnorm(SNPs_NR, 
                              mu=SNPs_mean, 
                              Sigma=(CoVarMatrix_SNPs))
    REff_SNPs[,1] = REff_SNPs[,1]/10
    REff_SNPs[,2] = REff_SNPs[,2]/100
    REff_SNPs[,3] = REff_SNPs[,3]/10
    
    #' Get filters for the different gc scenarios (correlation - shared SNPs; no correlation - distinct SNPs)
    if(gc[1]==-0.9 & gc[2]==0){
      filtM = c(1:20)
      filtS = c(1:20)
      filtV = c(21:30)
    }else if(gc[1]==0 & gc[2]==0){
      filtM = c(1:10)
      filtS = c(11:20)
      filtV = c(21:30)
    }else if(gc[1]==-0.9 & gc[2]==0.5){
      filtM = c(1:30)
      filtS = c(1:30)
      filtV = c(1:30)
    }

    #' Get scores
    AS1 = G[,filtM] %*% REff_SNPs[filtM,1]
    if(SNPs_Correct_ASs==T) AS1 = ((AS1 - mean(AS1))/sd(AS1))*SNPs_factor[1]
    AS2 = exp(0.5*G[,filtS] %*% REff_SNPs[filtS,2])
    if(SNPs_Correct_ASs==T) AS2 = AS2*SNPs_factor[2]
    AS3 = G[,filtV] %*% REff_SNPs[filtV,3]
    if(SNPs_Correct_ASs==T) AS3 = abs(AS3)*SNPs_factor[3]
    
    #' Plot scores
    if((do_plotting == T & save_data_perSim==T) | s %in% counter){
      plotData = data.table(AS = c(AS1,AS2,AS3),
                            type = c(rep("AS1",n_samples),rep("AS2",n_samples),rep("AS3",n_samples)))
      plotData_mean = plotData[,median(AS),by=type]
      
      ggp0 = ggplot(plotData, aes(x=AS)) +
        facet_wrap(~type,scales = "free") +
        geom_histogram(aes(y=after_stat(density)), colour="black", fill="white",bins=20)+
        geom_density(alpha=.2, fill="#FF6666") +
        geom_vline(data = plotData_mean, aes(xintercept=V1),
                   color="blue", linetype="dashed")
      
      png(filename = paste0(outdir_sim, "02_AlleleScores.png"),
          width = 1400, height = 700, res=125)
      print(ggp0)
      dev.off()
      
    }
    
  }
  
  #' ## Step 3: Simulate exposure X ####
  #' ***
  {
    #' Get random effects
    CoVarMatrix = diag(x=X_var_R,nrow = 2)
    REff = MASS::mvrnorm(n_samples, 
                         mu=X_mean_R, 
                         Sigma=(CoVarMatrix)^2)
    
    #' Get exposure
    myTabX_long= data.table(ID = rep(1:n_samples, each=n_times_default),
                            scan = rep(1:n_times_default, times=n_samples),
                            b0i = rep(REff[,1], each=n_times_default),
                            b1i = rep(REff[,2], each=n_times_default),
                            beta01 = rep(AS1, each=n_times_default),
                            beta11 = rep(AS2, each=n_times_default),
                            tau01 = rep(AS3, each =n_times_default))
    myTabX_long[,ID := paste0("ID",ID)]
    
    for(i in 1:n_times_default){
      myTabX_long[scan==i,age := rnorm(n=n_samples,mean=X_age_bl + X_scanSteps*i,sd=0.5)]
    }
    myTabX_long[,age2 := age^2]
    
    for(i in 1:n_samples){
      #i=1
      myTau01 = AS3[i]
      myTabX_long[ID==paste0("ID",i),epsilon := rnorm(n=n_times_default,mean=0,sd=myTau01)]
      myTabX_long[ID==paste0("ID",i),epsilon2 := rnorm(n=n_times_default,mean=X_mean_error,sd=mean(AS3))]
    }
    
    Xsim_MS = with(myTabX_long, (X_mean_int + b0i + beta01 + 
                                   (X_mean_age + b1i + beta11)*age + 
                                   X_mean_age2*age2 + epsilon2))
    Xsim_MSV = with(myTabX_long, (X_mean_int + b0i + beta01 + 
                                    (X_mean_age + b1i + beta11)*age + 
                                    X_mean_age2*age2 + epsilon))
    Xsim_MV = with(myTabX_long, ((X_mean_int + b0i + beta01 + 
                                    (X_mean_age + b1i + mean(AS2))*age + 
                                    X_mean_age2*age2 + epsilon)) )
    
    myTabX_long[,X_MS := Xsim_MS]
    myTabX_long[,X_MSV := Xsim_MSV]
    myTabX_long[,X_MV := Xsim_MV]
    
    #' Get outcome template
    myTabY_temp = copy(myTabX_long)
    myTabY_temp = myTabY_temp[scan==1,]
    myNames = c("ID","beta01","beta11","tau01")
    colsOut<-setdiff(colnames(myTabY_temp),myNames)
    myTabY_temp[,get("colsOut"):=NULL]
    
    #' Get filters for the different sample scenarios 
    if(scenario_name == "sample_POPS"){
      myTabX_long = myTabX_long[scan %in% c(5,7,9),]
    }else if(scenario_name == "sample_UKB"){
      currentNR = dim(myTabX_long)[1]
      aimedNR = n_samples*n_times_setting
      randomMissingPoints = sample(currentNR,aimedNR,replace = F)
      myTabX_long = myTabX_long[randomMissingPoints,]
      n_samples = length(unique(myTabX_long$ID))
      filt = is.element(rownames(G),unique(myTabX_long$ID))
      G = G[filt,]
      filt = is.element(myTabY_temp$ID,myTabX_long$ID)
      myTabY_temp = myTabY_temp[filt,]
      # dummy = myTabX_long[,.N,by=ID]
      # hist(dummy$N)
      
    }
    
    #' Check AS association (should roughly be 1) - if we know everything
    # summary(gamlss(X_MS ~ b0i + beta01 + age + b1i:age + beta11:age + age2 + random(x = as.factor(ID)),
    #                 sigma.formula = ~ tau01 + age ,
    #                 data = myTabX_long, family = "NO"))
    # summary(gamlss(X_MSV ~ b0i + beta01 + age + b1i:age + beta11:age + age2 + random(x = as.factor(ID)),
    #                sigma.formula = ~ tau01 + age ,
    #                data = myTabX_long, family = "NO"))
    # summary(gamlss(X_MV ~ b0i + beta01 + age + b1i:age + beta11:age + age2 + random(x = as.factor(ID)),
    #                sigma.formula = ~ tau01 + age ,
    #                data = myTabX_long, family = "NO"))
    
    #' Check AS association (should roughly be 1) - if we don't know the random effect and only include RI
    # summary(gamlss(X_MS ~ beta01 + age + beta11:age + age2 + random(x = as.factor(ID)),
    #                sigma.formula = ~ tau01 + age ,
    #                data = myTabX_long, family = "NO"))
    # summary(gamlss(X_MSV ~ beta01 + age + beta11:age + age2 + random(x = as.factor(ID)),
    #                sigma.formula = ~ tau01 + age ,
    #                data = myTabX_long, family = "NO"))
    # summary(gamlss(X_MV ~ beta01 + age + beta11:age + age2 + random(x = as.factor(ID)),
    #                sigma.formula = ~ tau01 + age ,
    #                data = myTabX_long, family = "NO"))
    
    
    #' Get some plots
    if((do_plotting == T & save_data_perSim==T) | s %in% counter){
      plotData = melt(myTabX_long,id.vars=names(myTabX_long)[c(1,2,5:9)],
                      measure.vars=c("X_MSV","X_MS","X_MV"),
                      variable.name="value_type",value.name="value")
      
      ggp1  = ggplot(plotData, aes(x=age, y=value, col=beta01, fill = beta01, group=ID)) +
        facet_wrap(~value_type,scales = "fixed")+
        geom_line(show.legend = TRUE,aes(alpha=0.01)) + 
        geom_point(shape = 21,colour = "black")+
        labs(x="Simulated Age", y="Simulated Exposure", 
             color="AS mean") +
        theme(legend.position = "none") + theme_classic()+ 
        guides(alpha="none",fill="none",color="none")
      
      png(filename = paste0(outdir_sim, "03_TimeVsExposures_Age.png"),
          width = 1400, height = 700, res=125)
      print(ggp1)
      dev.off()
      
    }
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTabX_long, 
           file = paste0(outdir_sim, "03_ExposureData.RData"))
    }
    
  }
  
  #' ## Step 4: Get SNP-X association ####
  #' ***
  {
    #' Get SNP association for the three different exposures, depending on GAMLSS model
    myTab_SNPAssocs_X_MSV = GetAssociation(data = myTabX_long, 
                                           method = GAMLSS_model, 
                                           genotypes = G, 
                                           dep_var_name = "X_MSV",
                                           time_var_name = "age",
                                           growth_var = "quadradic", 
                                           getIA = T)
    
    myTab_SNPAssocs_X_MS = GetAssociation(data = myTabX_long, 
                                          method = GAMLSS_model, 
                                          genotypes = G, 
                                          dep_var_name = "X_MS",
                                          time_var_name = "age",
                                          growth_var = "quadradic", 
                                          getIA = T)
    
    myTab_SNPAssocs_X_MV = GetAssociation(data = myTabX_long, 
                                          method = GAMLSS_model, 
                                          genotypes = G, 
                                          dep_var_name = "X_MV",
                                          time_var_name = "age",
                                          growth_var = "quadradic", 
                                          getIA = T)
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTab_SNPAssocs_X_MSV,myTab_SNPAssocs_X_MS,myTab_SNPAssocs_X_MV,
           file = paste0(outdir_sim, "04_AssociationData.RData"))
    }
    
    #' Check 1: are the SNPs distinct or not? 
    myTab_SNPAssocs_X_MSV[,table(exposure, pval<5e-8)]
    myTab_SNPAssocs_X_MS[,table(exposure, pval<5e-8)]
    myTab_SNPAssocs_X_MV[,table(exposure, pval<5e-8)]
    
    #' Check 2: are the t-values/betas correlated or not?
    if(GAMLSS_model == "gamlssNoIA"){
      c4 = cor.test(myTab_SNPAssocs_X_MSV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MSV[exposure=="var",beta])
      c5 = cor.test(myTab_SNPAssocs_X_MS[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MS[exposure=="var",beta])
      c6 = cor.test(myTab_SNPAssocs_X_MV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MV[exposure=="var",beta])
      c13 = cor.test(AS1,AS3)
      
      myCorTable = data.table(data1 = c(rep("GAMLSS_mean",3),"AS_mean"),
                              data2 = c(rep("GAMLSS_var",3),"AS_var"),
                              exposure = c("MSV","MS","MV","AS"),
                              cor = c(c4$estimate,c5$estimate,c6$estimate,
                                      c13$estimate),
                              pval = c(c4$p.value,c5$p.value,c6$p.value,
                                       c13$p.value))
    }else if(GAMLSS_model == "gamlssNoVar"){
      c1 = cor.test(myTab_SNPAssocs_X_MSV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MSV[exposure=="slope",beta])
      c2 = cor.test(myTab_SNPAssocs_X_MS[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MS[exposure=="slope",beta])
      c3 = cor.test(myTab_SNPAssocs_X_MV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MV[exposure=="slope",beta])
      c12 = cor.test(AS1,AS2)
      
      myCorTable = data.table(data1 = c(rep("GAMLSS_mean",3),"AS_mean"),
                              data2 = c(rep("GAMLSS_slope",3),"AS_slope"),
                              exposure = c("MSV","MS","MV","AS"),
                              cor = c(c1$estimate,c2$estimate,c3$estimate,
                                      c12$estimate),
                              pval = c(c1$p.value,c2$p.value,c3$p.value,
                                       c12$p.value))
      
    }else{
      c1 = cor.test(myTab_SNPAssocs_X_MSV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MSV[exposure=="slope",beta])
      c2 = cor.test(myTab_SNPAssocs_X_MS[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MS[exposure=="slope",beta])
      c3 = cor.test(myTab_SNPAssocs_X_MV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MV[exposure=="slope",beta])
      c4 = cor.test(myTab_SNPAssocs_X_MSV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MSV[exposure=="var",beta])
      c5 = cor.test(myTab_SNPAssocs_X_MS[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MS[exposure=="var",beta])
      c6 = cor.test(myTab_SNPAssocs_X_MV[exposure=="mean",beta],
                    myTab_SNPAssocs_X_MV[exposure=="var",beta])
      c7 = cor.test(myTab_SNPAssocs_X_MSV[exposure=="slope",beta],
                    myTab_SNPAssocs_X_MSV[exposure=="var",beta])
      c8 = cor.test(myTab_SNPAssocs_X_MS[exposure=="slope",beta],
                    myTab_SNPAssocs_X_MS[exposure=="var",beta])
      c9 = cor.test(myTab_SNPAssocs_X_MV[exposure=="slope",beta],
                    myTab_SNPAssocs_X_MV[exposure=="var",beta])
      c12 = cor.test(AS1,AS2)
      c13 = cor.test(AS1,AS3)
      c23 = cor.test(AS2,AS3)
      
      myCorTable = data.table(data1 = c(rep("GAMLSS_mean",6),rep("GAMLSS_slope",3),
                                        "AS_mean","AS_mean","AS_slope"),
                              data2 = c(rep("GAMLSS_slope",3),rep("GAMLSS_var",6),
                                        "AS_slope","AS_var","AS_var"),
                              exposure = c(rep(c("MSV","MS","MV"),3),rep("AS",3)),
                              cor = c(c1$estimate,c2$estimate,c3$estimate,
                                      c4$estimate,c5$estimate,c6$estimate,
                                      c7$estimate,c8$estimate,c9$estimate,
                                      c12$estimate,c13$estimate,c23$estimate),
                              pval = c(c1$p.value,c2$p.value,c3$p.value,
                                       c4$p.value,c5$p.value,c6$p.value,
                                       c7$p.value,c8$p.value,c9$p.value,
                                       c12$p.value,c13$p.value,c23$p.value))
    }
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myCorTable, file = paste0(outdir_sim, "04_CorrelationSNPEffects.RData"))
    }
    
  }
  
  #' ## Step 5: Simulate outcome Y ####
  #' ***
  {
    myTabY = copy(myTabY_temp)
    myMean = myTabX_long[scan==max(scan),mean(age)]
    myMean = round(myMean,0)+10
    myTabY[,age := rnorm(n=n_samples,mean=myMean,sd=1)]
    
    myTabY[,Y1 := rnorm(n=n_samples,mean = Y_mean,sd = Y_var)]
    myTabY[,Y2 := Y_theta[1] * beta01 + rnorm(n=n_samples,mean = Y_mean,sd = sd(beta01))]
    myTabY[,Y3 := Y_theta[2] * beta11 * age + 
             rnorm(n=n_samples,mean = Y_mean,sd = sd(beta11))]
    myTabY[,Y4 := Y_theta[3] * log(tau01) + rnorm(n=n_samples,mean = Y_mean,sd = sd(tau01))]
    myTabY[,Y5 := Y_theta[1] * beta01 + Y_theta[2] * beta11 * age + 
             rnorm(n=n_samples,mean = Y_mean,sd = sd(beta11)*sd(beta01))]
    myTabY[,Y6 := Y_theta[1] * beta01 + Y_theta[3] * log(tau01) + 
             rnorm(n=n_samples,mean = Y_mean,sd = sd(tau01)*sd(beta01))]
    myTabY[,Y7 := Y_theta[3] * log(tau01) + Y_theta[2] * beta11 * age + 
             rnorm(n=n_samples,mean = Y_mean,sd = sd(beta11)*sd(tau01))]
    myTabY[,Y8 := Y_theta[1] * beta01 + Y_theta[2] * beta11 * age + Y_theta[3] * log(tau01) +
             rnorm(n=n_samples,mean = Y_mean,sd = sd(beta11)*sd(beta01)*sd(tau01))]
    
    #' Check AS association
    # summary(lm(Y1 ~ beta01 + beta11 + tau01, data=myTabY))
    # summary(lm(Y2 ~ beta01 + beta11 + tau01, data=myTabY))
    # summary(lm(Y3 ~ beta01 + beta11 + tau01, data=myTabY))
    # summary(lm(Y4 ~ beta01 + beta11 + tau01, data=myTabY))
    # summary(lm(Y5 ~ beta01 + beta11 + tau01, data=myTabY))
    # summary(lm(Y6 ~ beta01 + beta11 + tau01, data=myTabY))
    # summary(lm(Y7 ~ beta01 + beta11 + tau01, data=myTabY))
    # summary(lm(Y8 ~ beta01 + beta11 + tau01, data=myTabY))
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTabY,
           file = paste0(outdir_sim, "05_OutcomeData.RData"))
    }
  }
  
  #' ## Step 6: Get SNP-Y association ####
  #' ***
  {
    #' Get associations
    myVars = names(myTabY)[grepl("Y",names(myTabY))]
    
    myTab_SNPAssocs = foreach(i = 1:length(myVars))%do%{
      #i=1
      myMethod = "linReg" 
      if(grepl("_bin",myVars[i])==T) myMethod = "glm"
      
      myTab_GY1 = GetAssociation(data = myTabY, 
                                 method = myMethod, 
                                 genotypes = G, 
                                 dep_var_name = myVars[i],
                                 time_var_name = "age",
                                 growth_var = set_growth,
                                 getIA = F)
      
      setnames(myTab_GY1,"exposure","outcome")
      myTab_GY1[,outcome := myVars[i]]
      myTab_GY1
    }
    myTab_SNPAssocs_Y = rbindlist(myTab_SNPAssocs)
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTab_SNPAssocs_Y,
           file = paste0(outdir_sim, "06_AssociationData.RData"))
    }
  }
  
  #' ## Step 7: Run MVMR ####
  #' ***
  {
    #' Do age correction on GAMLSS results
    mean_age_outcome = myTabY[,mean(age)]
    myTab_SNPAssocs_X_MSV[exposure=="slope", beta := beta*mean_age_outcome]
    myTab_SNPAssocs_X_MSV[exposure=="slope", SE := SE*mean_age_outcome]
    myTab_SNPAssocs_X_MS[exposure=="slope", beta := beta*mean_age_outcome]
    myTab_SNPAssocs_X_MS[exposure=="slope", SE := SE*mean_age_outcome]
    myTab_SNPAssocs_X_MV[exposure=="slope", beta := beta*mean_age_outcome]
    myTab_SNPAssocs_X_MV[exposure=="slope", SE := SE*mean_age_outcome]
    
    #' Run analysis
    MR_analysis_MSV = MVMR_jp_sim(data_GX_long = myTab_SNPAssocs_X_MSV,
                                  data_GY_long = myTab_SNPAssocs_Y, 
                                  filterBadSNPs = MR_filterBadSNPs,
                                  filterBadSNPs_threshold = MR_filterBadSNPs_treshold,
                                  sampleSize_GX = dim(myTabX_long)[1],
                                  sampleSize_GY = dim(myTabY)[1], 
                                  MVMR_method=MVMR_model)
    
    MR_analysis_MS = MVMR_jp_sim(data_GX_long = myTab_SNPAssocs_X_MS,
                                 data_GY_long = myTab_SNPAssocs_Y, 
                                 filterBadSNPs = MR_filterBadSNPs,
                                 filterBadSNPs_threshold = MR_filterBadSNPs_treshold,
                                 sampleSize_GX = dim(myTabX_long)[1],
                                 sampleSize_GY = dim(myTabY)[1], 
                                 MVMR_method=MVMR_model)
    
    MR_analysis_MV = MVMR_jp_sim(data_GX_long = myTab_SNPAssocs_X_MV,
                                 data_GY_long = myTab_SNPAssocs_Y, 
                                 filterBadSNPs = MR_filterBadSNPs,
                                 filterBadSNPs_threshold = MR_filterBadSNPs_treshold,
                                 sampleSize_GX = dim(myTabX_long)[1],
                                 sampleSize_GY = dim(myTabY)[1], 
                                 MVMR_method=MVMR_model)
    
    MR_analysis_MSV[,comment2 := "MSV"]
    MR_analysis_MS[,comment2 := "MS"]
    MR_analysis_MV[,comment2 := "MV"]
    
    MRTab = rbind(MR_analysis_MSV,MR_analysis_MS,MR_analysis_MV)
    
    #' Get some plots
    if((do_plotting == T & save_data_perSim==T) | s %in% counter){
      if(scenario_name == "GAMLSS_noSlope"){
        data_hlines = data.frame(exposure = unique(MRTab$exposure),
                                 mylines = Y_theta[c(1,3)])
      }else if(scenario_name == "GAMLSS_noVar"){
        data_hlines = data.frame(exposure = unique(MRTab$exposure),
                                 mylines = Y_theta[c(1,2)])
      }else{
        data_hlines = data.frame(exposure = unique(MRTab$exposure),
                                 mylines = Y_theta)
      }
      MRTab[,lowerCI95 := beta-1.96*SE]
      MRTab[,upperCI95 := beta+1.96*SE]
      
      ggp3 = ggplot(MRTab, aes(x=outcome, y=beta, col=comment2,
                                               ymin=lowerCI95, ymax=upperCI95)) + 
        facet_wrap(~exposure,scales = "free")+
        geom_hline(data = data_hlines, aes(yintercept = mylines),linetype="dashed", show.legend = FALSE) +
        geom_hline(yintercept = 0,linetype="solid")+
        geom_pointrange(position = position_dodge2(width = 0.6),size=0.6,linewidth=1) +
        xlab("Outcomes") +
        ylab("MVMR estimate") + 
        labs(color = "Exposures")
      
      png(filename = paste0(outdir_sim, "07_MREstimates.png"), 
          width = 1200, height = 600, res=125)
      print(ggp3)
      dev.off()
      
      ggp4 = ggplot(MRTab[outcome=="Y8"], aes(x=comment2, y=condFStats, col=comment2,
                               ymin=lowerCI95, ymax=upperCI95)) + 
        facet_wrap(~exposure,scales = "fixed")+
        geom_hline(yintercept = 10, linetype="dashed")+
        geom_point() +
        xlab("Exposures") +
        ylab("MVMR estimate") + 
        labs(color = "Exposures")
      
      png(filename = paste0(outdir_sim, "08_CondFStats.png"), 
          width = 1200, height = 600, res=125)
      print(ggp4)
      dev.off()
      
    }
    
  }
  
  MRTab[,n_sim := s]
  MRTab
}
SimTab = rbindlist(SimTab,fill=T)
save(SimTab, file = paste0(outfiles_dir,outfiles_prefix, "_SimAll.RData"))
save(variable_parameters, file = paste0(outfiles_dir,outfiles_prefix, "_VarParameters.RData"))
save(fixed_parameters, file = paste0(outfiles_dir,outfiles_prefix, "_FixParameters.RData"))

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
