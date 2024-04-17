#' ---
#' title: "Simulation Study"
#' subtitle: "Sensitivity 1 - shared SNP sets"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output: github_document
#' ---
#'
#' # Introduction ####
#' 
#' **In simulation**
#' 
#' - Time: age
#' - Growth: quadratic
#' - SNP sets: shared
#' 
#' **Regression model**: GAMLSS
#' 
#' - Time: scan
#' - Growth: quadratic
#' - SNP x Time - interaction: TRUE
#' - Sigma dependencies: time and SNP
#' 
#'
#' # Initialize ####
#' ***
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
  
}

#' # Parameter Settings ####
#' ***
#' ## Helper function ####
{
  source("../../helperfunctions/SimulateSNP.R")
  source("../../helperfunctions/HWETest.R")
  source("../../helperfunctions/GetAssociation_simulation_v2.R")
  source("../../helperfunctions/MVMR_jp_simulation.R")
  
}

#' ## General settings ####
#' These parameters should be fix for all simulations!
{
  n_samples = 3000
  n_times = 5
  n_sim = 100
  outfiles_dir = "results/"
  save_data_perSim = F
  do_plotting = F
  
  SNPs_NR = 20
  SNPs_EAF = 0.25
  SNPs_centering = F
  SNPs_Correct_ASs = T
  
  X_beta_mean = c(-0.2,0.02)
  X_beta_sd = c(0.05,0.05)
  X_mean_random = c(0,1)
  X_var_random = c(1,0.25)
  X_covar_random = 0
  X_error_mean = 0
  X_error_sd = 1
  
  Y_alpha = c(0.3,0.3)
  Y_mean_random = 0
  Y_var_random = 1
  
  MR_filterBadSNPs = T
  MR_filterBadSNPs_treshold = 1e-6
  MR_doCorrection = F
  
} 

#' ## Specific settings ####
#' These parameters change between the scenarios!
{
  outfiles_prefix = "Sim_SENS1_shared"
  
  set_time = "age"         # either age or scan

  #SNPs_classes = c(rep("A",SNPs_NR/2),rep("B",SNPs_NR/2))
  SNPs_classes = c(rep("A",SNPs_NR))
  tag = "sens1_shared"
  
  set_growth = "quadradic"     # either linear or quadradic
  
  AssocModelX = "gamlssIA"  # either linMixed or gamlssIA or gamlss
  
}

#' ## Create result directory ####
{
  if(dir.exists(paste0("../",outfiles_dir,"/",outfiles_prefix))==F){
    dir.create(paste0("../",outfiles_dir,"/",outfiles_prefix))
    message("Created results folder ",paste0("../",outfiles_dir,"/",outfiles_prefix))
  }else{
    message("Using pre-existing results folder ",paste0("../",outfiles_dir,"/",outfiles_prefix))
  }
  
}

#' ## Summarize main input ####
message("Fixed parameters:", 
        "\n - number of samples: ",n_samples,
        "\n - number of time points: ",n_times,
        "\n - number of simulations: ",n_sim,
        "\n - number of SNPs: ",SNPs_NR)

message("Scenario setting:", 
        "\n - tag: ",tag,
        "\n - age in growth function: ", set_time,
        "\n - number of SNP sets: ",length(unique(SNPs_classes))," (1=shared; 2=distinct)",
        "\n - growth function: ",set_growth,
        "\n - regression model for exposure: ",AssocModelX)

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
    outdir_sim = paste0("../",outfiles_dir,"/",outfiles_prefix,"/Simulation_",s,"/")
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
    if(save_data_perSim == T | s %in% counter) save(G, file = paste0(outdir_sim, "/01_GenotypeMatrix.RData"))
    
  }
  
  #' ## Step 2: Get allele score AS ####
  #' ***
  {
    #' Get SNP effect 
    if(length(unique(SNPs_classes))==2){
      CoVarMatrix_SNPs = diag(x=X_beta_sd,nrow = 2)
      REff_SNPs = MASS::mvrnorm(SNPs_NR, 
                                mu=X_beta_mean, 
                                Sigma=(CoVarMatrix_SNPs)^2)
      filtA = SNPs_classes == "A"
      filtB = SNPs_classes == "B"  
      AS_factor = c(2,2)  

    }else{
      CoVarMatrix_SNPs = diag(x=X_beta_sd,nrow = 2)
      CoVarMatrix_SNPs[1,2] = 0.05
      CoVarMatrix_SNPs[2,1] = 0.05
      REff_SNPs = MASS::mvrnorm(SNPs_NR, 
                                mu=X_beta_mean, 
                                Sigma=(CoVarMatrix_SNPs)^2)
      filtA = SNPs_classes == "A"
      filtB = SNPs_classes == "A"  
      AS_factor = c(1,3)  
      
    }
    
    #' Get scores
    AS1 = G[,filtA] %*% REff_SNPs[filtA,1]
    if(SNPs_Correct_ASs==T) AS1 = ((AS1 - mean(AS1))/sd(AS1))*AS_factor[1]
    AS2 = exp(0.5*G[,filtB] %*% REff_SNPs[filtB,2])
    if(SNPs_Correct_ASs==T) AS2 = AS2*AS_factor[2]
    
    #' Plot scores
    if((do_plotting == T & save_data_perSim==T) | s %in% counter){
      plotData = data.table(AS = c(AS1,AS2),
                            type = c(rep("AS1",n_samples),rep("AS2",n_samples)))
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
    CoVarMatrix = diag(x=X_var_random,nrow = 2)
    REff = MASS::mvrnorm(n_samples, 
                         mu=X_mean_random, 
                         Sigma=(CoVarMatrix)^2)
    
    #' Get exposure
    myTabX_long= data.table(ID = rep(1:n_samples, each=n_times),
                            scan = rep(1:n_times, times=n_samples),
                            beta00 = rep(REff[,1], each=n_times),
                            beta10 = rep(REff[,2], each=n_times),
                            beta01=rep(AS1, each=n_times),
                            beta11=rep(AS2, each=n_times))
    
    if(set_time=="scan"){
      myTabX_long[,age:= scan-1]
      
    }else{
      myTabX_long[scan==1,age := rnorm(n=n_samples,mean=4,sd=1)]
      myTabX_long[scan==2,age := rnorm(n=n_samples,mean=12,sd=1)]
      myTabX_long[scan==3,age := rnorm(n=n_samples,mean=20,sd=1)]
      myTabX_long[scan==4,age := rnorm(n=n_samples,mean=28,sd=1)]
      myTabX_long[scan==5,age := rnorm(n=n_samples,mean=36,sd=1)]
      
    }
    myTabX_long[,age2 := age^2]
    
    if(set_growth == "linear"){
      Xsim = with(myTabX_long, (beta00 + beta01 + beta10*age + beta11*age) 
                  + rnorm(n=n_samples*n_times, mean=X_error_mean, sd=sqrt(X_error_sd)))
    }else{
      Xsim = with(myTabX_long, (beta00 + beta01 + beta10*age + beta11*age + 0.1*age2) 
                  + rnorm(n=n_samples*n_times, mean=X_error_mean, sd=sqrt(X_error_sd)))
    }
    
    myTabX_long[,A := Xsim]
    # myTabX_long[,Z := scale(A),by = scan]
    # myTabX_long[,C := pnorm(Z),by = scan]
    
    #' Check AS association
    # summary(lm(A ~ beta01 + age + beta11:age, data=myTabX_long, subset = scan==1))
    # summary(lm(A ~ beta01 + age + beta11:age, data=myTabX_long, subset = scan==2))
    # summary(lm(A ~ beta01 + age + beta11:age, data=myTabX_long, subset = scan==3))
    # summary(lm(A ~ beta01 + age + beta11:age, data=myTabX_long, subset = scan==4))
    # summary(lm(A ~ beta01 + age + beta11:age, data=myTabX_long, subset = scan==5))
    # summary(lmer(A ~ beta01 + age + beta11:age + (1|ID), data=myTabX_long))
    
    #' Get some plots
    if((do_plotting == T & save_data_perSim==T) | s %in% counter){
      plotData = melt(myTabX_long,id.vars=names(myTabX_long)[c(1,2,5,6,7)],
                      measure.vars=c("A"),
                      #measure.vars=c("A", "C", "Z"),
                      variable.name="value_type",value.name="value")
      
      ggp1  = ggplot(plotData, aes(x=age, y=value, col=beta11, fill = beta11, group=ID)) +
        facet_wrap(~value_type,scales = "free_y")+
        geom_line(show.legend = TRUE,aes(alpha=0.01)) + 
        geom_point(shape = 21,colour = "black")+
        labs(x="Gestational Week", y="Simulated Exposure", 
             color="PRS") +
        theme(legend.position = "none") + theme_classic()+ 
        guides(alpha="none",fill="none")
      
      ggp2  = ggplot(plotData, aes(x=scan, y=value, col=beta11, fill = beta11, group=ID)) +
        facet_wrap(~value_type,scales = "free_y")+
        geom_line(show.legend = TRUE,aes(alpha=0.01)) + 
        geom_point(shape = 21,colour = "black")+
        labs(x="Scan Numbers", y="Simulated Exposure", 
             color="PRS") +
        theme(legend.position = "none") + theme_classic()+ 
        guides(alpha="none",fill="none")
      
      
      png(filename = paste0(outdir_sim, "03_TimeVsExposures_Age.png"),
          width = 1400, height = 700, res=125)
      print(ggp1)
      dev.off()
      
      png(filename = paste0(outdir_sim, "03_TimeVsExposures_Scan.png"),
          width = 1400, height = 700, res=125)
      print(ggp2)
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
    if(grepl("main",tag)){
      # Main model: time = age, quadratic growth, gamlssIA regression
      myTab_SNPAssocs_A = GetAssociation_v2(data = myTabX_long, 
                                         method = AssocModelX, 
                                         genotypes = G, 
                                         dep_var_name = "A",
                                         time_var_name = "age",
                                         growth_var = set_growth, getIA = T)
    }else if(grepl("sens1",tag)){
      # sensitivity model 1: time = scan, quadratic growth, gamlssIA regression
      myTab_SNPAssocs_A = GetAssociation_v2(data = myTabX_long, 
                                         method = AssocModelX, 
                                         genotypes = G, 
                                         dep_var_name = "A",
                                         time_var_name = "scan",
                                         growth_var = set_growth, getIA = T)
    }else if(grepl("sens2",tag)){
      # sensitivity model 2: time = age, linear growth, gamlssIA regression
      myTab_SNPAssocs_A = GetAssociation_v2(data = myTabX_long, 
                                         method = AssocModelX, 
                                         genotypes = G, 
                                         dep_var_name = "A",
                                         time_var_name = "age",
                                         growth_var = "linear", getIA = T)
    }else if(grepl("sens3",tag)){
      # sensitivity model 3: time = age, quadratic growth, gamlss without SNP x time interaction 
      myTab_SNPAssocs_A = GetAssociation_v2(data = myTabX_long, 
                                         method = AssocModelX, 
                                         genotypes = G, 
                                         dep_var_name = "A",
                                         time_var_name = "age",
                                         growth_var = set_growth, getIA = F)
    }else if(grepl("sens4",tag)){
      # sensitivity model 4: time = age, quadratic growth, gamlss with time-independent sigma function
      myTab_SNPAssocs_A = GetAssociation_v2(data = myTabX_long, 
                                         method = AssocModelX, 
                                         genotypes = G, 
                                         dep_var_name = "A",
                                         time_var_name = "age",
                                         growth_var = set_growth, getIA = T,sigmaCovars = "SNP")
    }else if(grepl("sens5",tag)){
      # sensitivity model 5: time = age, quadratic growth, gamlss with SNP-independent sigma function
      myTab_SNPAssocs_A = GetAssociation_v2(data = myTabX_long, 
                                            method = AssocModelX, 
                                            genotypes = G, 
                                            dep_var_name = "A",
                                            time_var_name = "age",
                                            growth_var = set_growth, getIA = T,sigmaCovars = "time")
    }else if(grepl("sens6",tag)){
      # sensitivity model 6: time = age, quadratic growth, linear mixed model
      myTab_SNPAssocs_A = GetAssociation_v2(data = myTabX_long, 
                                         method = "linMixed", 
                                         genotypes = G, 
                                         dep_var_name = "A",
                                         time_var_name = "age",
                                         growth_var = set_growth, getIA = T)
    }

    #' Check 1: are the SNPs distinct or not? 
    # myTab_SNPAssocs_A[,table(exposure, pval<5e-8)]
    # myTab_SNPAssocs_A[SNP<11, table(exposure, pval<5e-8)]
    # myTab_SNPAssocs_A[SNP>=11,table(exposure, pval<5e-8)]
    # myTab_SNPAssocs_A[SNP>=11,table(exposure, pval<1e-6)]
    
    #' Check 2: are the t-values/betas correlated or not?
    # cor.test(myTab_SNPAssocs_A[exposure=="mean",beta],
    #          myTab_SNPAssocs_A[exposure=="slope",beta])
    # cor.test(myTab_SNPAssocs_A[exposure=="mean",beta],
    #          myTab_SNPAssocs_A[exposure=="var",beta])

    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(myTab_SNPAssocs_A,
           #myTab_SNPAssocs_Z,myTab_SNPAssocs_C,
           file = paste0(outdir_sim, "04_AssociationData.RData"))
    }
  }
  
  #' ## Step 5: Simulate outcome Y ####
  #' ***
  {
    myTabY = data.table(ID = unique(myTabX_long$ID))
    if(set_time == "scan"){
      myTabY[,age := 6]
    }else{
      myTabY[,age := rnorm(n=n_samples,mean=40,sd=1)]
    }
    
    myTabY[,Y1 := rnorm(n=n_samples,mean = Y_mean_random,sd = Y_var_random)]
    myTabY[,Y2 := Y_alpha[1] * AS1 + rnorm(n=n_samples,mean = Y_mean_random,sd = sd(AS1))]
    # myTabY[,Y3 := Y_alpha[2] * AS2 + 
    #          rnorm(n=n_samples,mean = Y_mean_random,sd = sd(AS2))]
    # myTabY[,Y4 := Y_alpha[1] * AS1 + Y_alpha[2] * AS2 + 
    #          rnorm(n=n_samples,mean = Y_mean_random,sd = sd(AS2)*sd(AS1))]
    myTabY[,Y3 := Y_alpha[2] * AS2 * age + 
             rnorm(n=n_samples,mean = Y_mean_random,sd = sd(AS2))]
    myTabY[,Y4 := Y_alpha[1] * AS1 + Y_alpha[2] * AS2 * age + 
             rnorm(n=n_samples,mean = Y_mean_random,sd = sd(AS2)*sd(AS1))]
    
    y_dummy = rbinom(n = n_samples, size = 1, prob = plogis(scale(myTabY$Y1)))
    myTabY[,Y1_bin := y_dummy]
    y_dummy = rbinom(n = n_samples, size = 1, prob = plogis(scale(myTabY$Y2)))
    myTabY[,Y2_bin := y_dummy]
    y_dummy = rbinom(n = n_samples, size = 1, prob = plogis(scale(myTabY$Y3)))
    myTabY[,Y3_bin := y_dummy]
    y_dummy = rbinom(n = n_samples, size = 1, prob = plogis(scale(myTabY$Y4)))
    myTabY[,Y4_bin := y_dummy]
    
    #' Check AS association
    # summary(lm(myTabY$Y1 ~ AS1 + AS2))
    # summary(lm(myTabY$Y2 ~ AS1 + AS2))
    # summary(lm(myTabY$Y3 ~ AS1 + AS2))
    # summary(lm(myTabY$Y4 ~ AS1 + AS2))
    # summary(lm(myTabY$Y3 ~ AS1 + AS2:myTabY$age))
    # summary(lm(myTabY$Y4 ~ AS1 + AS2:myTabY$age))
    # 
    # summary(glm(myTabY$Y1_bin ~ AS1 + AS2),family="binomial")
    # summary(glm(myTabY$Y2_bin ~ AS1 + AS2),family="binomial")
    # summary(glm(myTabY$Y3_bin ~ AS1 + AS2),family="binomial")
    # summary(glm(myTabY$Y4_bin ~ AS1 + AS2),family="binomial")
    # summary(glm(myTabY$Y3_bin ~ AS1 + AS2:myTabY$age),family="binomial")
    # summary(glm(myTabY$Y4_bin ~ AS1 + AS2:myTabY$age),family="binomial")
    
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
    myVars = names(myTabY)[-c(1,2)]
    
    myTab_SNPAssocs = foreach(i = 1:length(myVars))%do%{
      #i=1
      myMethod = "linReg" 
      if(grepl("_bin",myVars[i])==T) myMethod = "glm"
      
      myTab_GY1 = GetAssociation_v2(data = myTabY, 
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
    #' Run analysis
    MR_analysis_A = MVMR_jp_sim(data_GX_long = myTab_SNPAssocs_A,
                                data_GY_long = myTab_SNPAssocs_Y, 
                                filterBadSNPs = MR_filterBadSNPs,
                                filterBadSNPs_threshold = MR_filterBadSNPs_treshold)
    MR_analysis_A[,type:=tag]
    MRTab = copy(MR_analysis_A)
      
    #' Do some correction
    MRTab[,beta_IVW2 := beta_IVW]
    MRTab[,SE_IVW2 := SE_IVW]
    
    mean_age_outcome = myTabY[,mean(age)]
    MRTab[exposure=="slope",beta_IVW2 := beta_IVW2/mean_age_outcome]
    MRTab[exposure=="slope",SE_IVW2 := SE_IVW2/mean_age_outcome]
    
    #' Get some plots
    if((do_plotting == T & save_data_perSim==T) | s %in% counter){
      data_hlines = data.frame(exposure = unique(MRTab$exposure)[1:2],
                               mylines = c(Y_alpha))
      ggp3 = ggplot(MRTab[!grepl("bin",outcome)], aes(x=outcome, y=beta_IVW2,shape=type )) + 
        facet_wrap(~exposure,scales = "free")+
        geom_hline(data = data_hlines, aes(yintercept = mylines, col="red"), show.legend = FALSE) +
        geom_hline(yintercept = 0,col="black")+
        geom_point()+
        geom_errorbar(aes(ymin=beta_IVW2 -1.96*SE_IVW2, ymax=beta_IVW2 +1.96*SE_IVW2), width=.2)
      
      png(filename = paste0(outdir_sim, "07_MREstimates.png"), 
          width = 1200, height = 600, res=125)
      print(ggp3)
      dev.off()
    }
    
    #' Save data
    if(save_data_perSim == T | s %in% counter){
      save(MRTab,
           file = paste0(outdir_sim, "07_MVMRResults.RData"))
    }
    
  }
  
  MRTab[,n_sim := s]
  MRTab
}
SimTab = rbindlist(SimTab,fill=T)
save(SimTab, file = paste0("../",outfiles_dir,"/",outfiles_prefix, "_SimAll.RData"))

#' # Session Info ####
#' ***

sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
