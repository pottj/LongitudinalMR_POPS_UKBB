MVMR_jp_POPS = function(data_exposure,exposure_name, data_outcome, outcome_name,getPlot,flag,GX_pval_treshold,corTab,corTab_threshold,sampleSize_GX=15000,random=F,getCondF=F,getUni=T){
  #debug
  # data_exposure=copy(myAssocs_X_long2)
  # data_outcome=copy(myAssocs_Y2)
  # exposure_name=myExposure
  # outcome_name=myOutcome
  # flag="main"
  # getPlot=F
  # GX_pval_treshold = 1
  # corTab = copy(LDTab)
  # corTab_threshold = 0.1
  # sampleSize_GX = 3*2996
  # random=T
  # getCondF=T
  # getUni = T
  
  data_GX = copy(data_exposure)
  data_GX = data_GX[phenotype == exposure_name,]
  
  # do some priority pruning
  mySNPs = unique(data_GX[pval<GX_pval_treshold,SNP])
  corTab= corTab[SNP1 %in% mySNPs & SNP2 %in% mySNPs,]
  corTab= corTab[value>corTab_threshold]
  corTab[,SNP2 := as.character(SNP2)]
  data_GX2 = copy(data_GX)
  data_GX2 = data_GX2[SNP %in% mySNPs]
  setorder(data_GX2,pval)
  data_GX3 = data_GX2[!duplicated(SNP),]
  
  data_GX3[,indep := NA]
  while(sum(is.na(data_GX3$indep))!=0){
    mySNPs2 = data_GX3[is.na(indep),SNP]
    mySNP2 = mySNPs2[1]
    cor2 = corTab[SNP1 %in% mySNP2 | SNP2 %in% mySNP2]
    if(dim(cor2)[1]==0){
      data_GX3[SNP == mySNP2,indep := T]
    }else{
      mySNPs3 = unique(c(cor2$SNP1,cor2$SNP2))
      mySNPs3 = mySNPs3[!is.element(mySNPs3,mySNP2)]
      data_GX3[SNP %in% mySNPs3,indep := F]
      data_GX3[SNP == mySNP2,indep := T]
    }
  }
  table(data_GX3$indep, data_GX3$type)
  
  # update my associated SNP list
  mySNPs = unique(data_GX3[pval<GX_pval_treshold & indep==T,SNP])
  
  if(length(mySNPs)>=5){
    data_GX = data_GX[SNP %in% mySNPs,]
    
    data_GX_wide = dcast(data_GX, dumID + SNP ~ type, value.var=c("beta","SE","tval","pval"))
    
    types = unique(data_GX$type)
    SNPs_per_type = c()
    for(i in 1:length(types)){
      dum1 = dim(data_GX[pval<GX_pval_treshold & type == types[i]])[1]
      SNPs_per_type = c(SNPs_per_type,dum1)
    }
    
    data_GY = copy(data_outcome)
    data_GY = data_GY[phenotype == outcome_name,]
    data_GY = data_GY[SNP %in% data_GX_wide$SNP,]
    matched = match(data_GX_wide$SNP,data_GY$SNP)
    data_GY = data_GY[matched,]
    stopifnot(data_GX_wide$SNP == data_GY$SNP)
    
    filt1 = grepl("beta",names(data_GX_wide))
    data_beta = copy(data_GX_wide)
    data_beta = data_beta[,filt1,with=F]
    filt2 = grepl("SE",names(data_GX_wide))
    data_SE = copy(data_GX_wide)
    data_SE = data_SE[,filt2,with=F]
    types = gsub("beta_","",names(data_GX_wide)[filt1])
    
    if(random==T){
      filt= rbinom(n=dim(data_GX_wide)[1],size = 1,prob = 0.05)
      mvmr_obj = mr_mvinput(bx = as.matrix(data_beta)[filt==1,],
                            bxse = as.matrix(data_SE)[filt==1,],
                            by = data_GY$beta_mean[filt==1], 
                            byse = data_GY$SE_mean[filt==1],
                            exposure = types,
                            outcome = outcome_name)
    }else{
      mvmr_obj = mr_mvinput(bx = as.matrix(data_beta),
                            bxse = as.matrix(data_SE),
                            by = data_GY$beta_mean, 
                            byse = data_GY$SE_mean,
                            exposure = types,
                            outcome = outcome_name)
    }

    if(getCondF==F){
      res2 = mr_mvivw(mvmr_obj) 
      res3 = data.table(setting = rep("multivariate",length(types)),
                        exposure = rep(exposure_name,length(types)),
                        exposure_type = c(res2@Exposure),
                        outcome = rep(res2@Outcome,length(types)),
                        NR_SNPs_total = rep(res2@SNPs,length(types)),
                        NR_SNPs_type = SNPs_per_type,
                        beta_IVW = c(res2@Estimate),
                        SE_IVW = c(res2@StdError),
                        pval_IVW = c(res2@Pvalue),
                        HeteroStat = rep(res2@Heter.Stat[1],length(types)),
                        HeteroStat_pval = rep(res2@Heter.Stat[2],length(types)))
      
    }else{
      res2 = mr_mvivw(mvmr_obj,nx = sampleSize_GX) 
      res3 = data.table(setting = rep("multivariate",length(types)),
                        exposure = rep(exposure_name,length(types)),
                        exposure_type = c(res2@Exposure),
                        outcome = rep(res2@Outcome,length(types)),
                        NR_SNPs_total = rep(res2@SNPs,length(types)),
                        NR_SNPs_type = SNPs_per_type,
                        beta_IVW = c(res2@Estimate),
                        SE_IVW = c(res2@StdError),
                        pval_IVW = c(res2@Pvalue),
                        condFstat = c(res2@CondFstat),
                        HeteroStat = rep(res2@Heter.Stat[1],length(types)),
                        HeteroStat_pval = rep(res2@Heter.Stat[2],length(types)))
    }
    
    if(getUni==T){
      dumTab1 = foreach(i = 1:length(types))%do%{
        # i=1
        filt1 = grepl(types[i],data_GX$type)
        SNPs = unique(data_GX[filt1 & pval<GX_pval_treshold,SNP])
        filt2 = is.element(data_GY$SNP,SNPs)
        
        if(length(SNPs)>1){
          mr_obj = mr_input(bx = as.matrix(data_beta)[filt2,i],
                            bxse = as.matrix(data_SE)[filt2,i],
                            by = data_GY$beta_mean[filt2],
                            byse = data_GY$SE_mean[filt2],
                            exposure = types[i],
                            outcome = outcome_name,)
          res2 = mr_ivw(mr_obj)
          res4 = data.table(setting = "univariate",
                            exposure = exposure_name,
                            exposure_type = res2@Exposure,
                            outcome = res2@Outcome,
                            NR_SNPs_total = dim(data_GY[filt2,])[1],
                            NR_SNPs_type = dim(data_GY[filt2,])[1],
                            beta_IVW = res2@Estimate,
                            SE_IVW = res2@StdError,
                            pval_IVW = res2@Pvalue,
                            condFstat = c(res2@Fstat),
                            HeteroStat = res2@Heter.Stat[1],
                            HeteroStat_pval = res2@Heter.Stat[2])
          res4
          
        }else{
          res4 = data.table(setting = "univariate",
                            exposure = exposure_name,
                            exposure_type = types[i],
                            outcome = outcome_name,
                            NR_SNPs_total = dim(data_GY[filt2,])[1],
                            NR_SNPs_type = dim(data_GY[filt2,])[1])
        }
        res4
      }
      res4 = rbindlist(dumTab1,fill=T)
    }

    if(getPlot==T & getUni==T){
      # make a plot
      plotData = copy(data_GX)
      matched = match(plotData$SNP,data_GY$SNP)
      plotData[,myY := data_GY[matched,beta_mean]]
      
      matched = match(plotData$SNP,data_GX3$SNP)
      plotData[,myColor := data_GX3[matched,type]]
      
      data_hline_multi = data.table(type = res3$exposure_type,
                                    line = res3$beta_IVW)
      data_hline_uni = data.table(type = res4$exposure_type,
                                  line = res4$beta_IVW)
      
      myPlot = ggplot(data = plotData, aes(x = beta, y = myY, color = as.factor(myColor))) + 
        geom_point() +
        facet_wrap(~type,scales = "free")+        
        geom_abline(data = data_hline_multi, aes(slope = line,intercept=0,linetype="multivariate")) +
        geom_abline(data = data_hline_uni, aes(slope = line,intercept=0,linetype="univariate")) +
        scale_linetype_manual("IVW estimate",values=c("univariate"=2,"multivariate"=1)) +
        labs(x=paste0("SNP effect on ",exposure_name), 
             y=paste0("SNP effect on ",outcome_name),
             color="Lowest pval in ")
      
      #filename1 = paste0("../figures/08_7_MVMR_top20_",flag,"_",exposure_name,"_",outcome_name,"_",tag,".tiff")
      #tiff(filename = filename1,width = 2250, height = 1125, res=200, compression = 'lzw')
      print(myPlot)
      #dev.off()
    }
    
    res = res3
    if(getUni==T) res = rbind(res3,res4,fill=T)
    res[,ID := flag]
    
  }else{
    res = data.table(exposure = exposure_name,
                     outcome = outcome_name,
                     NR_SNPs_total = length(mySNPs))
  }
  return(res)
  
}

