MR_jp_POPS = function(data_exposure,exposure_name, data_outcome, outcome_name,outcome_type,flag,GX_pval_treshold,corTab,corTab_threshold){
  #debug
  # data_exposure=copy(myAssocs_X_long2)
  # data_outcome=copy(myAssocs_Y2)
  # exposure_name=myExposure
  # outcome_name=myOutcome
  # flag="dumID"
  # outcome_type = "mean"
  # GX_pval_treshold = 0.05
  # corTab = copy(LDTab)
  # corTab_threshold = 0.1

  data_GX = copy(data_exposure)
  data_GX = data_GX[phenotype == exposure_name,]
  data_GX = data_GX[type == outcome_type,]
  
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
    
    data_GY = copy(data_outcome)
    data_GY = data_GY[phenotype == outcome_name,]
    data_GY = data_GY[SNP %in% data_GX$SNP,]
    matched = match(data_GX$SNP,data_GY$SNP)
    data_GY = data_GY[matched,]
    stopifnot(data_GX$SNP == data_GY$SNP)
    
    filt1 = grepl("beta",names(data_GX))
    data_betaX = copy(data_GX)
    data_betaX = data_betaX[,filt1,with=F]
    filt2 = grepl("SE",names(data_GX))
    data_SEX = copy(data_GX)
    data_SEX = data_SEX[,filt2,with=F]
    filt3 = grepl(paste0("beta_",outcome_type),names(data_GY))
    data_betaY = copy(data_GY)
    data_betaY = data_betaY[,filt3,with=F]
    filt4 = grepl(paste0("SE_",outcome_type),names(data_GY))
    data_SEY = copy(data_GY)
    data_SEY = data_SEY[,filt4,with=F]
    
    mr_obj = mr_input(bx = as.matrix(data_betaX)[,1],
                      bxse = as.matrix(data_SEX)[,1],
                      by = as.matrix(data_betaY)[,1],
                      byse = as.matrix(data_SEY)[,1],
                      exposure = exposure_name,
                      outcome = paste0(outcome_name," - ",outcome_type))
    res2 = mr_ivw(mr_obj)
    res4 = data.table(setting = "univariate",
                      exposure = exposure_name,
                      outcome = outcome_name,
                      type = outcome_type,
                      NR_SNPs_total = dim(data_GY)[1],
                      beta_IVW = res2@Estimate,
                      SE_IVW = res2@StdError,
                      pval_IVW = res2@Pvalue,
                      HeteroStat = res2@Heter.Stat[1],
                      HeteroStat_pval = res2@Heter.Stat[2])
    res4[,ID := flag]
    res4
    }else{
      res4 = data.table(setting = "univariate",
                     exposure = exposure_name,
                     outcome = outcome_name,
                     type = outcome_type,
                     NR_SNPs_total = length(mySNPs))
  }
  return(res4)
  
}

