Venn_jp_POPS_SNPs = function(data_exposure,exposure_name,top20,SNPSets="distinct",GX_pval_treshold,corTab,corTab_threshold){
  #debug
  # data_exposure=copy(myAssocs_X_gamlssIA)
  # exposure_name="efwcomb"
  # top20 = T
  # SNPSets = "distinct"
  # GX_pval_treshold = 0.05
  # corTab = copy(LDTab)
  # corTab_threshold = 0.1
  
  data_GX = copy(data_exposure)
  data_GX = data_GX[phenotype == exposure_name,]
  
  data_GX_long1 = melt(data_GX,
                       id.vars=c("SNP", "phenotype","sampleSize"),
                       measure.vars=c("beta_mean", "beta_slope", "beta_var" ),
                       variable.name="type",value.name="beta")
  data_GX_long1[,type := gsub("beta_","",type)]
  data_GX_long2 = melt(data_GX,
                       id.vars=c("SNP", "phenotype","sampleSize"),
                       measure.vars=c("SE_mean", "SE_slope", "SE_var" ),
                       variable.name="type",value.name="SE")
  data_GX_long2[,type := gsub("SE_","",type)]
  data_GX_long3 = melt(data_GX,
                       id.vars=c("SNP", "phenotype","sampleSize"),
                       measure.vars=c("tval_mean", "tval_slope", "tval_var" ),
                       variable.name="type",value.name="tval")
  data_GX_long3[,type := gsub("tval_","",type)]
  data_GX_long4 = melt(data_GX,
                       id.vars=c("SNP", "phenotype","sampleSize"),
                       measure.vars=c("pval_mean", "pval_slope", "pval_var" ),
                       variable.name="type",value.name="pval")
  data_GX_long4[,type := gsub("pval_","",type)]
  data_GX_long = cbind(data_GX_long1,data_GX_long2[,5],data_GX_long3[,5],data_GX_long4[,5])
  
  if(top20 == F){
    # do some priority pruning
    mySNPs = unique(data_GX_long[pval<GX_pval_treshold,SNP])
    corTab= corTab[SNP1 %in% mySNPs & SNP2 %in% mySNPs,]
    corTab= corTab[value>corTab_threshold]
    corTab[,SNP2 := as.character(SNP2)]
    data_GX2 = copy(data_GX_long)
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
    data_GX = data_GX_long[SNP %in% mySNPs,]
    
    SNPList1 = data_GX[pval<0.05 & type=="mean",SNP]
    SNPList2 = data_GX[pval<0.05 & type=="slope",SNP]
    SNPList3 = data_GX[pval<0.05 & type=="var",SNP]
    
    qq3 = venn3(x1=SNPList1,y1=SNPList2,z1=SNPList3,mylabels = c("mean","slope","var"))
    
  }else{
    # do some priority pruning
    mySNPs = unique(data_GX[,SNP])
    corTab= corTab[SNP1 %in% mySNPs & SNP2 %in% mySNPs,]
    corTab= corTab[value>corTab_threshold]
    corTab[,SNP2 := as.character(SNP2)]
    
    data_GX2 = copy(data_GX_long)
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
    mySNPs = unique(data_GX3[indep==T,SNP])
    data_GX = data_GX[SNP %in% mySNPs,]
    
    setorder(data_GX_long,pval)
    data_GX_long[,flag1 := FALSE]
    data_GX_long[,type2 := ""]
    types = unique(data_GX_long$type)
    
    if(SNPSets == "overlap"){
      # overlaps allowed
      for(i in 1:length(types)){
        # i=1
        mySNPs1 = data_GX_long[type == types[i],SNP]
        mySNPs1 = mySNPs1[1:20]
        data_GX_long[SNP %in% mySNPs1,flag1 := TRUE]
        data_GX_long[SNP %in% mySNPs1,type2 := paste0(type2,types[i],sep=", ")]
      }
      data_GX_long = data_GX_long[flag1==T,]
    }else if(SNPSets == "distinct"){
      # no overlaps, selecting the first 20 SNPs with no overlap
      dummy2 = data_GX_long[,min(pval),by=SNP]
      dummy2[,dumID2 := paste(SNP,V1,sep="__")]
      data_GX_long[,dumID2 := paste(SNP,pval,sep="__")]
      dummy3 = data_GX_long[dumID2 %in% dummy2$dumID2,]
      setorder(dummy3,pval)
      
      for(i in 1:length(types)){
        # i=1
        mySNPs1 = dummy3[type == types[i],SNP]
        mySNPs1 = mySNPs1[1:20]
        data_GX_long[SNP %in% mySNPs1,flag1 := TRUE]
        data_GX_long[SNP %in% mySNPs1,type2 := paste0(type2,types[i],sep=", ")]
      }
      data_GX_long = data_GX_long[flag1==T,]
    }
    
    # update my associated SNP list
    SNPList1 = data_GX_long[pval<0.05 & type=="mean",SNP]
    SNPList2 = data_GX_long[pval<0.05 & type=="slope",SNP]
    SNPList3 = data_GX_long[pval<0.05 & type=="var",SNP]
    
    qq3 = venn3(x1=SNPList1,y1=SNPList2,z1=SNPList3,mylabels = c("mean","slope","var"))
  }
  return(qq3)
}

