#' @title MVMR function to get all causal estimates over several outcomes
#' @description In this function, I test for each outcome an MVMR model (classical = as defined in the MendelianRandomization package, correction = using mrest_me function to correct for possible measurement error)
#' @param data_GX_long data.table with SNP associations with any exposure levels 
#' @param data_GY_long data.table with SNP associations with any outcome types 
#' @param filterBadSNPs A Boolean, default is T. Should bad SNPs (weak instruments) be filtered for the analysis?  
#' @param filterBadSNPs_threshold threshold value to filter bad SNPs in case of filterBadSNPs == T, default: 1e-06
#' @return data.table with columns for exposure, outcome, beta_IVW, SE_IVW, pval_IVW (1 = classic, 2 = correction)
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname MVMR_jp
#' @export
MVMR_jp_sim_gmm = function(data_GX_long, data_GY_long, filterBadSNPs = T, filterBadSNPs_threshold = 1e-6,sampleSize_GX=15000,sampleSize_GY){
  #debug
  # data_GX_long = copy(myTab_SNPAssocs_A)
  # data_GY_long = copy(myTab_SNPAssocs_Y)
  # filterBadSNPs = MR_filterBadSNPs
  # filterBadSNPs_threshold = MR_filterBadSNPs_treshold
  # sampleSize_GX = dim(myTabX_long)[1]
  # sampleSize_GY = dim(myTabY)[1]

  data_GX_wide = dcast(data_GX_long, SNP ~ exposure, value.var = c("beta","SE","tval","pval"))
  outcomes = unique(data_GY_long$outcome)
  exposures = unique(data_GX_long$exposure)
  NR_exposures = length(exposures)
  
  if(filterBadSNPs == T){
    filt = data_GX_long[, pval<filterBadSNPs_threshold]
    goodSNPs = data_GX_long[filt,SNP]
    goodSNPs = unique(goodSNPs)
    
  }else{
    goodSNPs = data_GX_long[,SNP]
    goodSNPs = unique(goodSNPs)
  }
  comment = paste("p-value threshold =",filterBadSNPs_threshold)

  if(length(goodSNPs)<=NR_exposures){
    filterBadSNPs_threshold = 0.05
    if(filterBadSNPs == T){
      filt = data_GX_long[, pval<filterBadSNPs_threshold]
      goodSNPs = data_GX_long[filt,SNP]
      goodSNPs = unique(goodSNPs)

    }else{
      goodSNPs = data_GX_long[,SNP]
      goodSNPs = unique(goodSNPs)
    }
    comment = paste("p-value threshold =",filterBadSNPs_threshold)
    
  }
  
  if(length(goodSNPs)<=NR_exposures){
    myTab_MR = data.table(exposure = rep(exposures,length(outcomes)),
                     outcome = rep(outcomes,each=length(exposures)),
                     NR_SNPs_total = rep(length(goodSNPs),length(outcomes)* length(exposures)),
                     comment = "not enough SNPs to do an MVMR")
    
  }else{
    data_GX_long = data_GX_long[SNP %in% goodSNPs,]
    data_GX_wide = data_GX_wide[SNP %in% goodSNPs,]
    data_GY_long = data_GY_long[SNP %in% goodSNPs,]
    
    types = unique(data_GX_long$exposure)
    SNPs_per_type = c()
    for(i in 1:length(types)){
      dum1 = dim(data_GX_long[pval<filterBadSNPs_threshold & exposure == types[i]])[1]
      SNPs_per_type = c(SNPs_per_type,dum1)
    }
    
    dumTab_MR = foreach(i = 1:length(outcomes))%do%{
      #i=1
      myOutcome = outcomes[i]
      myTab_GY = copy(data_GY_long)
      myTab_GY = myTab_GY[outcome == myOutcome,]
      
      filt1 = grepl("beta",names(data_GX_wide))
      data_beta = copy(data_GX_wide)
      data_beta = data_beta[,filt1,with=F]
      corTab = cor(data_beta)
      filt2 = grepl("SE",names(data_GX_wide))
      data_SE = copy(data_GX_wide)
      data_SE = data_SE[,filt2,with=F]
      exposures2 = gsub("beta_","",names(data_beta))
      
      if(NR_exposures>1){
        mvmr_obj = mr_mvinput(bx = as.matrix(data_beta),
                              bxse = as.matrix(data_SE),
                              by = myTab_GY$beta, 
                              byse = myTab_GY$SE,
                              exposure = exposures2,
                              outcome = myOutcome)
        
        result = tryCatch({
          res2 = mr_mvgmm(mvmr_obj,nx = sampleSize_GX,ny = sampleSize_GY,cor.x = corTab)    

        }, error = function(e) {
          res2 = mr_mvivw(mvmr_obj,nx = sampleSize_GX)
          return(res2)
        })
        res = data.table(exposure = c(result@Exposure),
                         outcome = rep(result@Outcome,NR_exposures),
                         NR_SNPs_total = rep(length(goodSNPs),NR_exposures),
                         NR_SNPs_type = SNPs_per_type,
                         beta_IVW = c(result@Estimate),
                         SE_IVW = c(result@StdError),
                         pval_IVW = c(result@Pvalue),
                         condFStats = c(result@CondFstat),
                         HeteroStat = rep(result@Heter.Stat[1],NR_exposures),
                         HeteroStat_pval = rep(result@Heter.Stat[2],NR_exposures),
                         comment = class(result)[1])

      }else{
        mr_obj = mr_input(bx = data_GX_wide$beta_mean,
                          bxse = data_GX_wide$SE_mean,
                          by = myTab_GY$beta,
                          byse = myTab_GY$SE,
                          exposure = exposures2,
                          outcome = myOutcome)
        res2 = mr_ivw(mr_obj)
        res = data.table(exposure = exposures2,
                         outcome = res2@Outcome,
                         NR_SNPs_total = rep(length(goodSNPs),NR_exposures),
                         NR_SNPs_type = SNPs_per_type,
                         beta_IVW = res2@Estimate,
                         SE_IVW = res2@StdError,
                         pval_IVW = res2@Pvalue,
                         HeteroStat = res2@Heter.Stat[1],
                         HeteroStat_pval = res2@Heter.Stat[2],
                         comment = comment)
        
      }
      res
    }

    myTab_MR = rbindlist(dumTab_MR)
    
  }
  
  return(myTab_MR)
}
