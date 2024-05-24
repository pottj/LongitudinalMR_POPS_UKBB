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
MR_jp_sim = function(data_GX_long, data_GY_long, filterBadSNPs = T, filterBadSNPs_threshold = 1e-6,sampleSize_GX=15000){
  #debug
  # data_GX_long = copy(myTab_SNPAssocs_A)
  # data_GY_long = copy(myTab_SNPAssocs_Y)
  # filterBadSNPs = MR_filterBadSNPs
  # filterBadSNPs_threshold = 1e-6
  # sampleSize_GX = n_samples * n_times
  
  outcomes = unique(data_GY_long$outcome)
  exposures = c("mean","slope")
  
  dumTab = foreach(i = 1:length(exposures))%do%{
    #i=1
    myExposure = exposures[i]
    data_GX = copy(data_GX_long)
    data_GX = data_GX[exposure == myExposure,]

    if(filterBadSNPs == T){
      filt = data_GX[, pval<filterBadSNPs_threshold]
      goodSNPs = data_GX[filt,SNP]
      goodSNPs = unique(goodSNPs)
      
    }else{
      goodSNPs = data_GX[,SNP]
      goodSNPs = unique(goodSNPs)
    }
    comment = paste("p-value threshold =",filterBadSNPs_threshold)
    data_GX = data_GX[SNP %in% goodSNPs,]
    data_GY = copy(data_GY_long)
    data_GY = data_GY[SNP %in% goodSNPs,]
    
    dumTab_MR = foreach(j = 1:length(outcomes))%do%{
      #j=1
      myOutcome = outcomes[j]
      myTab_GY = copy(data_GY)
      myTab_GY = myTab_GY[outcome == myOutcome,]
      
      mr_obj = mr_input(bx = as.numeric(data_GX$beta),
                          bxse = as.numeric(data_GX$SE),
                          by = as.numeric(myTab_GY$beta),
                          byse = as.numeric(myTab_GY$SE),
                          exposure = myExposure,
                          outcome = myOutcome)
      res2 = mr_ivw(mr_obj)
      res = data.table(exposure = myExposure,
                       outcome = res2@Outcome,
                       NR_SNPs_total = res2@SNPs,
                         beta_IVW = res2@Estimate,
                         SE_IVW = res2@StdError,
                         pval_IVW = res2@Pvalue,
                         HeteroStat = res2@Heter.Stat[1],
                         HeteroStat_pval = res2@Heter.Stat[2],
                         comment = comment)
        
      res
    }
    
    myTab_MR = rbindlist(dumTab_MR)
    myTab_MR
  }
  dumTab = rbindlist(dumTab)
  return(dumTab)
}
  
  
    
  return(myTab_MR)
}
