#' @title Get SNP associations with exposure or outcome
#' @description This function estimates the beta coefficients and standard errors for all SNP-phenotype association necessary in my simulation setup
#' @param data data.table with phenotype column
#' @param method parameter indicating what type of phenotype is given in the phenotype column and which regression model should be used. Must be either "linReg", "linMixed", "gamlss".
#' @param genotypes Genotype matrix
#' @param dep_var_name name of phenotype column 
#' @return data.table containing one row per SNP from the genotype matrix (in case of "linReg") or two rows per SNP (in case of "linMixed": intercept and slope; in case of "gamlss": mu and sigma). Columns are SNP number, exposure (as given in dep_var_name), beta, SE, tval, and pval. 
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname GetAssociation
#' @export
GetAssociation = function(data,method,genotypes,dep_var_name){
  # debug
  # genotypes = G
  # method = AssocModel
  # data = copy(myTabX_long)
  # dep_var_name = X_type
  
  stopifnot(method %in% c("linReg","linMixed","gamlss","gamlssIA"))
  
  # step 1: get number of SNPs to be tested
  SNPs_NR = dim(genotypes)[2]
  time_NR = unique(data$time)
  sample_NR = unique(data$ID)
  
  # step 2: loop per SNP
  if(method == "linReg"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = genotypes[,j]
      helper = data.table(ID = sample_NR,
                          SNP = Gj)
      data2 = copy(data)
      matched = match(data2$ID,helper$ID)
      data2[,mySNP := helper[matched,SNP]]
      data2[,myVar := get(dep_var_name)]
      
      mod1 = lm(myVar ~ mySNP, data = data2)
      #summary(mod1)
      data2[, mySNP := NULL]
      data2[, myVar := NULL]
      
      res = data.table(SNP = j,
                       exposure = dep_var_name,
                       beta = summary(mod1)$coef[2,1],
                       SE = summary(mod1)$coef[2,2],
                       tval = summary(mod1)$coef[2,3],
                       pval = summary(mod1)$coef[2,4])
      res
    }
    
  }else if(method == "linMixed"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = genotypes[,j]
      helper = data.table(ID = sample_NR,
                          SNP = Gj)
      data2 = copy(data)
      matched = match(data2$ID,helper$ID)
      data2[,mySNP := helper[matched,SNP]]
      data2[,myVar := get(dep_var_name)]
      
      modX_G = lmer(myVar ~ mySNP*scale(age) + age2 + (1|ID), data = data2)
      #modX_G = lmer(myVar ~ mySNP*time + age2 + (1|ID), data = data2)
      dummy1 = summary(modX_G)$coef
      dummy1 = dummy1[grepl("mySNP",rownames(dummy1)),]
      
      res = data.table(SNP = c(j,j),
                       exposure = c("mean","slope"),
                       beta = dummy1[,1],
                       SE = dummy1[,2],
                       tval = dummy1[,3])
      res[,pval := pnorm(-abs(beta/SE))*2]
      res
    }
    
  }else if(method == "gamlss"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = genotypes[,j]
      helper = data.table(ID = sample_NR,
                          SNP = Gj)
      data2 = copy(data)
      matched = match(data2$ID,helper$ID)
      data2[,mySNP := helper[matched,SNP]]
      data2[,myVar := get(dep_var_name)]
      
      modX_G = gamlss(myVar ~ mySNP + scale(time) + random(x = as.factor(ID)), 
                      sigma.formula = ~mySNP, 
                      data = data2, family = "NO")
      dummy2 = summary(modX_G)
      dummy2 = dummy2[grepl("mySNP",rownames(dummy2)),]
      
      res = data.table(SNP = c(j,j),
                       exposure = c("mean","var"),
                       beta = dummy2[,1],
                       SE = dummy2[,2],
                       tval = dummy2[,3],
                       pval = dummy2[,4])
      res
    }
    
  }else if(method == "gamlssIA"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
        #j=1
        Gj = genotypes[,j]
        helper = data.table(ID = sample_NR,
                            SNP = Gj)
        data2 = copy(data)
        matched = match(data2$ID,helper$ID)
        data2[,mySNP := helper[matched,SNP]]
        data2[,myVar := get(dep_var_name)]
        
        modX_G = gamlss(myVar ~ mySNP*scale(time) + random(x = as.factor(ID)), 
                        sigma.formula = ~mySNP, 
                        data = data2, family = "NO")
        dummy3 = summary(modX_G)
        dummy3 = dummy3[grepl("mySNP",rownames(dummy3)),]
        
        res = data.table(SNP = c(j,j,j),
                         exposure = c("mean","slope","var"),
                         beta = dummy3[,1],
                         SE = dummy3[,2],
                         tval = dummy3[,3],
                         pval = dummy3[,4])
        res
    }
    
  }
  myAssocs = rbindlist(modTab)
  
  # step 3: return results
  return(myAssocs)
  
}
