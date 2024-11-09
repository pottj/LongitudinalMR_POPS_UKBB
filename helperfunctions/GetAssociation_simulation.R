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
GetAssociation = function(data,method,genotypes,dep_var_name,time_var_name,growth_var,getIA=T){
  # debug
  # genotypes = G
  # method = "gamlssIARS"
  # data = copy(myTabX_long)
  # dep_var_name = "A"
  # time_var_name = "age"
  # growth_var = set_growth
  # getIA=T
  
  stopifnot(method %in% c("linReg","glm","linMixed","gamlssIA","gamlssIARS","gamlssNoIA","gamlssNoVar"))
  
  # step 1: get number of SNPs to be tested
  SNPs_NR = dim(genotypes)[2]
  data[,time := get(time_var_name)]
  if(!is.element(method,c("glm","linReg"))) time_NR = unique(data$time)
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
      
      if(getIA==T){
        mod1 = lm(myVar ~ mySNP:time, data = data2)
      }else{
        mod1 = lm(myVar ~ mySNP, data = data2)
      }
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
    
  }else if(method == "glm"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = genotypes[,j]
      helper = data.table(ID = sample_NR,
                          SNP = Gj)
      data2 = copy(data)
      matched = match(data2$ID,helper$ID)
      data2[,mySNP := helper[matched,SNP]]
      data2[,myVar := get(dep_var_name)]
      
      if(getIA==T){
        mod1 = glm(myVar ~ mySNP:time, data = data2,family="binomial")
      }else{
        mod1 = glm(myVar ~ mySNP, data = data2,family="binomial")
      }
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
      
      if(growth_var == "linear" & getIA==T){
        modX_G = lmer(myVar ~ mySNP + time + (1|ID) + mySNP:time, data = data2)
      }else if(growth_var != "linear" & getIA==T){
        modX_G = lmer(myVar ~ mySNP + time + (1|ID) + mySNP:time + age2, data = data2)
      }else if(growth_var == "linear" & getIA==F){
        modX_G = lmer(myVar ~ mySNP + time + (1|ID), data = data2)
      }else if(growth_var != "linear" & getIA==F){
        modX_G = lmer(myVar ~ mySNP + time + (1|ID) + age2, data = data2)
      }
      dummy1 = summary(modX_G)$coef
      dummy1 = dummy1[grepl("mySNP",rownames(dummy1)),]
      
      if(getIA==T){
        res = data.table(SNP = c(j,j),
                         exposure = c("mean","slope"),
                         beta = dummy1[,1],
                         SE = dummy1[,2],
                         tval = dummy1[,3])
        res[,pval := pnorm(-abs(beta/SE))*2]
      }else if(getIA==F){
        res = data.table(SNP = c(j),
                         exposure = c("mean"),
                         beta = dummy1[1],
                         SE = dummy1[2],
                         tval = dummy1[3])
        res[,pval := pnorm(-abs(beta/SE))*2]
      }
        
      
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
        
        if(growth_var == "linear"){
          modX_G = gamlss(myVar ~ mySNP*time + random(x = as.factor(ID)), 
                          sigma.formula = ~mySNP + time, 
                          data = data2, family = "NO")
        }else if(growth_var != "linear"){
          modX_G = gamlss(myVar ~ mySNP*time + age2 + random(x = as.factor(ID)), 
                          sigma.formula = ~mySNP + time, 
                          data = data2, family = "NO")
        }
        dummy1 = summary(modX_G)
        dummy1 = dummy1[grepl("mySNP",rownames(dummy1)),]
        
        res = data.table(SNP = c(j,j,j),
                         exposure = c("mean","slope","var"),
                         beta = dummy1[,1],
                         SE = dummy1[,2],
                         tval = dummy1[,3],
                         pval = dummy1[,4])
        res
    }
    
  }else if(method == "gamlssIARS"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = genotypes[,j]
      helper = data.table(ID = sample_NR,
                          SNP = Gj)
      data2 <<- copy(data)
      matched = match(data2$ID,helper$ID)
      data2[,mySNP := helper[matched,SNP]]
      data2[,myVar := get(dep_var_name)]
      
      message("working on SNP ",j)
      result = tryCatch({
        modX_G = gamlss(formula = myVar ~ mySNP * time + age2 + re(random=~1|as.factor(ID)),
                        sigma.formula = ~ mySNP + time + re(random=~1|as.factor(ID)),
                        data=data2, family = "NO")    
        dummy1 = summary(modX_G)
        dummy1 = dummy1[grepl("mySNP",rownames(dummy1)),]
        
        res = data.table(SNP = c(j,j,j),
                         exposure = c("mean","slope","var"),
                         beta = dummy1[,1],
                         SE = dummy1[,2],
                         tval = dummy1[,3],
                         pval = dummy1[,4])
      }, error = function(e) {
        res = data.table(SNP = c(j,j,j),
                         exposure = c("mean","slope","var"),
                         beta = c(0,0,0),
                         SE = c(1,1,1),
                         tval = c(0,0,0),
                         pval = c(1,1,1))
        return(res)
      })
      message("finished SNP ",j)
      
      result
    }
    
  }else if(method == "gamlssNoIA"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = genotypes[,j]
      helper = data.table(ID = sample_NR,
                          SNP = Gj)
      data2 = copy(data)
      matched = match(data2$ID,helper$ID)
      data2[,mySNP := helper[matched,SNP]]
      data2[,myVar := get(dep_var_name)]
      
      if(growth_var == "linear"){
        modX_G = gamlss(myVar ~ mySNP + time + random(x = as.factor(ID)), 
                        sigma.formula = ~mySNP + time, 
                        data = data2, family = "NO")
      }else if(growth_var != "linear"){
        modX_G = gamlss(myVar ~ mySNP + time + age2 + random(x = as.factor(ID)), 
                        sigma.formula = ~mySNP + time, 
                        data = data2, family = "NO")
      }
      dummy1 = summary(modX_G)
      dummy1 = dummy1[grepl("mySNP",rownames(dummy1)),]
      
      res = data.table(SNP = c(j,j),
                       exposure = c("mean","var"),
                       beta = dummy1[,1],
                       SE = dummy1[,2],
                       tval = dummy1[,3],
                       pval = dummy1[,4])
      res
    }
    
  }else if(method == "gamlssNoVar"){
    modTab = foreach(j = 1:SNPs_NR)%do%{
      #j=1
      Gj = genotypes[,j]
      helper = data.table(ID = sample_NR,
                          SNP = Gj)
      data2 = copy(data)
      matched = match(data2$ID,helper$ID)
      data2[,mySNP := helper[matched,SNP]]
      data2[,myVar := get(dep_var_name)]
      
      if(growth_var == "linear"){
        modX_G = gamlss(myVar ~ mySNP*time + random(x = as.factor(ID)), 
                        sigma.formula = ~ time, 
                        data = data2, family = "NO")
      }else if(growth_var != "linear"){
        modX_G = gamlss(myVar ~ mySNP*time + age2 + random(x = as.factor(ID)), 
                        sigma.formula = ~ time, 
                        data = data2, family = "NO")
      }
      dummy1 = summary(modX_G)
      dummy1 = dummy1[grepl("mySNP",rownames(dummy1)),]
      
      res = data.table(SNP = c(j,j),
                       exposure = c("mean","slope"),
                       beta = dummy1[,1],
                       SE = dummy1[,2],
                       tval = dummy1[,3],
                       pval = dummy1[,4])
      res
    }
    
  }
  myAssocs = rbindlist(modTab)
  
  # step 3: return results
  return(myAssocs)
  
}
