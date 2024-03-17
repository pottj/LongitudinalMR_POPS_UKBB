#' ---
#' title: "Check candidate SNPs from PGS"
#' subtitle: "Longitudinal MVMR in POPS"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' In the real data application of my longitudinal MVMR analyses, I want to use SNPs affecting (estimated) birth weight.
#' 
#' Here, I load the list of variants contributing to the PGS for birth weight (downloaded from the PGS Catalog on 13/11/2023).
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","24-",tag)
tag = gsub("-","",tag)

#' # Load and prep PGS data ####
#' ***
PGS_data = fread(paste0(pathData,"/PGS001226_hmPOS_GRCh38.txt.gz"))
head(PGS_data)
PGS_data[,table(chr_name)]
PGS_data[,table(chr_name == hm_chr)]
PGS_data = PGS_data[chr_name %in% c(1:22),]
PGS_data[chr_name != hm_chr,]
PGS_data = PGS_data[chr_name == hm_chr,]

PGS_data[,dumID1 := paste0("chr",chr_name,":",hm_pos,":",effect_allele,":",other_allele)]
PGS_data[,dumID2 := paste0("chr",chr_name,":",hm_pos,":",other_allele,":",effect_allele)]


#' # Loop per chromosome ####
#' ***
#' 
dumTab = foreach(i=1:22)%do%{
  #i=22
  message("Working on chromosome ",i," ...")
  
  # Step 1: load SNP information 
  POPS_data = fread(paste0(POPS_imputed_data,"chr",i,".info.gz"))
  
  # Step 2: match SNPs to PGS data
  filt1 = is.element(POPS_data$SNP, PGS_data$dumID2)
  filt2 = is.element(POPS_data$SNP, PGS_data$dumID1)
  filt = filt1 | filt2
  POPS_data = POPS_data[filt,]
  PGS_data2 = copy(PGS_data)
  matched = match(POPS_data$SNP, PGS_data2$dumID2)
  matched2 = match(POPS_data$SNP, PGS_data2$dumID1)
  matched[is.na(matched)] = matched2[!is.na(matched2)]
  PGS_data2 = PGS_data2[matched,]
  
  # Step 3: save results
  data = cbind(POPS_data,PGS_data2)
  data
  
}
matchedData = rbindlist(dumTab)
head(matchedData)

outFile1 = paste0("../results/01_Prep_01_SNPList_",tag,".RData")
save(matchedData, file = outFile1)

outFile2 = paste0("../results/01_Prep_01_SNPList_IDonly_",tag,".txt")
write.table(matchedData$SNP, file = outFile2,col.names = F,row.names = F, quote = F)

outFile3 = paste0(POPS_SNP_data,"01_Prep_01_SNPList_PGS_IDonly_",tag,".txt")
write.table(matchedData$SNP, file = outFile3,col.names = F,row.names = F, quote = F)

outFile4 = paste0(POPS_SNP_data,"01_Prep_01_SNPList_PGS_EA_",tag,".txt")
dummy = matchedData[,c("SNP","effect_allele","effect_weight")]
write.table(matchedData[,c("SNP","effect_allele","effect_weight")], 
            file = outFile4, col.names = F,row.names = F, quote = F,  sep="\t")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

