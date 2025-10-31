#' ---
#' title: "Check POPS PCs"
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
#' I want to use Jasmines PCs. According to her mail, she used a pruned set of 80,396 directly genotyped SNPs. All SNPs have MAF >5%, have a call rate >95%, and are autosomal. LD was calculated via sliding window of 50,000 bp per step, with LD r2 <0.1. The software used was PCAIR, which was adjusted for genetic relationship matrix (GRM).
#' 
# 
# snpset_POP <- snpgdsLDpruning(gds_POP_f, 
#                               method="corr", 
#                               slide.max.bp=50e5,
#                               ld.threshold=sqrt(0.1), 
#                               verbose=FALSE, 
#                               maf=.05,
#                               missing.rate=.05)
#
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
POPS_PCS = fread(paste0(POPS_SNP_data,"PCs_Fetal_3484_GBRLabel.txt"))
POPS_PCS[,1:6]

myGD_files = list.files(path = POPS_phenotypes_filtered,pattern = "01_Prep_02")
myGD_file = myGD_files[grepl("filtered",myGD_files)]
loaded2 = load(paste0(POPS_phenotypes_filtered,myGD_file))
loaded2
psam

#' # Match data ####
#' ***
table(is.element(psam$IID,POPS_PCS$IID))
matched = match(psam$IID,POPS_PCS$IID)
POPS_PCS = POPS_PCS[matched,]
table(psam$IID == POPS_PCS$IID)

names(POPS_PCS) = c("IID","ancestry",paste0("PC",1:32))

POPS_PCS[,sex := psam[,SEX]]
POPS_PCS[sex==1,sex2 := "male"]
POPS_PCS[sex==2,sex2 := "female"]
POPS_PCS[,table(sex2,is.na(ancestry))]
POPS_PCS[,FID := psam$FID]

#' # Check ancestry ####
#' ***
#' ## PC1 vs PC2
i=3
j=4
max_all_PC1 = max(POPS_PCS[,c(i),with=F])
max_all_PC2 = max(POPS_PCS[,c(j),with=F])
min_all_PC1 = min(POPS_PCS[,c(i),with=F])
min_all_PC2 = min(POPS_PCS[,c(j),with=F])
png(file=paste0("../results/_figures/01_Prep_03_PCs/POPS_PCA_ancestry.png"),width=900,height=600)
par(mfrow=c(1,1))

plot(0,0,col="white",
     xlim=c(min_all_PC1,max_all_PC1),
     ylim=c(min_all_PC2,max_all_PC2),
     main=paste0("PC ",i-2," vs PC ",j-2),xlab=paste0("PC ",i-2),ylab=paste0("PC ",j-2))
lines(POPS_PCS[is.na(ancestry),c(i,j),with=F],col=alpha("purple",0.5),
      type="p",pch=19,cex=1)
lines(POPS_PCS[!is.na(ancestry),c(i,j),with=F],col=alpha("darkgreen",0.5),
      type="p",pch=19,cex=1)
legend(min_all_PC1, max_all_PC2, legend=c("not GBR", "GBR"),
       col=c("purple", "darkgreen"), pch=19, cex=1, title = "Ancestry")

dev.off()

#' OK, the GBR subgroup is rather homogeneous, while the non-GBR samples are heterogeneous. That is what I was expecting. 
#' 
#' ## First 5 PCs
#' 
png(file=paste0("../results/_figures/01_Prep_03_PCs/POPS_PCA_ancestry_panel.png"),width=1920,height=1200)
par(mfrow=c(5,5))
par(mar=c(2,2,2,1))

for(i in 1:5){
  #i=1
  i=i+2
  for(j in 1:5){
    #j=2
    j=j+2
    max_all_PC1 = max(POPS_PCS[,c(i),with=F])
    max_all_PC2 = max(POPS_PCS[,c(j),with=F])
    min_all_PC1 = min(POPS_PCS[,c(i),with=F])
    min_all_PC2 = min(POPS_PCS[,c(j),with=F])
    
    plot(0,0,col="white",
         xlim=c(min_all_PC1,max_all_PC1),
         ylim=c(min_all_PC2,max_all_PC2),
         main=paste0("PC ",i-2," vs PC ",j-2),xlab=paste0("PC ",i-2),ylab=paste0("PC ",j-2))
    lines(POPS_PCS[is.na(ancestry),c(i,j),with=F],col=alpha("purple",0.5),
          type="p",pch=19,cex=1)
    lines(POPS_PCS[!is.na(ancestry),c(i,j),with=F],col=alpha("darkgreen",0.5),
          type="p",pch=19,cex=1)
    # legend(min_all_PC1, max_all_PC2, legend=c("not GBR", "GBR"),
    #        col=c("purple", "darkgreen"), pch=19, cex=1, title = "Ancestry")
    
  }
  
}
dev.off()

#' # Check sex ####
#' ***
#' Now I exclude all non-GBR samples and check if the remaining samples are rather homogeneous.
#'  
#' ## PC1 vs PC2
POPS_PCS2 = POPS_PCS[!is.na(ancestry)]

i=3
j=4
max_all_PC1 = max(POPS_PCS2[,c(i),with=F])
max_all_PC2 = max(POPS_PCS2[,c(j),with=F])
min_all_PC1 = min(POPS_PCS2[,c(i),with=F])
min_all_PC2 = min(POPS_PCS2[,c(j),with=F])
png(file=paste0("../results/_figures/01_Prep_03_PCs/POPS_PCA_sex.png"),width=900,height=600)
par(mfrow=c(1,1))

plot(0,0,col="white",
     xlim=c(min_all_PC1,max_all_PC1),
     ylim=c(min_all_PC2,max_all_PC2),
     main=paste0("PC ",i-2," vs PC ",j-2),xlab=paste0("PC ",i-2),ylab=paste0("PC ",j-2))
lines(POPS_PCS2[!is.na(ancestry) & sex == 2,c(i,j),with=F],col=alpha("darkred",0.5),
      type="p",pch=19,cex=1)
lines(POPS_PCS2[!is.na(ancestry) & sex == 1,c(i,j),with=F],col=alpha("steelblue",0.5),
      type="p",pch=19,cex=1)
legend("bottomleft", inset = 0.05, legend=c("women", "men"),
       col=c("darkred", "steelblue"), pch=19, cex=1, title = "Sex")
dev.off()

#' OK, this looks good: men and women should be similar.
#' 
#' ## First 5 PCs
#' 
png(file=paste0("../results/_figures/01_Prep_03_PCs/POPS_PCA_sex_panel.png"),width=1920,height=1200)
par(mfrow=c(5,5))
par(mar=c(2,2,2,1))

for(i in 1:5){
  #i=1
  i=i+2
  for(j in 1:5){
    #j=2
    j=j+2
    max_all_PC1 = max(POPS_PCS2[,c(i),with=F])
    max_all_PC2 = max(POPS_PCS2[,c(j),with=F])
    min_all_PC1 = min(POPS_PCS2[,c(i),with=F])
    min_all_PC2 = min(POPS_PCS2[,c(j),with=F])
    
    plot(0,0,col="white",
         xlim=c(min_all_PC1,max_all_PC1),
         ylim=c(min_all_PC2,max_all_PC2),
         main=paste0("PC ",i-2," vs PC ",j-2),xlab=paste0("PC ",i-2),ylab=paste0("PC ",j-2))
    lines(POPS_PCS2[!is.na(ancestry) & sex == 2,c(i,j),with=F],col=alpha("darkred",0.5),
          type="p",pch=19,cex=1)
    lines(POPS_PCS2[!is.na(ancestry) & sex == 1,c(i,j),with=F],col=alpha("steelblue",0.5),
          type="p",pch=19,cex=1)
    
  }
  
}

dev.off()

#' # Save data ####
#' ***
save(POPS_PCS, file = paste0(POPS_phenotypes_filtered,"/01_Prep_03_PCs.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
