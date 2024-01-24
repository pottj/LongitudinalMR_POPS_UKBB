#' ---
#' title: "Get POPS PCs"
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
#' I want to get the principal components for all my samples of the POP study. Here, I create all the plink calls, which I will then copy-paste into a slurm script to be run. 
#' I will keep only those samples which are in my phenotype data, and only those SNPs with rsID on the autosomes (not critical - they will be pruned anyway).
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../SourceFile_HPC.R")
.libPaths()

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("202.-","23-",tag)
tag = gsub("-","",tag)

#' # Load data ####
#' ***
#' ## Filter SNPs ####
#' I only want the autosomes, and only those with an rsID
#' 
rslist<-fread(paste0(POPS_genotyped_data,"POPs_Final-Child.bim"),sep="\t",stringsAsFactors=F)
dim(rslist)
head(rslist)
myTab = copy(rslist)
myTab = myTab[V1 %in% c(1:22),]
table(is.na(myTab$V6))
table(grepl("rs",myTab$V2))
myTab2 = myTab[!grepl("rs",V2),]
myTab = myTab[grepl("rs",V2),]
dummy<-as.character(myTab$V2)
write.table(dummy,paste0(file="../temp/POPS_genotypedSNPs_filtered_",tag,".txt"),
            quote=F,row.names=F,col.names=F)

#' ## Filter samples ####
#' I only want those with phenotypic data
#' 
fam.data = fread(paste0(POPS_genotyped_data,"POPs_Final-Child.fam"),sep="\t",stringsAsFactors=F)
load("../data/IndividualLevelData/02_Outcome.RData")
table(is.element(fam.data$V1,myTab_Y$POPSID))
table(is.element(myTab_Y$POPSID,fam.data$V1))

fam.data.restr = fam.data[V1 %in% myTab_Y$POPSID]

write.table(fam.data.restr,file=paste0("../temp/POPS_genotypedSamples_filtered_",tag,".txt"),
            quote=F,row.names=F,col.names=F)

#' # PLINK calls ####
#' ***
#' ## Step 1: LD pruning ####
call1<-paste0("plink2",
             " --bfile ",POPS_genotyped_data,"POPs_Final-Child",
             " --extract ../temp/POPS_genotypedSNPs_filtered_",tag,".txt",
             " --keep ../temp/POPS_genotypedSamples_filtered_",tag,".txt",
             " --maf 0.05",
             " --indep-pairwise 50 5 0.1",
             " --out ../temp/POPS_PCA_pruning_filter_",tag)
call1

#' ## Step 2: PCA ####
#' 
call2<-paste0("plink2",
             " --bfile ",POPS_genotyped_data,"POPs_Final-Child",
             " --extract ../temp/POPS_PCA_pruning_filter_",tag,".prune.in",
             " --keep ../temp/POPS_genotypedSamples_filtered_",tag,".txt",
             " --pca",
             " --out ../temp/POPS_PCA_",tag)
call2

#' ## Step 3: Evaluate PCA ####
#' 
pca2values<-fread(paste0("../temp/POPS_PCA_",tag,".eigenval"))$V1
pca2vector<-fread(paste0("../temp/POPS_PCA_",tag,".eigenvec"),
                       stringsAsFactors=F,sep="\t")


#' Plot 1: Eigenvalue per rank
#' 
plotData = data.table(rank = 1:10,
                      eigenvalue = pca2values,
                      explVar = pca2values/sum(pca2values))
plotData[,cumSum := cumsum(explVar)]

plotData[cumSum<0.85,]

png(file=paste0("../temp/POPS_PCA_Eigenvalues_",tag,".png"),width=600,height=600)
ggplot(plotData, aes(x=rank)) +
  
  geom_line( aes(y=eigenvalue), linewidth=1.5, color="darkred") + 
  geom_line( aes(y=cumSum*pca2values[1]), linewidth=1.5, color="steelblue") +
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Eigenvalue",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./pca2values[1],name="Cumulative Sums of Eigenvalues")
  ) + 
  
  theme(
    axis.title.y = element_text(color = "darkred", size=13),
    axis.title.y.right = element_text(color = "steelblue", size=13)
  ) +
  scale_x_continuous(breaks=seq(0,10,1)) +
  
  ggtitle("POPS eigenvalues")

dev.off()

#' Plot 2: Eigenvalue 1 vs 2
#' 
names(pca2vector) = c("ID","ID2",paste0("PC",1:10))
table(pca2vector$ID == fam.data.restr$V1)
pca2vector[,sex := fam.data.restr$V5]
pca2vector[sex==1,sex2 := "male"]
pca2vector[sex==2,sex2 := "female"]

max_all = max(pca2vector[,3:12])
min_all = min(pca2vector[,3:12])

png(file=paste0("../temp/POPS_PCA_panel_",tag,".png"),width=1920,height=1200)
par(mfrow=c(5,5))
par(mar=c(2,2,2,1))

for(i in 1:5){
  #i=1
  i=i+2
  for(j in 1:5){
    #j=2
    j=j+2
    # max_all = max(pca2vector[,c(i,j),with=F])
    # min_all = min(pca2vector[,c(i,j),with=F])
    plot(0,0,col="white",
         xlim=c(min_all,max_all),
         ylim=c(min_all,max_all),
         main=paste0("PC",i-2," vs PC",j-2))
    lines(pca2vector[sex2=="female",c(i,j),with=F],col=alpha("darkred",0.5),
          type="p",pch=19,cex=1)
    lines(pca2vector[sex2=="male",c(i,j),with=F],col=alpha("steelblue",0.5),
          type="p",pch=19,cex=1)
    
  }

}

dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
