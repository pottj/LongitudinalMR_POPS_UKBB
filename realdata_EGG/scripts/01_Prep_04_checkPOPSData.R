#' ---
#' title: "Check POPS data"
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
#' I receive the POPS data from Ulla (see email from 26/09/2023). This individual level data cannot be shared via github or any other repository. Hence, I stored it under /data/IndividualLevelData/, which is ignored by git.
#' 
#' Here, I want to load the data and check out the time-(in)dependent trajectories for the growth parameters: 
#' 
#' - mothers weight
#' - estimated fetal weight (EFW)
#' - head circumference (HC)
#' - femur length (FL)
#' - biparietal diameter (BPD)
#' - abdominal circumference (AC)
#' 
#' The data was already filtered by Ulla Sovio for. Exclusion criteria were: 
#' 
#' - No biometry available at any time point,
#' - No information on the mode of delivery,
#' - Preterm birth,
#' - Non-cephalic presentation at delivery,
#' - Prelabor CD,
#' - Antepartum stillbirth, or
#' - Preexisting diabetes.
#' 
#' Gestational diabetes might affect fetal growth and will be used as covariable only. 
#' 
#' 
#' # Initialize ####
#' ***
#' 
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

tag = "EGG"

#' # Load POS data ####
#' ***
#' ## Load genetic data ####
myPC_files = list.files(path = POPS_phenotypes,pattern = "01_Prep_03")
myPC_file = myPC_files[grepl(tag,myPC_files)]
loaded1 = load(paste0(POPS_phenotypes,myPC_file))
loaded1

myGD_files = list.files(path = POPS_phenotypes,pattern = "01_Prep_02")
myGD_file = myGD_files[grepl("filtered_EGG",myGD_files)]
loaded2 = load(paste0(POPS_phenotypes,myGD_file))
loaded2
pvar = pvar2
geno_mat = geno_mat2

myScores = list.files(path="../results/",pattern = "01_Prep_02_SNPList")
myScore = myScores[grepl("filtered_EGG",myScores)]
loaded3 = load(paste0("../results/",myScore))
SNPList = SNPList_filtered

stopifnot(psam$FID == POPS_PCS$FID)
stopifnot(pvar$ID == SNPList$SNP)
stopifnot(pvar$ID == colnames(geno_mat))

filt = grepl("OK",SNPList$comment_AF) & grepl("OK",SNPList$comment_R2) & grepl("OK",SNPList$comment_HWE)
table(filt)

PGS_BW = geno_mat[,filt] %*% SNPList[filt,beta2]
POPS_PCS[,PGS_BW := PGS_BW]

#' ## Load phenotype data ####
data = data.table(read_xlsx(paste0(POPS_phenotypes,"/POPS_data_extract_JP_26sep2023.xlsx")))
names(data)
names(data)[grepl("an_scan2_",names(data))]

table(is.element(data$POPSID,psam$FID))

#' The gestational date is not separated by "_", so I will change that column name
names(data)[c(15,40,65)]
names(data)[c(15,40,65)] = c("an_scan2_ga","an_scan3_ga","an_scan4_ga")

#' # Reorganize data ####
#' ***
#' Right now, I have all data in one table. I want to split the data a bit: 
#' 
#' - longitudinal exposure data: data per scan (26 parameters)
#' - outcome data: data at birth (sex, GA, birth weight, caesarean section)
#' 
#' ## Longitudinal exposure data ####
myVariables = names(data)[grepl("an_scan2",names(data))]

dumTab = foreach(i = 1:length(myVariables))%do%{
  #i=12
  myVar = myVariables[i]
  myVars = c(myVar,gsub("scan2","scan3",myVar),gsub("scan2","scan4",myVar))
  dummy = myVars
  dummy
}
names(dumTab) = gsub("an_scan2_","",myVariables)

myTab_X <- melt(data,
                id.vars=names(data)[c(1,87,90,96)],
                measure.vars=dumTab,
                variable.name="scan")
myTab_X = myTab_X[,c(1:5,15,6:14,16:31)]

#' ## Outcome data ####
#' ***
myTab_Y = copy(data)
myTab_Y = myTab_Y[,!grepl("an_scan",names(data)),with=F]
myTab_Y = myTab_Y[,c(1:3,20:22,4:19,23:24)]

#' Check if there are any pre-term CS
myTab_Y[,table(pn_elcs)]

#' Check if there are non-cephalic presentations at birth 
myTab_Y[,table(is.na(pt_presdel))]
myTab_Y[,table(pt_presdel)]

#' Okay, a lot of entries are NA. Those without NA look ok, I guess. But this I will check with Ulla and Gordon one more time. (non-cephalics should have been filtered by Ulla before she gave me the data) 
#' 
#' Check gestational diabetes 
myTab_Y[,table(pn_gdm_diet_medication,is.na(pn_gdm_ga))]
myTab_Y[,summary(pn_gdm_ga)]

#' Okay, there are 2 women with GDM with no GA of diagnosis. Maybe not that critical
#' 
#' # Filter data ####
#' ***
#' ## Hard filters ####
#' 
#' Filter samples if there is
#' 
#' - no genetic data available
#' - mismatch between genotyped sex and database sex
#' 
#' Filter time points if they are out of GA bounds:
#' 
#' - scan 2: GA must be $\in [18,23)$
#' - scan 3: GA must be $\in [26,31)$
#' - scan 4: GA must be $\in [34,39)$
#' 
myTab_Y = myTab_Y[POPSID %in% POPS_PCS$FID,]
myTab_X = myTab_X[POPSID %in% POPS_PCS$FID,]

matched = match(myTab_Y$POPSID, POPS_PCS$FID)
POPS_PCS = POPS_PCS[matched,]
table(myTab_Y$POPSID == POPS_PCS[,FID])
table(myTab_Y$pn_sex,POPS_PCS[,sex])
filt = myTab_Y$pn_sex == "MALE" & POPS_PCS$sex == 2
table(filt)
myTab_Y = myTab_Y[!filt,]
POPS_PCS = POPS_PCS[!filt,]
myTab_X = myTab_X[POPSID %in% POPS_PCS$FID,]

myTab_X[,ga_flag := F]
myTab_X[scan==1 & ga<18 | scan==1 & ga>=23,ga_flag := T]
myTab_X[scan==2 & ga<26 | scan==2 & ga>=31,ga_flag := T]
myTab_X[scan==3 & ga<34 | scan==3 & ga>=39,ga_flag := T]
myTab_X[ga_flag == T,]
myTab_X = myTab_X[ga_flag == F,]
myTab_X[,attend := NULL]
myTab_X[,ga_flag := NULL]

#' ## Soft filters ####
#' 
#' Later on, I want to be able to filter for 
#' 
#' - ancestry (GBR only vs all ancestry, PC adjusted)
#' - complete time points (3 time points available vs at least 2 time points available)
#' 
stopifnot(myTab_Y$POPSID == POPS_PCS$FID)
myTab_Y = cbind(myTab_Y,POPS_PCS[,c(2:7,38),with=F])
table(is.na(myTab_Y$ancestry))

matched = match(myTab_X$POPSID,POPS_PCS$FID)
stopifnot(myTab_X$POPSID == POPS_PCS[matched,FID])
myTab_X = cbind(myTab_X,POPS_PCS[matched,2:7,with=F])

# okay, I want the number of time points per exposure in myTab_Y - this will make sample filtering later easier

myVariables = names(myTab_X)[8:30]

for(i in 1:length(myVariables)){
  #i=1
  myVariable = myVariables[i]
  myTab_X[,myVar := get(myVariable)]
  dummy = myTab_X[!is.na(myVar),.N,by=POPSID]
  # print(myVariable)
  # print(table(dummy$N))
  matched = match(myTab_Y$POPSID,dummy$POPSID)
  myTab_Y[,myVar := dummy[matched,N]]
  myTab_Y[is.na(myVar),myVar:=0]
  setnames(myTab_Y,"myVar",myVariable)
  myTab_X[,myVar := NULL]
}
dummy = myTab_Y[,c(32:54),with=F]
test = rowSums(dummy)
table(test)

#' So there are 
#' 
#' - 23 samples, which only have one time point for all 23 exposures, 
#' - 2229 samples, which have no missing time point for all 23 exposures,
#' - 570 samples, which have missing time points for hc and bpd 
#' 
#' What happens when I only test the 11 primary and secondary exposures?
names(myTab_Y)[c(42,47,48)]
dummy2 = myTab_Y[,c(42,47,48),with=F]
test2 = rowSums(dummy2)
table(test2)
table(test2,is.na(myTab_Y$ancestry))

#' OK, when using the primary exposure *efwcomb* and its Z-score and centiles, I will have a sample size of 
#' 
#' - all samples with at least 2 time points: 2996 --> main analysis?!
#' - GBR samples with at least 2 time points: 2594
#' - all samples with all 3 time points: 2885 
#' - GBR samples with all 3 time points: 2509
#' 

names(myTab_Y)[c(36:39,49,50,53,54)]
dummy2 = myTab_Y[,c(36:39,49,50,53,54),with=F]
test2 = rowSums(dummy2)
table(test2)

#' OK, for the secondary exposures, the sample sizes will vary per variable. 
#'
#' # Some plots ####
#' ***
#' ## Plot exposure data 
myTab_X[pn_sex == "MALE" & pn_emcsall == 0 ,myCol := "boy - \n normal"]
myTab_X[pn_sex == "MALE" & pn_emcsall == 1 ,myCol := "boy - \n eCS"]
myTab_X[pn_sex == "FEMALE" & pn_emcsall == 0 ,myCol := "girl - \n normal"]
myTab_X[pn_sex == "FEMALE" & pn_emcsall == 1 ,myCol := "girl - \n eCS"]

myVariables = names(myTab_Y)[c(42,47,48,36:39,49,50,53,54)]

for(i in 1:(length(myVariables))){
  #i=1
  message("Working on variable: ",myVariables[i])
  myVar = myVariables[i]
  
  myTab_X2 = copy(myTab_X)
  myTab_X2[,value := get(myVar)]
  
  goodSamples = myTab_Y[get(myVar)>=2,POPSID]
  myTab_X2 = myTab_X2[POPSID %in% goodSamples,]
  myTab_X2 = myTab_X2[!is.na(ga),]
  myTab_X2 = myTab_X2[!is.na(value),]
  SampleSize = length(unique(myTab_X2$POPSID))
  
  ggp1  = ggplot(myTab_X2, aes(x=ga, y=value, col=myCol, group=POPSID)) +
    geom_point()+
    geom_line(show.legend = TRUE,aes(alpha=0.05)) + 
    labs(x="Gestational Week",y=myVar, color="Babys sex \n delivery",
         title = paste0(myVar,", sample size: ",SampleSize," (with at least 2 measurements)")) +
    theme(legend.position = "none") + theme_classic()+ 
    scale_colour_manual(values = c("darkblue","steelblue","darkred","firebrick1"))+
    guides(alpha="none")
    
  png(file=paste0("../results/_figures/01_Prep_04_POPS/Pheno_",myVar,"_filtered.png"),
      width=900,height=600)
  print(ggp1)
  dev.off()

}

#' ## Plot outcome data 
#' I want histograms for pn_bw, BW_SDS_Br1990 and BW_Centile_Br1990
#' 
plotData = copy(myTab_Y)
plotData = plotData[,c(1,10,9,18,12,15,16,31)]

plotData_long =  melt(plotData,
                      id.vars=c("POPSID","pn_ga_wk","pn_sex","pn_emcsall"),
                      measure.vars=c("pn_bw","BW_Centile_Br1990","BW_SDS_Br1990","PGS_BW"),
                      variable.name="birthweight",
                      value.name="measurement")

ggp2  = ggplot(plotData_long, aes(x=measurement, color=pn_sex)) +
  facet_wrap(~birthweight, ncol = 2,scales = "free") +
  geom_histogram(fill="white",alpha=0.5, position="identity")+
  # labs(x="Gestational Week",y=myVar, color="Babys sex \n delivery",
  #      title = paste0(myVar,", sample size: ",SampleSize," (with at least 2 measurements)")) +
  theme(legend.position = "none") + theme_classic()+ 
  # scale_colour_manual(values = c("darkblue","steelblue","darkred","firebrick1"))+
  guides(fill="none")

png(file=paste0("../results/_figures/01_Prep_04_POPS/Pheno_birthweight_filtered.png"),
    width=1500,height=700)
print(ggp2)
dev.off()


#' # Save data ####
#' ***
#' 
save(myTab_Y,file = paste0(POPS_phenotypes,"/01_Prep_04_Outcome_",tag,".RData"))
save(myTab_X,file = paste0(POPS_phenotypes,"/01_Prep_04_Exposure_",tag,".RData"))

psam = psam[FID %in% POPS_PCS$FID]
filt = is.element(rownames(geno_mat),psam$FID)
geno_mat = geno_mat[filt,]
dim(geno_mat)
save(psam,pvar,geno_mat, 
     file = paste0(POPS_phenotypes,"/01_Prep_04_SNPData_",tag,".RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

