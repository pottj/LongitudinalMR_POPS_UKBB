#' ---
#' title: "Check outcome data"
#' subtitle: "Longitudinal MR"
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
#' I want to test CAD (sex-combined, data taken from Aragam et al. (2022)) and coronary atherosclerosis (sex-stratified, Neale lab (2021?)) as outcomes. 
#' 
#' Here, I load the data and reduce them to the IVs for TC. I check the alleles and allele frequency. Every base position should be in hg19 / GRCh37 (our BSU UKB data is in hg19, Aragam data was downloaded using the GRCh37 build, Neale lab is using UKB, again hg19).  
#' 
#' - Neale lab: Variant identifier in the form "chr:pos:ref:alt", where "ref" is aligned to the forward strand of GRCh37 and "alt" is the effect allele (use this to join with variant annotation file).
#' - Plink: ALT is effect allele
#' 
#' **ALT == EA, REF == OA** in both Neale and my data. Aragam uses explicitly *effect_allele* as column name. 
#' 
#' # Init ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile_HPC.R")
.libPaths()

#' # Load HR SNP data ####
#' ***
#' I need the 151 SNPs, and the four data sets.
#' 
load(paste0(UKB_phenotypes_filtered,"/01_Prep_06_SNPData_filtered_sens.RData"))
pvar2[,chrPosEAOA := paste(CHR,POS,ALT,REF,sep=":")]
pvar2[,chrPosOAEA := paste(CHR,POS,REF,ALT,sep=":")]

#' # Load UKB CAD sex-combined ####
#' ***
pathData = gsub("~","../../../",pathData)
UKB_CAD_comb = fread(paste0(pathData,"I9_CORATHER.gwas.imputed_v3.both_sexes.tsv.bgz"))
UKB_CAD_comb
table(is.element(pvar2$chrPosEAOA,UKB_CAD_comb$variant))
table(is.element(pvar2$chrPosOAEA,UKB_CAD_comb$variant))

UKB_CAD_comb = UKB_CAD_comb[variant %in% pvar2$chrPosEAOA]
pvar3 = copy(pvar2)
pvar3 = pvar3[chrPosEAOA %in% UKB_CAD_comb$variant]

UKB_CAD_comb[,EA := gsub(".*:","",variant)]
UKB_CAD_comb[,OA := gsub("[0123456789]*:","",variant)]
UKB_CAD_comb[,OA := substr(OA,1,1)]
UKB_CAD_comb[,table(OA == EA)]
UKB_CAD_comb[EA == minor_allele,EAF := minor_AF]
UKB_CAD_comb[EA != minor_allele,EAF := 1-minor_AF]

table(pvar3$chrPosEAOA == UKB_CAD_comb$variant)
plot(pvar3$EAF,UKB_CAD_comb$EAF)
plot(pvar3$MAF,UKB_CAD_comb$minor_AF)

UKB_CAD_comb[, EAF2 :=1-EAF]
UKB_CAD_comb[, beta2 := beta*(-1)]
UKB_CAD_comb[, tstat2 := tstat*(-1)]
plot(pvar3$EAF,UKB_CAD_comb$EAF2)

UKB_CAD_comb[,rsID := pvar3$ID]
UKB_CAD_comb[,phenotype := "UKB_CAD_combined"]

#' which columns do I want? 
names(UKB_CAD_comb)
myAssocs_Y1 = copy(UKB_CAD_comb)
myAssocs_Y1 = myAssocs_Y1[,c(1,19,13,14,16,20,6,17,10,18,12)]
names(myAssocs_Y1) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' # Load UKB CAD females ####
#' ***
UKB_CAD_females = fread(paste0(pathData,"I9_CORATHER.gwas.imputed_v3.female.tsv.bgz"))
UKB_CAD_females
table(is.element(pvar2$chrPosEAOA,UKB_CAD_females$variant))
table(is.element(pvar2$chrPosOAEA,UKB_CAD_females$variant))

UKB_CAD_females = UKB_CAD_females[variant %in% pvar2$chrPosEAOA]
pvar3 = copy(pvar2)
pvar3 = pvar3[chrPosEAOA %in% UKB_CAD_females$variant]

UKB_CAD_females[,EA := gsub(".*:","",variant)]
UKB_CAD_females[,OA := gsub("[0123456789]*:","",variant)]
UKB_CAD_females[,OA := substr(OA,1,1)]
UKB_CAD_females[,table(OA == EA)]
UKB_CAD_females[EA == minor_allele,EAF := minor_AF]
UKB_CAD_females[EA != minor_allele,EAF := 1-minor_AF]

table(pvar3$chrPosEAOA == UKB_CAD_females$variant)
plot(pvar3$EAF,UKB_CAD_females$EAF)
plot(pvar3$MAF,UKB_CAD_females$minor_AF)

UKB_CAD_females[, EAF2 :=1-EAF]
UKB_CAD_females[, beta2 := beta*(-1)]
UKB_CAD_females[, tstat2 := tstat*(-1)]
plot(pvar3$EAF,UKB_CAD_females$EAF2)

UKB_CAD_females[,rsID := pvar3$ID]
UKB_CAD_females[,phenotype := "UKB_CAD_females"]

#' which columns do I want? 
names(UKB_CAD_females)
myAssocs_Y2 = copy(UKB_CAD_females)
myAssocs_Y2 = myAssocs_Y2[,c(1,19,13,14,16,20,6,17,10,18,12)]
names(myAssocs_Y2) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' # Load UKB CAD males ####
#' ***
UKB_CAD_males = fread(paste0(pathData,"I9_CORATHER.gwas.imputed_v3.male.tsv.bgz"))
UKB_CAD_males
table(is.element(pvar2$chrPosEAOA,UKB_CAD_males$variant))
table(is.element(pvar2$chrPosOAEA,UKB_CAD_males$variant))

UKB_CAD_males = UKB_CAD_males[variant %in% pvar2$chrPosEAOA]
pvar3 = copy(pvar2)
pvar3 = pvar3[chrPosEAOA %in% UKB_CAD_males$variant]

UKB_CAD_males[,EA := gsub(".*:","",variant)]
UKB_CAD_males[,OA := gsub("[0123456789]*:","",variant)]
UKB_CAD_males[,OA := substr(OA,1,1)]
UKB_CAD_males[,table(OA == EA)]
UKB_CAD_males[EA == minor_allele,EAF := minor_AF]
UKB_CAD_males[EA != minor_allele,EAF := 1-minor_AF]

table(pvar3$chrPosEAOA == UKB_CAD_males$variant)
plot(pvar3$EAF,UKB_CAD_males$EAF)
plot(pvar3$MAF,UKB_CAD_males$minor_AF)

UKB_CAD_males[, EAF2 :=1-EAF]
UKB_CAD_males[, beta2 := beta*(-1)]
UKB_CAD_males[, tstat2 := tstat*(-1)]
plot(pvar3$EAF,UKB_CAD_males$EAF2)

UKB_CAD_males[,rsID := pvar3$ID]
UKB_CAD_males[,phenotype := "UKB_CAD_males"]

#' which columns do I want? 
names(UKB_CAD_males)
myAssocs_Y3 = copy(UKB_CAD_males)
myAssocs_Y3 = myAssocs_Y3[,c(1,19,13,14,16,20,6,17,10,18,12)]
names(myAssocs_Y3) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' # Load Aragam CAD sex-combined ####
#' ***
CADdata = fread(paste0(pathData,"GCST90132314_buildGRCh37.tsv"))
CADdata
pvar2[,chrPosEAOA2 := paste0(CHR,":",POS,"_",ALT,"_",REF)]
pvar2[,chrPosOAEA2 := paste0(CHR,":",POS,"_",REF,"_",ALT)]
table(is.element(pvar2$chrPosEAOA2,CADdata$markername),
      is.element(pvar2$chrPosOAEA2,CADdata$markername))

CADdata = CADdata[markername %in% pvar2$chrPosEAOA2 | markername %in% pvar2$chrPosOAEA2]
pvar3 = copy(pvar2)
pvar3 = pvar3[chrPosEAOA2 %in% CADdata$markername | chrPosOAEA2 %in% CADdata$markername]

CADdata[,EA := toupper(effect_allele)]
CADdata[,OA := toupper(other_allele)]
CADdata[,table(OA == EA)]
CADdata[,EAF := effect_allele_frequency]

table(pvar3$chrPosEAOA2 == CADdata$markername)
table(pvar3$chrPosOAEA2 == CADdata$markername)
filt = pvar3$ALT == CADdata$EA
plot(pvar3$EAF[filt],CADdata$EAF[filt])
plot(pvar3$EAF[!filt],CADdata$EAF[!filt])

CADdata[!filt, EAF :=1-EAF]
CADdata[!filt, beta := beta*(-1)]
plot(pvar3$EAF,CADdata$EAF)

CADdata[,rsID := pvar3$ID]
CADdata[,phenotype := "Agaram_CAD_combined"]
CADdata[,tval := beta/standard_error]

#' which columns do I want? 
names(CADdata)
myAssocs_Y0 = copy(CADdata)
myAssocs_Y0 = myAssocs_Y0[,c(10,26,23,24,25,27,21,8,9,28,1)]
names(myAssocs_Y0) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 
myAssocs_Y0[,SNP := pvar3$chrPosEAOA]

#' # Load FinnGen + UKB CAD sex-combined ####
#' ***
pathData = gsub("~","../../../",pathData)
FinnGen_UKB = fread(paste0(pathData,"ukbb_summary_stats_I9_CORATHER_meta_out.tsv.gz"))
FinnGen_UKB
table(is.element(pvar2$chrPosEAOA,FinnGen_UKB$SNP))
table(is.element(pvar2$chrPosOAEA,FinnGen_UKB$SNP))
table(is.element(pvar2$ID,FinnGen_UKB$rsid))

FinnGen_UKB = FinnGen_UKB[rsid %in% pvar2$ID]
pvar3 = copy(pvar2)
pvar3 = pvar3[ID %in% FinnGen_UKB$rsid]
stopifnot(pvar3$ID == FinnGen_UKB$rsid)

table(pvar3$REF == FinnGen_UKB$ALT)
table(pvar3$ALT == FinnGen_UKB$ALT)
plot(pvar3$EAF,FinnGen_UKB$UKBB_af_alt)
plot(pvar3$EAF,FinnGen_UKB$FINNGEN_af_alt)

filt = pvar3$ALT != FinnGen_UKB$ALT
FinnGen_UKB[, EAF2 :=UKBB_af_alt]
FinnGen_UKB[filt, EAF2 :=1-UKBB_af_alt]
FinnGen_UKB[, beta2 := UKBB_beta]
FinnGen_UKB[filt, beta2 := UKBB_beta*(-1)]
FinnGen_UKB[,tstat := beta2/UKBB_sebeta]
plot(pvar3$EAF,FinnGen_UKB$EAF2)
FinnGen_UKB[,phenotype := "FinnGen_UKB_CAD_combined"]

#' which columns do I want? 
names(FinnGen_UKB)
myAssocs_Y4 = copy(FinnGen_UKB)
myAssocs_Y4 = myAssocs_Y4[,c(5,22,4,3,23,26,16,24,13,25,14)]
names(myAssocs_Y4) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' ### Combine data
myAssocs_Y = rbind(myAssocs_Y0,myAssocs_Y1,myAssocs_Y2,myAssocs_Y3,myAssocs_Y4)

#' ## save as temporary files
save(myAssocs_Y, file = paste0("../results/01_Prep_07_CADsummaryStats_sens.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
