#' ---
#' title: "Check outcome data (AF)"
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
#' I want to test AF (sex-combined, data taken from Miyazawa et al. (2023)) and AF (Diagnoses - main ICD10: I48 Atrial fibrillation and flutter) (sex-stratified, Neale lab (2021?)) as outcomes. 
#' 
#' Here, I load the data and reduce them to the IVs for HR. I check the alleles and allele frequency. Every base position should be in hg19 / GRCh37 (our BSU UKB data is in hg19, Miyazawa data was downloaded using the GRCh37 build, Neale lab is using UKB, again hg19).  
#' 
#' - Neale lab: Variant identifier in the form "chr:pos:ref:alt", where "ref" is aligned to the forward strand of GRCh37 and "alt" is the effect allele (use this to join with variant annotation file).
#' - Plink: ALT is effect allele
#' 
#' **ALT == EA, REF == OA** in both Neale and my data. Miyazawa uses explicitly *effect_allele* as column name. 
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
load(paste0(UKB_phenotypes_filtered,"/01_Prep_04_SNPData_filtered.RData"))
pvar2[,chrPosEAOA := paste(CHR,POS,ALT,REF,sep=":")]
pvar2[,chrPosOAEA := paste(CHR,POS,REF,ALT,sep=":")]

#' # Load UKB AF sex-combined ####
#' ***
UKB_AF_comb = fread(paste0(pathData,"I48.gwas.imputed_v3.both_sexes.tsv.bgz"))
UKB_AF_comb
table(is.element(pvar2$chrPosEAOA,UKB_AF_comb$variant))
table(is.element(pvar2$chrPosOAEA,UKB_AF_comb$variant))

UKB_AF_comb = UKB_AF_comb[variant %in% pvar2$chrPosEAOA]
pvar3 = copy(pvar2)
pvar3 = pvar3[chrPosEAOA %in% UKB_AF_comb$variant]

UKB_AF_comb[,EA := gsub(".*:","",variant)]
UKB_AF_comb[,OA := gsub("[0123456789]*:","",variant)]
UKB_AF_comb[,OA := substr(OA,1,1)]
UKB_AF_comb[,table(OA == EA)]
UKB_AF_comb[EA == minor_allele,EAF := minor_AF]
UKB_AF_comb[EA != minor_allele,EAF := 1-minor_AF]

table(pvar3$chrPosEAOA == UKB_AF_comb$variant)
plot(pvar3$EAF,UKB_AF_comb$EAF)
plot(pvar3$MAF,UKB_AF_comb$minor_AF)

UKB_AF_comb[, EAF2 :=1-EAF]
UKB_AF_comb[, beta2 := beta*(-1)]
UKB_AF_comb[, tstat2 := tstat*(-1)]
plot(pvar3$EAF,UKB_AF_comb$EAF2)

UKB_AF_comb[,rsID := pvar3$ID]
UKB_AF_comb[,phenotype := "UKB_AF_combined"]

#' which columns do I want? 
names(UKB_AF_comb)
myAssocs_Y1 = copy(UKB_AF_comb)
myAssocs_Y1 = myAssocs_Y1[,c(1,19,13,14,16,20,6,17,10,18,12)]
names(myAssocs_Y1) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' # Load UKB AF females ####
#' ***
UKB_AF_females = fread(paste0(pathData,"I48.gwas.imputed_v3.female.tsv.bgz"))
UKB_AF_females
table(is.element(pvar2$chrPosEAOA,UKB_AF_females$variant))
table(is.element(pvar2$chrPosOAEA,UKB_AF_females$variant))

UKB_AF_females = UKB_AF_females[variant %in% pvar2$chrPosEAOA]
pvar3 = copy(pvar2)
pvar3 = pvar3[chrPosEAOA %in% UKB_AF_females$variant]

UKB_AF_females[,EA := gsub(".*:","",variant)]
UKB_AF_females[,OA := gsub("[0123456789]*:","",variant)]
UKB_AF_females[,OA := substr(OA,1,1)]
UKB_AF_females[,table(OA == EA)]
UKB_AF_females[EA == minor_allele,EAF := minor_AF]
UKB_AF_females[EA != minor_allele,EAF := 1-minor_AF]

table(pvar3$chrPosEAOA == UKB_AF_females$variant)
plot(pvar3$EAF,UKB_AF_females$EAF)
plot(pvar3$MAF,UKB_AF_females$minor_AF)

UKB_AF_females[, EAF2 :=1-EAF]
UKB_AF_females[, beta2 := beta*(-1)]
UKB_AF_females[, tstat2 := tstat*(-1)]
plot(pvar3$EAF,UKB_AF_females$EAF2)

UKB_AF_females[,rsID := pvar3$ID]
UKB_AF_females[,phenotype := "UKB_AF_females"]

#' which columns do I want? 
names(UKB_AF_females)
myAssocs_Y2 = copy(UKB_AF_females)
myAssocs_Y2 = myAssocs_Y2[,c(1,19,13,14,16,20,6,17,10,18,12)]
names(myAssocs_Y2) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' # Load UKB AF males ####
#' ***
UKB_AF_males = fread(paste0(pathData,"I48.gwas.imputed_v3.male.tsv.bgz"))
UKB_AF_males
table(is.element(pvar2$chrPosEAOA,UKB_AF_males$variant))
table(is.element(pvar2$chrPosOAEA,UKB_AF_males$variant))

UKB_AF_males = UKB_AF_males[variant %in% pvar2$chrPosEAOA]
pvar3 = copy(pvar2)
pvar3 = pvar3[chrPosEAOA %in% UKB_AF_males$variant]

UKB_AF_males[,EA := gsub(".*:","",variant)]
UKB_AF_males[,OA := gsub("[0123456789]*:","",variant)]
UKB_AF_males[,OA := substr(OA,1,1)]
UKB_AF_males[,table(OA == EA)]
UKB_AF_males[EA == minor_allele,EAF := minor_AF]
UKB_AF_males[EA != minor_allele,EAF := 1-minor_AF]

table(pvar3$chrPosEAOA == UKB_AF_males$variant)
plot(pvar3$EAF,UKB_AF_males$EAF)
plot(pvar3$MAF,UKB_AF_males$minor_AF)

UKB_AF_males[, EAF2 :=1-EAF]
UKB_AF_males[, beta2 := beta*(-1)]
UKB_AF_males[, tstat2 := tstat*(-1)]
plot(pvar3$EAF,UKB_AF_males$EAF2)

UKB_AF_males[,rsID := pvar3$ID]
UKB_AF_males[,phenotype := "UKB_AF_males"]

#' which columns do I want? 
names(UKB_AF_males)
myAssocs_Y3 = copy(UKB_AF_males)
myAssocs_Y3 = myAssocs_Y3[,c(1,19,13,14,16,20,6,17,10,18,12)]
names(myAssocs_Y3) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 

#' # Load Miyazawa AF sex-combined ####
#' ***
AFdata = fread(paste0(pathData,"GCST90204201_buildGRCh37.tsv"))
AFdata
table(is.element(pvar2$ID,AFdata$variant_id))

AFdata = AFdata[variant_id %in% pvar2$ID]
pvar3 = copy(pvar2)
pvar3 = pvar3[ID %in% AFdata$variant_id]

AFdata[,EA := effect_allele]
AFdata[,OA := other_allele]
AFdata[,table(OA == EA)]
AFdata[,EAF := effect_allele_frequency]

table(pvar3$ID == AFdata$variant_id)
table(pvar3$ALT == AFdata$OA,
      pvar3$REF == AFdata$EA)
plot(pvar3$EAF,AFdata$EAF)

AFdata[, EAF :=1-EAF]
AFdata[, beta := beta*(-1)]
plot(pvar3$EAF,AFdata$EAF)
abline(0,1)

AFdata[,phenotype := "Miyazawa_AF_combined"]
AFdata[,tval := beta/standard_error]
AFdata[,sampleSize := 1244730]

#' which columns do I want? 
names(AFdata)
myAssocs_Y0 = copy(AFdata)
myAssocs_Y0 = myAssocs_Y0[,c(1,1,10:13,15,7,8,14,9)]
names(myAssocs_Y0) = c("SNP","rsID","EA","OA","EAF","phenotype","sampleSize","beta_mean","SE_mean","tval_mean","pval_mean") 
myAssocs_Y0[,SNP := pvar3$chrPosEAOA]
myAssocs_Y0[,pval_mean := as.numeric(pval_mean)]


#' ### Combine data
myAssocs_Y = rbind(myAssocs_Y0,myAssocs_Y1,myAssocs_Y2,myAssocs_Y3)

#' ## save as temporary files
save(myAssocs_Y, file = paste0("../results/01_Prep_07_AFsummaryStats.RData"))

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
