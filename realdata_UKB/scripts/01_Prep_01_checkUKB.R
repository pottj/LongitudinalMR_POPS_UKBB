#' ---
#' title: "Check LDLC variability"
#' subtitle: "Longitudinal MVMR - LDLC"
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
#' 
#' 1) Load UKB data from baseline and follow-up
#' 2) Load UKB GP data for TC and lipid lowering medication
#' 3) Merge data sets 
#' 4) Checks
#' 5) Save 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()

source("../../SourceFile.R")
.libPaths()

#' # UKB data from BL and FU ####
#' ***
#' 
#' ## Load data ####
#' 
#' - **exposure**: TC at baseline (BL) and follow-up (FU)
#' - **covariables**: ID, sex, date (BL + FU), menopausal status (BL + FU), medication (BL + FU), ancestry, age at recruitment (BL), genetic sex, genetic PCs, genetic kinship
#' - **outcome** (CAD): main diagnoses, secondary diagnoses, any diagnoses (BL only) 
#' 
{
  #' Load the first 20 individuals to get overview of all columns available
  myTab_20 = fread(UKB_phenotypes, 
                   header=TRUE, sep="\t",nrows = 20)
  myAnnot = data.table(colNm = names(myTab_20))
  myAnnot[,colNR := 1:18506]
  
  #' Get relevant parameters
  exposure = c("f.30690.0.0","f.30690.1.0")
  table(is.element(exposure,myAnnot$colNm))
  
  covars = c("f.eid", "f.31.0.0","f.53.0.0","f.53.1.0","f.2724.0.0","f.2724.1.0",
             paste0("f.20003.0.",c(0:47)),paste0("f.20003.1.",c(0:47)),
             "f.21000.0.0","f.21022.0.0","f.22001.0.0",paste0("f.22009.0.",1:10), "f.22021.0.0")
  table(is.element(covars,myAnnot$colNm))
  
  CAD_diagnoses = c(paste0("f.41202.0.",c(0:79)), 
                    paste0("f.41204.0.",c(0:209)),
                    paste0("f.41270.0.",c(0:258)))
  table(is.element(CAD_diagnoses,myAnnot$colNm))
  
  #' Reduce annotation to those columns
  myAnnot = myAnnot[colNm %in% c(exposure,covars,CAD_diagnoses),]
  x = myAnnot[,colNR]
  
  #' Load the UKB data - restricted to these columns
  myTab_cross <- fread(UKB_phenotypes, header=TRUE, sep="\t",select = x)
  
  #' Rename column names
  names(myTab_cross)
  names(myTab_cross) = c("ID","sex","date_BL","date_FU","menopause_BL","menopause_FU",
                         paste("medication_BL",1:48,sep="_"),
                         paste("medication_FU",1:48,sep="_"),
                         "ancestry","age","geneticSex",
                         paste("PC",1:10,sep="_"),
                         "kinship","TC_BL","TC_FU",
                         paste("diagnoses_main",1:80,sep="_"),
                         paste("diagnoses_sec",1:210,sep="_"),
                         paste("diagnoses_any",1:259,sep="_"))
  
  #' Save data
  save(myTab_cross, file = paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_unfiltered.RData"))
  
}

#' ## Filter data ####
#' 
#' - White British ancestry
#' - No pairwise kinship
#' - Consent still active
#' - Genetic sex = data base sex
#' 
{
  # load(paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_unfiltered.RData"))
  dim(myTab_cross)
  
  #' Filter for White British People without any kinship
  myTab = copy(myTab_cross) 
  myTab = myTab[ancestry == 1001,]
  myTab = myTab[kinship == 0,]
  
  #' Filter for consent
  ToExclude = fread(gsub("ukb672224.tab","withdraw98032_20250818.csv",UKB_phenotypes))
  table(is.element(myTab$ID,ToExclude$V1))
  myTab = myTab[!is.element(ID,ToExclude$V1),]
  
  #' Filter for same sex information
  myTab[,table(sex,geneticSex)]
  myTab = myTab[sex == geneticSex,]
  
  #' Check dimension after filtering and save
  dim(myTab)
  save(myTab, file = paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_filtered.RData"))
  write.table(myTab$ID,file = paste0(UKB_phenotypes_filtered,"/01_Prep_01_UKB_filtered_samples.txt"), 
              col.names = F, row.names = F, quote = F)
  
}

#' ## Collapse columns
#' 
#' - extract medication information (ATC code C10 - lipid modifying agents, coding from Wu, Y., Byrne, E.M., Zheng, Z. et al. Genome-wide association study of medication-use and associated disease in the UK Biobank. Nat Commun 10, 1891 (2019). https://doi.org/10.1038/s41467-019-09572-5; Supplemental Data 1)
#' - extract CAD information (ICD-10 codes I20-I25, coding table 19 from UKB)
#' 
{
  # load(paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_filtered.RData"))
  dim(myTab)
  
  #' Get medication coding
  codingTable = data.table(read_xlsx(paste0(pathData,"/Wu_2019_SupplementalData1_modified.xlsx"),sheet=1))
  lipidLowering = codingTable[grepl("C10",ATC_Code),Coding]
  
  #' Collapse medication BL
  myMeds = names(myTab)[grep("medication_BL",names(myTab))]
  myTab[,lipidLow_BL := 0]
  for(i in 1:length(myMeds)){
    #i=1
    myTab[get(myMeds[i]) %in% lipidLowering,lipidLow_BL := 1]
  }
  myTab[,get("myMeds"):=NULL]
  myTab[,table(lipidLow_BL,sex)]
  
  #' Collapse medication FU
  myMeds = names(myTab)[grep("medication_FU",names(myTab))]
  myTab[,lipidLow_FU := 0]
  for(i in 1:length(myMeds)){
    #i=1
    myTab[get(myMeds[i]) %in% lipidLowering,lipidLow_FU := 1]
  }
  myTab[,get("myMeds"):=NULL]
  myTab[,table(lipidLow_FU,sex)]
  
  #' Get CAD coding
  codingTable = fread(paste0(pathData,"/UKB_coding19.tsv"))
  codingTable = codingTable[grepl("I20",meaning) | grepl("I21",meaning) | 
                              grepl("I22",meaning) | grepl("I23",meaning) | 
                              grepl("I24",meaning) | grepl("I25",meaning),]
  codingTable = codingTable[selectable == "Y",]
  CAD_codes = codingTable[,coding]
  
  #' Collapse diagnoses
  myDiag1 = names(myTab)[grep("diagnoses_main",names(myTab))]
  myDiag2 = names(myTab)[grep("diagnoses_sec",names(myTab))]
  myDiag3 = names(myTab)[grep("diagnoses_any",names(myTab))]
  myDiag = c(myDiag1,myDiag2,myDiag3)
  myTab[,CAD := 0]
  
  for(i in 1:length(myDiag)){
    #i=1
    myTab[get(myDiag[i]) %in% CAD_codes,CAD := 1]
  }
  myTab[,get("myDiag"):=NULL]
  myTab[,table(CAD,sex)]
  
  #' Check dimension after filtering and save
  dim(myTab)
  save(myTab, file = paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_filtered_collapsed.RData"))
  # load(paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_filtered_collapsed.RData"))
  
}

#' # Get GP data ####
#' ***
#' 
#' ## Total cholesterol ####
#' 
#' Codings from Denaxas S, Shah AD, Mateen BA, Kuan V, Quint JK, Fitzpatrick N, Torralbo A, Fatemifar G, Hemingway H. A semi-supervised approach for rapidly creating clinical biomarker phenotypes in the UK Biobank using different primary care EHR and clinical terminology systems. JAMIA Open. 2020 Dec 5;3(4):545-556. doi: 10.1093/jamiaopen/ooaa047; Supplemental Table 2 - restricted to cholesterol
#' 
{
  #' Load primary care data of TC and filter for IDs in my subset 
  UKB_primaryCare = gsub("phenotypes/ukb672224.tab","primary_care/gp_clinical.txt",UKB_phenotypes)
  myTab_GP_TC = fread(UKB_primaryCare)
  myTab_GP_TC = myTab_GP_TC[eid %in% myTab$ID,]
  myTab_GP_TC[,length(unique(eid))]
  
  #' Load coding table for primary care data (Denaxas et al, Supplementary Table 2)
  codingTable = fread(paste0(pathData,"/Denaxas_2020_SupplementalTable2_modified.txt"))
  codingTable[,table(duplicated(readcode))]
  
  #' Filter for reads from the TC coding
  myTab_GP_TC = myTab_GP_TC[read_3 %in% codingTable$readcode | read_2 %in% codingTable$readcode,]
  myTab_GP_TC[,table(read_3)]
  myTab_GP_TC[,table(read_2)]
  
  #' Check that there is information on ID, date and TC value
  names(myTab_GP_TC)
  myTab_GP_TC[,table(is.na(eid))]
  
  myTab_GP_TC[,table(is.na(event_dt))]
  myTab_GP_TC[,table(event_dt=="")]
  myTab_GP_TC[event_dt==""]
  myTab_GP_TC = myTab_GP_TC[event_dt!=""]
  
  myTab_GP_TC[,table(is.na(value1))]
  myTab_GP_TC[,table(value1=="")]
  myTab_GP_TC[,table(is.na(value2))]
  myTab_GP_TC[,table(value2=="")]
  myTab_GP_TC[,table(is.na(value3))]
  myTab_GP_TC[,table(value3=="")]
  
  myTab_GP_TC[value2 != "",table(value1)]
  myTab_GP_TC[value3 != "",table(value1)]
  myTab_GP_TC[,table(value3=="",value2=="")]
  myTab_GP_TC[(value1 == "" | value1=="OPR003") & value2 != "",value1 := value2]
  myTab_GP_TC[(value1 == "" | value1=="OPR003") & value3 != "",value1 := value3]
  
  myTab_GP_TC = myTab_GP_TC[value1 != ""]
  myTab_GP_TC[,value1 := as.numeric(value1)]
  myTab_GP_TC = myTab_GP_TC[!is.na(value1),]
  myTab_GP_TC[,table(value1==0)]
  myTab_GP_TC = myTab_GP_TC[value1!=0]
  myTab_GP_TC[,table(value1<0)]
  
  #' Check distribution of TC values - similar to Ko et al., I filter everything above 30
  hist(myTab_GP_TC$value1)
  hist(myTab_GP_TC[value1<30,value1])
  myTab_GP_TC = myTab_GP_TC[value1<30,]
  
  #' Now reduce to relevant columns: ID, date, TC value
  myTab_GP_TC = myTab_GP_TC[,c(1,3,6)]
  names(myTab_GP_TC) = c("ID","date","TC")
  
  #' Now I can recreate the time since UKB baseline
  myTab_GP_TC[,day := substr(date,1,2)]
  myTab_GP_TC[,month := substr(date,4,5)]
  myTab_GP_TC[,year := substr(date,7,10)]
  myTab_GP_TC[,date2 := paste(year,month,day,sep = "-")]
  myTab_GP_TC[,date2 := as.Date(date2)]
  matched = match(myTab_GP_TC$ID,myTab$ID)
  table(is.na(matched))
  myTab_GP_TC[,date_UKBBaseline := myTab[matched,date_BL]]
  myTab_GP_TC[,age_UKBBaseline := myTab[matched,age]]
  myTab_GP_TC[,date_UKBBaseline := as.Date(date_UKBBaseline)]
  myTab_GP_TC[,dateDif := date2 - date_UKBBaseline]
  myTab_GP_TC[,dateDifYears := as.numeric(dateDif)/365]
  myTab_GP_TC[,age := age_UKBBaseline + dateDifYears]
  myTab_GP_TC[dateDifYears>0,flag := "after baseline"]
  myTab_GP_TC[dateDifYears<0,flag := "before baseline"]
  myTab_GP_TC[dateDifYears==0,flag := "same as baseline"]
  myTab_GP_TC[,table(flag)]
  
  #' Reorder data 
  names(myTab_GP_TC)
  myTab_GP_TC = myTab_GP_TC[,c(1,7,3,12,13)]
  myTab_GP_TC[,table(age<0)]
  myTab_GP_TC[age<0,table(date2)]
  myTab_GP_TC = myTab_GP_TC[age>0,]
  myTab_GP_TC[,min(age)]
  myTab_GP_TC[,max(age)]
  hist(myTab_GP_TC$age)
  myTab_GP_TC = myTab_GP_TC[age>=30 & age<=80,]
  hist(myTab_GP_TC$age)
  
  #' Check duplicates
  myTab_GP_TC[,dumID := paste(ID,date2,sep="__")]
  myTab_GP_TC[,table(duplicated(dumID))]
  dupIDs = myTab_GP_TC[duplicated(dumID),dumID]
  myTab_GP_TC = myTab_GP_TC[!is.element(dumID,dupIDs),]
  myTab_GP_TC[,dumID := NULL]
  
  #' Check dimension after filtering and save
  dim(myTab_GP_TC)
  save(myTab_GP_TC, file = paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_GP_TC.RData"))
  # load(paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_GP_TC.RData"))
  
}

#' ## Lipid lowering treatment #### 
#' 
#' Lipid-Regulating Drugs: BNF code starting with "02.12"
#' 
{
  #' Load primary care data of lipid lowering medication and filter for IDs in my subset 
  UKB_primaryCare_scripts = gsub("phenotypes/ukb672224.tab","primary_care/gp_scripts.txt",UKB_phenotypes)
  myTab_GP_scripts = fread(UKB_primaryCare_scripts)
  myTab_GP_scripts = myTab_GP_scripts[eid %in% myTab$ID,]
  myTab_GP_scripts[,length(unique(eid))]
  
  #' Filter for reads from the medication coding
  myTab_GP_scripts = myTab_GP_scripts[substr(bnf_code,1,5) == "02.12",]
  
  #' Check that there is information on ID and date
  myTab_GP_scripts[,table(is.na(issue_date))]
  myTab_GP_scripts[,table(issue_date=="")]
  myTab_GP_scripts[issue_date==""]
  myTab_GP_scripts = myTab_GP_scripts[issue_date!=""]
  
  myTab_GP_scripts[,day := substr(issue_date,1,2)]
  myTab_GP_scripts[,month := substr(issue_date,4,5)]
  myTab_GP_scripts[,year := substr(issue_date,7,10)]
  myTab_GP_scripts[,date2 := paste(year,month,day,sep = "-")]
  myTab_GP_scripts[,date2 := as.Date(date2)]
  myTab_GP_scripts[,table(is.na(date2))]
  myTab_GP_scripts[,table(date2=="")]
  myTab_GP_scripts[year %in% c(1901,1902,1903),]
  myTab_GP_scripts[year %in% c(2037),]
  myTab_GP_scripts = myTab_GP_scripts[!grepl("1902",date2)]
  myTab_GP_scripts = myTab_GP_scripts[!grepl("2037",date2)]
  
}

#' ## Merge GP data information ####
#' 
#' I will assume that after the first prescription they are always on treatment
{
  FirstScript = myTab_GP_scripts[,min(date2),by = eid]
  FirstScript = FirstScript[eid %in% myTab_GP_TC$ID,]
  myTab_GP_TC[!is.element(ID,FirstScript$eid),statin := 0]
  
  for(i in 1:dim(FirstScript)[1]){
    #i=1
    myTab_GP_TC[ID == FirstScript[i,eid] & date2<=FirstScript[i,V1],statin := 0]
    myTab_GP_TC[ID == FirstScript[i,eid] & date2>FirstScript[i,V1],statin := 1]
    
  }
  myTab_GP_TC = myTab_GP_TC[,c(1,2,3,4,6,5)]
  names(myTab_GP_TC)
  setnames(myTab_GP_TC,"date2","date")
  myTab_GP_TC[,table(is.na(statin))]
  myTab_GP_TC[,table(statin)]
  myTab_GP_TC[,table(statin,flag)]
  
  #' Check dimension after filtering and save
  dim(myTab_GP_TC)
  save(myTab_GP_TC, file = paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_GP_TC_statins.RData"))
  # load(paste0(UKB_phenotypes_filtered,"01_Prep_01_UKB_GP_TC_statins.RData"))
  
}


#' # Merge baseline and GP data ####
#' ***
#' 
{
  #' Modify BL
  myTab_Baseline = copy(myTab)
  myTab_Baseline = myTab_Baseline[ID %in% myTab_GP_TC$ID,]
  myTab_Baseline = myTab_Baseline[,c(1,3,21,8,23)]
  myTab_Baseline[,flag:="UKB baseline"]
  names(myTab_Baseline) = names(myTab_GP_TC)
  myTab_Baseline[,date := as.Date(date)]
  myTab_Baseline = myTab_Baseline[!is.na(TC),]
  
  #' Modify FU
  myTab_FollowUp = copy(myTab)
  myTab_FollowUp = myTab_FollowUp[ID %in% myTab_GP_TC$ID & !is.na(TC_FU),]
  myTab_FollowUp[,dateDif := date_FU - date_BL]
  myTab_FollowUp[,dateDifYears := as.numeric(dateDif)/365]
  myTab_FollowUp[,age_FU := age + dateDifYears]
  myTab_FollowUp = myTab_FollowUp[,c(1,4,22,28,24)]
  myTab_FollowUp[,flag:="UKB follow-up"]
  names(myTab_FollowUp) = names(myTab_GP_TC)
  myTab_FollowUp[,date := as.Date(date)]
  
  #' Combine the 3 data sets
  myTab_long = rbind(myTab_Baseline,myTab_FollowUp,myTab_GP_TC)
  setorder(myTab_long,ID,date)
  
  #' Add time-fixed covariables (sex and genetic PCs)
  matched = match(myTab_long$ID,myTab$ID)
  table(is.na(matched))
  stopifnot(myTab_long$ID == myTab[matched,ID])
  myTab_long[,sex := myTab[matched,sex]]
  myTab_long[sex==0,sex := 2]
  myTab_long = cbind(myTab_long,myTab[matched,10:19])
  
  #' Add menopausal state as additional group variable
  ID_premenopausal = myTab[menopause_BL == 0 & age<50 & sex==0,ID]
  ID_postmenopausal = myTab[menopause_BL == 1 & age>=60 & sex==0,ID]
  myTab_long[sex==1,group_BL := "men"]
  myTab_long[sex==2,group_BL := "women"]
  myTab_long[sex==2 & ID %in% ID_postmenopausal,group_BL := "postmenopausal women"]
  myTab_long[sex==2 & ID %in% ID_premenopausal,group_BL := "premenopausal women"]
  myTab_long[,table(group_BL)]
  
  #' Check dimension and save
  dim(myTab_long)
  save(myTab_long, file = paste0(UKB_phenotypes_filtered,"01_Prep_01_BL_FU_GP_merged.RData"))
  # load(paste0(UKB_phenotypes_filtered,"01_Prep_01_BL_FU_GP_merged.RData"))
  
}

#' # Checks ####
#' *** 
#' - remove duplicated time points (GP vs BL, GP vs FU)
#' - distance between time points >30 days
#' - number of observations per individual
#' - range of TC values
#' - statin treatment assumed constant after initial prescription
#' 
#' ## Check 1: duplicated time points ####
#' 
{
  myTab_long[,dumID := paste(ID,date,sep="__")]
  myTab_long[,table(duplicated(dumID))]
  dupIDs = myTab_long[duplicated(dumID),dumID]
  myTab_long[dumID %in% dupIDs,]
  filt = is.element(myTab_long$dumID,dupIDs) & myTab_long$flag %in% c("same as baseline","after baseline") 
  table(filt)
  myTab_long = myTab_long[!filt,]
  myTab_long[,table(duplicated(dumID))]
  myTab_long[,dumID := NULL]
  dim(myTab_long)
}

#' ## Check 2: distance between time points ####
#' 
{
  test = myTab_long[,.N,by=ID]
  test_IDs = test[N>2,ID]
  
  dumTab = foreach(i = 1:length(test_IDs))%do%{
    #i=1
    dumTab = copy(myTab_long)
    dumTab = dumTab[ID == test_IDs[i],]
    age = dumTab[1,age]
    dumTab[,ageDistOK := T]
    
    for(j in 2:dim(dumTab)[1]){
      #j=2
      test_age = dumTab[j,age]
      age_dif = test_age - age
      if(age_dif<0.1){
        dumTab[j,ageDistOK := F]
        age = dumTab[j-1,age]
      }else{
        dumTab[j,ageDistOK := T]
        age = dumTab[j,age]
      }
    }
    dumTab[,table(ageDistOK)]
    dumTab = dumTab[ageDistOK==T,]
    dumTab
  }
  myTab_long = rbindlist(dumTab)
}

#' ## Check 3: number of time points ####
#' 
#' In the previous check, I remove time points if they were too close to each other. Hence I might now have again individuals with less then 3 time points. 
#' 
{
  TP_overall = myTab_long[,.N,by="ID"]
  hist(TP_overall$N)
  summary(TP_overall$N)
  myTab_long = myTab_long[ID %in% TP_overall[N>2,ID],]
  TP_overall = myTab_long[,.N,by=c("ID","group_BL")]
  TP_overall[,table(group_BL)]
}

#' ## Check 4: range of TC values ####
{
  hist(myTab_long$TC)
  summary(myTab_long[sex==1,TC])
  summary(myTab_long[sex==2,TC])
  
}

#' ## Check 5: lipid lowering medication ####
#' 
#' After the first prescription, the status of lipid lowering medication should always be 1. I have not yet checked that with the BL and FU info of UKB
#' 
{
  test1 = myTab_long[statin==0,max(date),by=ID]
  test2 = myTab_long[statin==1,min(date),by=ID]
  matched = match(test2$ID,test1$ID)
  test2[,date_noStatin := test1[matched,V1]]
  test2 = test2[!is.na(date_noStatin),]
  test2 = test2[date_noStatin>V1,]
  
  for(i in 1:dim(test2)[1]){
    #i=1
    myTab_long[ID == test2[i,ID] & date<test2[i,V1],statin := 0]
    myTab_long[ID == test2[i,ID] & date>=test2[i,V1],statin := 1]
    
  }
  
  test1 = myTab_long[statin==0,max(date),by=ID]
  test2 = myTab_long[statin==1,min(date),by=ID]
  matched = match(test2$ID,test1$ID)
  test2[,date_noStatin := test1[matched,V1]]
  test2 = test2[!is.na(date_noStatin),]
  test2 = test2[date_noStatin>V1,]
  test2
}

#' ## Check 6: sensitivity data set 1 ####
#' 
#' - no lipid lowering medication
#' - age between 40 and 70
#' - first value per year only
{
  myTab_long2 = copy(myTab_long)
  myTab_long2 = myTab_long2[statin==0,]
  
  #' ceiling of age and filter for 1 value per year
  myTab_long2[,age2:= ceiling(age)]
  myTab_long2[,dumID2 := paste0(ID,"_",age2)]
  myTab_long2 = myTab_long2[!duplicated(dumID2),]
  
  #' filter for age group
  myTab_long2[,min(age2)]
  myTab_long2[,max(age2)]
  myTab_long2 = myTab_long2[age2>=40 & age2<=70,]
  
  #' check number of samples and number of time points
  test = myTab_long2[,.N,by=ID]
  hist(test$N)
  test[,table(N>2)]
  myTab_long2 = myTab_long2[ID %in% test[N>2,ID]]
  
  myTab_long[,dumID := paste(ID,date,sep="__")]
  myTab_long2[,dumID := paste(ID,date,sep="__")]
  table(is.element(myTab_long$dumID, myTab_long2$dumID))
  table(is.element(myTab_long2$dumID, myTab_long$dumID))
  myTab_long[,sens1 := F]
  myTab_long[dumID %in% myTab_long2$dumID,sens1 := T]
  
}

#' ## Check 7: sensitivity data set 2 ####
#' 
#' - all time points after UKB BL
#' 
{
  myTab_long2 = copy(myTab_long)
  myTab_long2 = myTab_long2[flag != "before baseline",]
  
  #' check number of samples and number of time points
  test = myTab_long2[,.N,by=ID]
  hist(test$N)
  test[,table(N>2)]
  myTab_long2 = myTab_long2[ID %in% test[N>2,ID]]
  
  myTab_long2[,dumID := paste(ID,date,sep="__")]
  table(is.element(myTab_long$dumID, myTab_long2$dumID))
  table(is.element(myTab_long2$dumID, myTab_long$dumID))
  myTab_long[,sens2 := F]
  myTab_long[dumID %in% myTab_long2$dumID,sens2 := T]
  
}

#' ## Check 8: sensitivity data set 3 ####
#' 
#' - all time points before UKB BL
#' 
{
  myTab_long2 = copy(myTab_long)
  myTab_long2 = myTab_long2[flag != "after baseline",]
  myTab_long2 = myTab_long2[flag != "UKB FU",]
  
  #' check number of samples and number of time points
  test = myTab_long2[,.N,by=ID]
  hist(test$N)
  test[,table(N>2)]
  myTab_long2 = myTab_long2[ID %in% test[N>2,ID]]
  
  myTab_long2[,dumID := paste(ID,date,sep="__")]
  table(is.element(myTab_long$dumID, myTab_long2$dumID))
  table(is.element(myTab_long2$dumID, myTab_long$dumID))
  myTab_long[,sens3 := F]
  myTab_long[dumID %in% myTab_long2$dumID,sens3 := T]
  
}

#' 
#' # Summary ####
#' ***
#' I want a short summary about sample size, number of time points, sex, statin use, age, and TC. 
#' 
#' In addition, I want trajectory plots for the main and sensitivity data. 
#' 
#' ## Summary table ####
{
  test1 = myTab_long[,.N,by=c("ID","sex")]
  test2 = myTab_long[statin==1,.N,by=ID]
  test3 = myTab_long[statin==1,min(age),by=ID]
  
  tab1a = data.table(setting = "MAIN", 
                     sampleSize = dim(test1)[1],
                     NR_measurements_median = as.numeric(summary(test1$N)[3]),
                     NR_measurements_1stQ = as.numeric(summary(test1$N)[2]),
                     NR_measurements_3rdQ = as.numeric(summary(test1$N)[5]),
                     TC_mean = mean(myTab_long[,TC]),
                     TC_SD = sd(myTab_long[,TC]),
                     TC_BL_mean = mean(myTab_long[flag == "UKB baseline",TC]),
                     TC_BL_SD = sd(myTab_long[flag == "UKB baseline",TC]),
                     Female_absolute = dim(test1[sex==2,])[1],
                     Female_percent = dim(test1[sex==2,])[1]/dim(test1)[1],
                     Age_mean = mean(myTab_long[,age]),
                     Age_SD = sd(myTab_long[,age]),
                     Age_BL_mean = mean(myTab_long[flag == "UKB baseline",age]),
                     Age_BL_SD = sd(myTab_long[flag == "UKB baseline",age]),
                     Statin_absolut = dim(test2)[1],
                     Statin_percent = dim(test2)[1]/dim(test1)[1],
                     Age_1st_statin_mean = mean(test3$V1),
                     Age_1st_statin_SD = sd(test3$V1) )
  
  test1 = myTab_long[sens1==T,.N,by=c("ID","sex")]
  
  tab1b = data.table(setting = "SENSITIVITY - no statin", 
                     sampleSize = dim(test1)[1],
                     NR_measurements_median = as.numeric(summary(test1$N)[3]),
                     NR_measurements_1stQ = as.numeric(summary(test1$N)[2]),
                     NR_measurements_3rdQ = as.numeric(summary(test1$N)[5]),
                     TC_mean = mean(myTab_long[sens1==T,TC]),
                     TC_SD = sd(myTab_long[sens1==T,TC]),
                     TC_BL_mean = mean(myTab_long[sens1==T & flag == "UKB baseline",TC]),
                     TC_BL_SD = sd(myTab_long[sens1==T & flag == "UKB baseline",TC]),
                     Female_absolute = dim(test1[sex==2,])[1],
                     Female_percent = dim(test1[sex==2,])[1]/dim(test1)[1],
                     Age_mean = mean(myTab_long[sens1==T,age]),
                     Age_SD = sd(myTab_long[sens1==T,age]),
                     Age_BL_mean = mean(myTab_long[sens1==T & flag == "UKB baseline",age]),
                     Age_BL_SD = sd(myTab_long[sens1==T & flag == "UKB baseline",age]),
                     Statin_absolut = NA,
                     Statin_percent = NA,
                     Age_1st_statin_mean = NA,
                     Age_1st_statin_SD = NA)
  
  test1 = myTab_long[sens2==T,.N,by=c("ID","sex")]
  test2 = myTab_long[sens2==T & statin==1,.N,by=ID]
  test3 = myTab_long[sens2==T & statin==1,min(age),by=ID]
  
  tab1c = data.table(setting = "SENSITIVITY - observations after BL", 
                     sampleSize = dim(test1)[1],
                     NR_measurements_median = as.numeric(summary(test1$N)[3]),
                     NR_measurements_1stQ = as.numeric(summary(test1$N)[2]),
                     NR_measurements_3rdQ = as.numeric(summary(test1$N)[5]),
                     TC_mean = mean(myTab_long[sens2==T,TC]),
                     TC_SD = sd(myTab_long[sens2==T,TC]),
                     TC_BL_mean = mean(myTab_long[sens2==T & flag == "UKB baseline",TC]),
                     TC_BL_SD = sd(myTab_long[sens2==T & flag == "UKB baseline",TC]),
                     Female_absolute = dim(test1[sex==2,])[1],
                     Female_percent = dim(test1[sex==2,])[1]/dim(test1)[1],
                     Age_mean = mean(myTab_long[sens2==T,age]),
                     Age_SD = sd(myTab_long[sens2==T,age]),
                     Age_BL_mean = mean(myTab_long[sens2==T & flag == "UKB baseline",age]),
                     Age_BL_SD = sd(myTab_long[sens2==T & flag == "UKB baseline",age]),
                     Statin_absolut = dim(test2)[1],
                     Statin_percent = dim(test2)[1]/dim(test1)[1],
                     Age_1st_statin_mean = mean(test3$V1),
                     Age_1st_statin_SD = sd(test3$V1) )
  
  test1 = myTab_long[sens3==T,.N,by=c("ID","sex")]
  test2 = myTab_long[sens3==T & statin==1,.N,by=ID]
  test3 = myTab_long[sens3==T & statin==1,min(age),by=ID]
  
  tab1d = data.table(setting = "SENSITIVITY - observations after BL", 
                     sampleSize = dim(test1)[1],
                     NR_measurements_median = as.numeric(summary(test1$N)[3]),
                     NR_measurements_1stQ = as.numeric(summary(test1$N)[2]),
                     NR_measurements_3rdQ = as.numeric(summary(test1$N)[5]),
                     TC_mean = mean(myTab_long[sens3==T,TC]),
                     TC_SD = sd(myTab_long[sens3==T,TC]),
                     TC_BL_mean = mean(myTab_long[sens3==T & flag == "UKB baseline",TC]),
                     TC_BL_SD = sd(myTab_long[sens3==T & flag == "UKB baseline",TC]),
                     Female_absolute = dim(test1[sex==2,])[1],
                     Female_percent = dim(test1[sex==2,])[1]/dim(test1)[1],
                     Age_mean = mean(myTab_long[sens3==T,age]),
                     Age_SD = sd(myTab_long[sens3==T,age]),
                     Age_BL_mean = mean(myTab_long[sens3==T & flag == "UKB baseline",age]),
                     Age_BL_SD = sd(myTab_long[sens3==T & flag == "UKB baseline",age]),
                     Statin_absolut = dim(test2)[1],
                     Statin_percent = dim(test2)[1]/dim(test1)[1],
                     Age_1st_statin_mean = mean(test3$V1),
                     Age_1st_statin_SD = sd(test3$V1) )
  
  
  tab1 = rbind(tab1a,tab1b,tab1c,tab1d)
  tab1
  write.table(tab1,file = "../results/01_Prep_01_summary.txt")
}

#' ## Trajectory plot ####
{
  myIDs = sample(myTab_long$ID,15)
  
  ggplot1 = ggplot(myTab_long[ID %in% myIDs], aes(x=age, y=TC, group=ID, colour=as.factor(ID), shape=as.factor(statin))) +
    geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
    geom_line() + 
    geom_point()+
    labs(x="Age (in years)",
         y="Total cholesterol (in mmol/l)", 
         color="IDs",
         shape="Statin\ntreatment",
         title = paste0("Trajectories of TC levels (example with 15 samples)")) +
    scale_shape_manual(values=c(21,24),
                       labels = c("no","yes"))+
    # scale_colour_manual(values = c("darkred","steelblue"),
    #                     labels = c("women","men"))+
    theme_classic() + 
    guides(color = "none")
  ggplot1
  
  filename = paste0("../results/_figures/01_Trajectory/MAIN_setting_15samples.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(ggplot1)
  dev.off()
  
  myIDs = sample(myTab_long[sens1==T,ID],15)
  
  ggplot2 = ggplot(myTab_long[ID %in% myIDs & sens1==T], aes(x=age, y=TC, group=ID, colour=as.factor(ID), shape=as.factor(statin))) +
    geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
    geom_line() + 
    geom_point()+
    labs(x="Age (in years)",
         y="Total cholesterol (in mmol/l)", 
         color="IDs",
         shape="Statin\ntreatment",
         title = paste0("Trajectories of TC levels (sensitivity setting, example with 15 samples)")) +
    scale_shape_manual(values=c(21,24),
                       labels = c("no","yes"))+
    # scale_colour_manual(values = c("darkred","steelblue"),
    #                     labels = c("women","men"))+
    theme_classic() + 
    guides(color = "none")
  ggplot2
  
  filename = paste0("../results/_figures/01_Trajectory/SENSITIVITY_setting_15samples.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(ggplot2)
  dev.off()
}

#' ## Boxplots per sex ####
{
  plotData = copy(myTab_long)
  plotData = plotData[sens1==T,]
  
  plotData[,age2:= ceiling(age)]
  plotData[,dumID := paste0(ID,"_",age2)]
  plotData = plotData[!duplicated(dumID),]
  
  plotData[,min(age2)]
  plotData[,max(age2)]
  
  plotData[age2<=45, age3 := "40-45"]
  plotData[age2<=50 & age2>45, age3 := "46-50"]
  plotData[age2<=55 & age2>50, age3 := "51-55"]
  plotData[age2<=60 & age2>55, age3 := "56-60"]
  plotData[age2<=65 & age2>60, age3 := "61-65"]
  plotData[age2<=70 & age2>65, age3 := "66-70"]
  plotData[age2>70, age3 := "71-79"]
  
  plotData[sex==2,sex2 := "females"]
  plotData[sex==1,sex2 := "males"]
  
  ggplot3 = ggplot(plotData, aes(x=as.factor(age3), y=TC, fill=sex2)) +
    geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
    geom_boxplot(position=position_dodge(1)) +
    scale_fill_manual(values=c("darkred","steelblue"),
                      labels = c("women","men"))+
    ylim(2, 12) +
    labs(x="Age (in years)",
         y="Total cholesterol (in mmol/l)", 
         fill="Sex",
         title = paste0("Boxplots of TC levels")) +
    theme_classic() 
  ggplot3
  
  filename = paste0("../results/_figures/01_Boxplots/SENSITIVITY_setting_ageCategories.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(ggplot3)
  dev.off()
  
  ggplot4 = ggplot(plotData, aes(x=as.factor(age2), y=TC, fill=sex2)) +
    geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
    geom_boxplot(position=position_dodge(1), show.legend = FALSE) +
    scale_fill_manual(values=c("darkred","steelblue"),
                      labels = c("women","men"))+
    ylim(2, 12) +
    labs(x="Age (in years)",
         y="Total cholesterol (in mmol/l)", 
         fill="Sex",
         title = paste0("Boxplots of TC levels")) +
    theme_classic() 
  ggplot4
  filename = paste0("../results/_figures/01_Boxplots/SENSITIVITY_setting_agePerYear.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(ggplot4)
  dev.off()
  
  plotData = plotData[group_BL != "women",]
  plotData[grepl("postmenopausal",group_BL) & flag == "before baseline", TC := NA]
  plotData[grepl("premenopausal",group_BL) & flag == "after baseline", TC := NA]
  plotData = plotData[!is.na(TC),]
  
  ggplot5 = ggplot(plotData, aes(x=as.factor(age3), y=TC, fill=group_BL)) +
    geom_hline(yintercept = 8,linetype="dotted", linewidth=0.5)+
    geom_boxplot(position=position_dodge(1)) +
    scale_fill_manual(values=c("steelblue","darkred","salmon"),
                      labels = c("men","postmenopausal\nwomen","premenopausal\nwomen"))+
    ylim(2, 12) +
    labs(x="Age (in years)",
         y="Total cholesterol (in mmol/l)", 
         fill="Sex and \nmenopausal status",
         title = paste0("Boxplots of TC levels")) +
    theme_classic() 
  ggplot5
  filename = paste0("../results/_figures/01_Boxplots/SENSITIVITY_setting_ageCategories_menopause.png")
  png(filename = filename,width = 1800, height = 1000, res=200)
  plot(ggplot5)
  dev.off()
  
  
}

#' # Save data ####
#' ***
myTab_long
myTab_long[,ageDistOK := NULL]
myTab_long[,dumID := NULL]

save(myTab_long, file = paste0(UKB_phenotypes_filtered,"01_Prep_01_BL_FU_GP_merged_filtered.RData"))
myTab_IDs = myTab_long[,.N,by=ID]
write.table(myTab_IDs$ID,file = paste0(UKB_phenotypes_filtered,"/01_Prep_01_BL_FU_GP_merged_filtered_samples.txt"), 
            col.names = F, row.names = F, quote = F)

#' # SessionInfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

