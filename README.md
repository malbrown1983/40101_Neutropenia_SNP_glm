# 40101_Neutropenia_SNP_glm
glm of all the SNPs in the 40101 neutropenia data

setwd("~/Documents/PostDoc/Neutropenia_Data/Gen_Phen_files+Script/")

library(GenABEL)
library(gdata)
library(survival)
library(ggplot2)
library(data.table)

options(width=70)

geno.dat = "~/Desktop/plink/40101/40101-907QCv2.raw"
tools::md5sum(geno.dat)

##fread is like read.table but faster for larger data sets and the file is about a gb
genoAll <- fread(geno.dat, header = TRUE, data.table=FALSE)

##The first checks if the first character is numeric and, if so, adds an “a” at the beginning.  
##The second goes through and replaces all dashes with an underscore.  After you do those, 
##redefine snplist1
colnames(genoAll)[grepl("^[0-9]", colnames(genoAll))] = paste("a", colnames(genoAll)[grepl("^[0-9]", colnames(genoAll))], sep="")
colnames(genoAll) = gsub("-", "_", colnames(genoAll))

phenofile="ac_arm_anc_40101jan32014_fmEdits.csv"
tools::md5sum(phenofile)

pheno=read.csv(phenofile)
table(pheno$anc_event,exclude=NULL)

##CALGB40101-AC-sampleinfo.csv is a post QC generated csv file that has validate, genetic_cauc
##and consent columns. 1=yes 0=no I assume the validate column is patient information passing
##QC. This allows us to combine the dataframe with a phenofile and avoid loading a genofile.
##Now post QC and genetic Europeans can be subsetted with the phenofile info.
pdat=read.csv("CALGB40101-AC-sampleinfo.csv")
pdat=pdat[,c("expid", "patid","validate","consent","genetic_cauc")]
pdat=pdat[pdat$validate==1,]

##Subsetting cohort data

###value ancevent 
####   1= 'Neutropenic Event' 
####   0= 'Control' 
####  21='Inadequate AE Documentation' 
####  22='Less than 4 cycles of treatment' 
####  23='Only Drug given on 1st cycle' 
####  24='Wrong Drug Given' 
####  25='No treatment drug information' 
####  26='Concurrent non-protocol chemotherapy';

alldat=merge(pheno[,c("RACE_ID","patid","race","ethnic","racethnic","anc_event","event_cycle","ANC","wbc","wbcflag","bili","bilirubin","cc","GFR","G_CSF_cyc1","prdi","age_reg","ER_PgR","PERFORMANCE_ID","nodeneg2","tumor2","HERNeg","BMI","BMI_Cat","X_Con_NP_bio_resp_mod","relative_dose","TT033","OH018")],pdat,all.x=TRUE)
dim(alldat)

##All post QC patients
validdat=alldat[!is.na(alldat$validate)&alldat$validate==1&!is.na(alldat$consent)&alldat$consent==1,]
dim(validdat)

##All genetic Europeans
geneticdat=validdat[validdat$genetic_cauc==1,]
dim(geneticdat)


##The following is for the other races
blackdat=validdat[validdat$RACE_ID==3,]
dim(blackdat)

asiandat=validdat[validdat$RACE_ID==4,]
dim(asiandat)

pacifdat=validdat[validdat$RACE_ID==5,]
dim(pacifdat)

nativedat= validdat[validdat$RACE_ID==6,]
dim(nativedat)


##For all 537848 SNPs
rownames(geneticdat)<- geneticdat$expid
rownames(genoAll)<- genoAll$FID

geneticdat$ANCevent=ifelse(geneticdat$anc_event==1,1,ifelse(geneticdat$anc_event==0,0,NA))
geneticdat$AGE=ifelse(geneticdat$age_reg<65,0,1)
geneticdatgenoAll= merge(genoAll, geneticdat, by=0)
### merge by row names (by=0 or by="row.names")
###Since we can not list all 575575 separately in the formula we read in each SNP from geno file 
###with a function

##Had to run these commands first before the function worked
###rsid= snplist1
###sdat=validdatgenoAll

##for genetic Europeans
gwas.analysis<-function(rsid,sdat){
    #sdat = geneticdatgenoAll[,c("ANCevent", rsid)]
	fmod0<-as.formula(paste('ANCevent~',rsid))
	mod0<-summary(glm(fmod0,data=sdat, family="binomial"))

	if(nrow(mod0$coefficients) > 1){	# check if coefficient estimated for rsid, otherwise, skip to else and set pval = NA

        pval = coef(mod0)[rsid,"Pr(>|z|)"]
        
    }
    else{	# if all genotypes the same, set p-value to NA
        pval = NA
    }
    
    out = c(rsid, pval)
    return(out)
}


snplist1 = colnames(genoAll)[7:ncol(genoAll)]

output1 = t(sapply(snplist1[1:100], function(id) gwas.analysis(id,geneticdatgenoAll[,c("ANCevent", id)])))

write.table(output1,"~/Desktop/40101_Model_Results/SNP_mod_CEU_results1-10000.xls",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

save(output1, file=paste("40101Cglm_757_537848.RData",sep=""))

