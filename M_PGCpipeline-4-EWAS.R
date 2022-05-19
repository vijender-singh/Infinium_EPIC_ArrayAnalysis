#########################################################################
### PGC Pipeline - Script 4 - EWAS
#########################################################################
#' Clean
rm(list=ls())
gc()

#' Load packages
load_pkgs <- c("CpGassoc", "data.table", "tibble", "feather")
lapply(load_pkgs, require, character.only = TRUE)

## Load Methylation Data
beta<-fread("mother_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_ptsd.csv", data.table = F)
rownames(beta)<-beta$V1
beta<-beta[,-1]

## Load phenotype file with MethylationID (combination of Sentrix ID and Sentrix Position) as the 1st column
pheno <- read.csv("mother_Pheno_QC_04_mPC_Cell.csv", row.names = 1)

cpg <- beta[, colnames(beta) %in% row.names(pheno)]
cpg <- cpg[, order(colnames(cpg))]
pheno <- pheno[rownames(pheno) %in% colnames(cpg), ]
pheno <- pheno[order(rownames(pheno)), ]
table(rownames(pheno) == colnames(cpg)) # should be TRUE

#pheno$PTSD<-as.character(pheno$Sample_Group)

# pheno$PTSD[pheno$PTSD=="PTSD_M"] <- 1
# pheno$PTSD[pheno$PTSD=="NOPTSD_M"] <- 0


# Saving the clean phenotype file
# This phenotype file includes only the samples that passed the QC and tested in EWAS
write.csv(pheno, file = "mother_Pheno_QC_EWAS.csv",row.names = T)

## Define variables
study <- "Mother_PTSD" # name of the study, e.g. "GTP", "DNHS" etc.
ptsdVar <- "PTSD" # name of the ptsd variable, coded as: cases = 1 and controls = 0
ptsd <- pheno[, ptsdVar, FALSE]

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
covar<-data.frame(pheno[,c("Epi","Fib","B","NK","CD4T","CD8T","Mono","Mage")])

## Run EWAS with CpGAssoc
#test <- cpg.assoc(cpg, pheno[, ptsdVar], covar, logit.transform = T, large.data=TRUE)
test <- cpg.assoc(cpg, pheno[, ptsdVar], logit.transform = T, large.data=TRUE)
assoc <- test$results
eff <- test$coefficients
results <- cbind(assoc,eff)
write.csv(results,file="mother_results.csv")
save(results, file = paste(study,"mother_EWAS.Rdata", sep = "_")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP
signfResults<-results[results$FDR < 0.05,]
write.csv(signfResults,file="mother_results_-lt_FDR_05.csv")
