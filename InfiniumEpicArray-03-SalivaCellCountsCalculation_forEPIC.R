BiocManager::install("EpiDISH") # install the package

# Load package and the reference panel
library(EpiDISH)
data(centEpiFibIC.m) # loading this panel takes some time...
data(centBloodSub.m) # loading blood panel, also takes some time...

# Load your data
#beta <- read.csv("") # rownames should be cpgs and colnames should be subjects

beta <- fread("mother_noob_qcd_crossReactiveProbesRemoved.csv", data.table = F)
rownames(beta) <- beta$V1
beta <- beta[,-1]

frac.m <- hepidish(beta.m = beta, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')

# Output matrix (frac.m) contains epithelial (Epi), fibroblast (Fib), B, NK, CD4T, CD8T, Mono, Neutro, Eosino
# Then for EWAS you can control for Epi, Fib, B, NK, CD4T, CD8T, Mono
phen <- read.csv("mother_Pheno_QC_03_mPC.csv",stringsAsFactors=FALSE,header=TRUE)
phen <- merge(phen,frac.m,by.x = 1,by.y = "row.names",all.x = T)
write.csv(phen, file = "mother_Pheno_QC_04_mPC_Cell.csv",row.names = F) # Save the phenotype file with mPCs

