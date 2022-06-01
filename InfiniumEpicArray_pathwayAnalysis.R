#BiocManager::install("missMethyl") # install the package
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

rm(list=ls())
gc()

# Load summary statistics from EWAS and order based on p-values
beta<-read.csv("mother_results.csv",as.is = T)
beta<-beta[order(beta$P.value),]
rownames(beta)<-beta$X
beta<-beta[-1]

## Select CpGs based on your criteria, choose one option
topCpGs<-subset(beta,FDR<0.05) ## option 1: select CpGs with fdr < 0.05 (Works if you have > 100 CpGs)
#topCpGs<-subset(beta,fdr<0.2) ## option 2: select CpGs with fdr < 0.2 (If fdr < 0.05 is strict, that's another alternative)
#topCpGs<-beta[1:500,] ## option 3: select top 100 CpGs (Can also try 500 and 1000)
#topCpGs<-subset(beta,p<0.05) ## option 4: select CpGs with p < 0.05 (If no FDR significant CpGs, that's another alternative)

dim(topCpGs)
#[1] 760  11

## Get the names of selected CpGs
sigCpGs<-topCpGs$CPG.Labels

## GO pathway
gst<-gometh(sig.cpg = sigCpGs, all.cpg = beta$CpGs, collection = "GO", array.type = "EPIC", prior.prob = T, plot.bias = F)
gsa<-subset(gst,FDR<0.05) ## FDR significant pathways
write.csv(gst,file = "Mother_GOanalysis_results.R")

## KEGG
gst.kegg<-gometh(sig.cpg = sigCpGs, all.cpg = beta$CpGs, collection = "KEGG", array.type = "EPIC", prior.prob = T, plot.bia=F)
gsa.kegg<-subset(gst.kegg,FDR<0)
write.csv(gst.kegg,file = "Mother_KEGGanalysis_results.R")

