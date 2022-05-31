# EpicArrayAnalysis
R based Illumina Infinium Epicarray analysis
## Introduction

Illuminas Infinium MethylationEPIC arrays is a genome wide methylation analysis method which quantitaively detects methylation of over 850K methylation sites at single nucleotide resolution. The sites under investigation encompass CpG islands, open chromatin sites, enhancer regions and other regulatory regions of genome.

During the processing of the DNA sample under investigation, DNA is treated with bisulphite where unmethylated Cytosines, e.g 5mC and 5hmC, get converted to Uracil while methylated Cytosines remain unchanged.  In the array technology site specific probes are used that can differentiate between methylated and unmethylated versions and final read out is relative fluorescent signals from the methylated vs. unmethylated sites. The outline of the process is illustrated below to demonstrate the concept. 

![Infinium array Process Image](./img/image1.png)

It is crucial to understand the layout of the arrays as well as this information can help to identify and remove chip based confounding factors. The layout of the arrays may differ between what is showned here vs what you might be using, so be aware of it.   Each chip has a chipID number and there are 8 slots (here) that allows assay for 8 samples.  Each slot is labelled as R0NC01, where 0N is the row number, since there is only one column so it is always C01.

### Data Format:
The scanner will produce .jpg image and .idat files.  IDAT is Illumina's properietary format to store summary intensities for each probe.  So for analysis we will be using IDAT files. The files are named as
```
ChipID_Numbe_SlotID_col.idat
204792210063_R01C01_Grn.idat  -> Green Channel Intensities
204792210063_R01C01_Red.idat  -> Red Channel Intensities
```
### Note

In order to do a complete analysis we will also need metadata of each sample. While processing samples for a case-control study of few laboratory samples, one may not have a metadata, thats fine, but use whatever information to its maximum.  However in my case I have 192 human samples and I would like to estimate the changes in methylation so I am collecting as much metadata as I can including `Gender, Age, Race, smokingStatus` other relevant information.

Cell type Deconvolution :  If the samples have cells of mixed type e.g blood or saliva samples,  one would like to apply deconvolution to account for cell types. My samples were Saliva samples so for EWAS I will be controling for Epithelial, Fibroblast, B, NK, CD4T, CD8T, Monocytes. 

Ancestry(race): For methylation principal components, Comp.2 and Comp.3 are more accurate in assessing ancestry. One can also test the correlation between self-reported race and mPCs to see which components are better predictors. One can also use a PC plot to visualize.

## Analysis

Step1: Quality Check of the probes

This step primarly collate all the datasets together, and then apply filters to remove bad samples and bad probes.

a. The required packages are loaded and missing are installed installed using the `install_needed_packages.R` script.
b. 

