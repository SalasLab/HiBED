# HiBED
### Publication DOI: [10.3389/fnins.2023.1198243](https://www.frontiersin.org/articles/10.3389/fnins.2023.1198243/full)
### Bioconductor Package: https://bioconductor.org/packages/devel/data/experiment/html/HiBED.html
Hierarchical deconvolution for extensive cell type resolution in the human brain using DNA methylation.
The HiBED deconvolution estimates proportions up to 7 cell types (GABAergic neurons, glutamatergic neurons, astrocytes, microglial cells, oligodendrocytes, endothelial cells, and stromal cells) in bulk brain tissues.

![Figure1](https://user-images.githubusercontent.com/32206453/224516354-75e2b4bd-102f-4c11-be84-e40f36daf5f0.png)

## Installation
```
devtools::install_github("SalasLab/HiBED")
```

## Load package
```
library(HiBED)
```

## Deconvolution function
```
?HiBED_deconvolution
```

## Example
```
library(FlowSorted.Blood.EPIC)
library(FlowSorted.DLPFC.450k)
library(minfi)
Mset<-preprocessRaw(FlowSorted.DLPFC.450k)
Examples_Betas<-getBeta(Mset)
HiBED_result<-HiBED_deconvolution(Examples_Betas, h=2)
head(HiBED_result)
```

## EPICv2 Solution
#### if you're having trouble dowloading v2 annotation files, you can find them at 
####  https://www.dropbox.com/sh/rbxjhq9zalqq58e/AAABR8kKegXKVMeNJV8a2lJRa?dl=0
```
library(HiBED)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
```
#### dir should be the folder containing EPICv2 IDATs
```
v2_RGset = read.metharray.exp("dir",recursive = TRUE) 
annotation(v2_RGset)["array"] = "IlluminaHumanMethylationEPICv2" #Update annotation files for v2
annotation(v2_RGset)["annotation"] = "20a1.hg38"
v2_MSet <-preprocessNoob(v2_RGset)
v2_Betas<-getBeta(v2_MSet)
v2_Betas<- sesame::betasCollapseToPfx(v2_Betas)
HiBED_deconvolution(v2_Betas, h=2)
```


