# HiBED
Hierarchical deconvolution for extensive cell type resolution in the human brain using DNA methylation
The HiBED deconvolution estimates proportions up to 7 cell types (GABAergic neurons, glutamatergic neurons, astrocytes, microglial cells, oligodendrocytes, endothelial cells, and stromal cells) in bulk brain tissues

![Figure1](https://user-images.githubusercontent.com/32206453/224516354-75e2b4bd-102f-4c11-be84-e40f36daf5f0.png)

## Installation
```
devtools::install_github("SalasLab/HiBED")
```

## Load library 
```
library(HiBMED)

load("HiBED_Libraries")
```

## Deconvolution function
```
?HiBED_deconvolution
