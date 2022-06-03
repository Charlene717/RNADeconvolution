## CellR
## Ref: https://github.com/adoostparast/CellR

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Installation #####
# devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))
# install.packages("devtools")

devtools::install_github("adoostparast/CellR")

library(CellR)

## Ref:https://www.geeksforgeeks.org/working-with-different-versions-of-an-r-package/
library(versions)
devtools::install_version(package ='Seurat',  version = package_version('2.3.0'),lib="./SeuratV230/")


##### A real example to run #####
Bulk <- read.table(paste0("./Demo_Data/CellR_Sample_Data/Bulk.txt"),header=TRUE)

Single <- read.table(paste0("./Demo_Data/CellR_Sample_Data/FCortex.txt"),header=TRUE)

GTEx <- read.table(paste0("./Demo_Data/CellR_Sample_Data/GTExExpressionSymbol.txt"),header=TRUE)

Cells <- read.table(paste0("./Demo_Data/CellR_Sample_Data/Cells.txt"),header=TRUE)

# Output1 <- CellR::Deconvolution(Bulk,Single,GTEx,Cells,3,200,12,1,2500,1)
Output1 <- Deconvolution(Bulk,Single,GTEx,Cells,3,200,12,1,2500,1)

## Ch
trace(Deconvolution,edit = TRUE)


##### Estimating cell-specific gene expression profiles using CellR #####
Data <- read.table("./Demo_Data/CellR_Sample_Data/Bulk_data.txt",header=TRUE)

Proportion <- read.table("./Demo_Data/CellR_Sample_Data/Proportions.txt", header=TRUE)

Frequency <- read.table("./Demo_Data/CellR_Sample_Data/Frequency.txt", header=TRUE)

Reference <- read.table("./Demo_Data/CellR_Sample_Data/Reference.txt",header=TRUE)

Output2 <- Expression_estimate(Data, Proportion, Frequency, Reference, 1)



##### Export RData #####
save.image("D:/Dropbox/##_GitHub/##_Charlene/RNADeconvolution/CellR_Example.RData")
