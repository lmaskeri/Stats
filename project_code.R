## Load necessary packages
#Install EDASeq package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EDASeq")
#Load package
library(EDASeq)
#Install glmnet
install.packages('glmnet')
library(glmnet)

## Load data
data <- read.csv('/homes/anovak9/stats/data.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(data)[1] <- 'sample'
str(data)


#Need to figure out if they used RPM, RPKM, FPKM, or TPM to quantify counts of a gene



## Normalization