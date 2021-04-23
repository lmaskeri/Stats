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
library(dplyr)


## Load data
data <- read.csv('/homes/anovak9/stats/data.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(data)[1] <- 'sample'
str(data)
summary(data)
data[1:5,1:5]

#Data has 801 samples, 20531 genes
dim(data)
#Store number of zero counts per gene
count <- function(x){
  sum(x == 0)
}
zeros <- apply(data[2:20532], MARGIN = 2, FUN = count)
zeros <- as.data.frame(zeros)
zeros <- as.numeric(zeros$zeros)
summary(zeros)
hist(zeros)
zeros <- zeros %>% arrange(zeros)
#keep 75% after rank 

## Normalization