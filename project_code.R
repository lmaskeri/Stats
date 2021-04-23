## Load necessary packages
install.packages('glmnet')
library(glmnet)
library(dplyr)

## Load data
data <- read.csv('/homes/anovak9/stats/data.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(data)[1] <- 'sample'
str(data)
summary(data)


## Pre-processing/Normalization

#Data has 801 samples, 20531 genes
dim(data)
#Store number of zero counts per gene
count <- function(x){
  sum(x == 0)
}
zeros <- apply(data[2:20532], MARGIN = 2, FUN = count)
#Format zero count DF
zeros <- as.numeric(zeros)
zeros <- as.data.frame(zeros)
genes <- colnames(data)[2:20532]
rownames(zeros) <- genes
colnames(zeros)[1] <- 'Count'
#Look at info about zero counts
summary(zeros)
hist(zeros$Count)
#Arrange DF so counts in ascending order
zeros <- zeros %>% arrange(Count)
#Keep top 75% of rows in zeros
keep <- as.integer(20531*0.75)
keep
zeros2 <- zeros %>% slice(1:keep)
dim(zeros2)
genes_to_keep <- rownames(zeros2)
#Subset data DF based on genes left in zeros DF
data2 <- data[,colnames(data) %in% genes_to_keep]
dim(data2)
#Apply log2 transform 
data3 <- apply(data2, MARGIN = 2, FUN = log2)
#Replace all -Inf values (aka missing values) with -1
length(which(data3 == -Inf))
data3[data3 == -Inf] <- -1
length(which(data3 == -Inf))
length(which(data3 == -1))

## Ordination



## Training



## Testing



## Evaluations








