#### Load necessary packages#### 
# install.packages('glmnet')
# install.packages("umap") # click no to compiling from source
# install.packages("Rtsne")
# install.packages("caret")
library(glmnet)
library(dplyr)
library(umap)
library(Rtsne)
library(ggplot2)
library(caret)

## Load data
# data <- read.csv('/homes/anovak9/stats/data.csv', header = TRUE, stringsAsFactors = FALSE)
data <- read.csv('C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Loyola/Classes/STAT 437/Final-Project/TCGA-PANCAN-HiSeq-801x20531/data.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(data)[1] <- 'sample'
label <- read.csv('C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Loyola/Classes/STAT 437/Final-Project/TCGA-PANCAN-HiSeq-801x20531/labels.csv') #stringsAsFactors = FALSE
colnames(label)[1] <- 'sample'
# str(data)
# summary(data)


#### Pre-processing/Normalization #### 

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
#Replace all -Inf values (aka values that were originally 0) with -1 
length(which(data3 == -Inf))
data3[data3 == -Inf] <- -1
length(which(data3 == -Inf))
length(which(data3 == -1))
#Transform back into a dataframe
data3 <- data.frame(data3)
#Add the sample column back
data4 <- cbind(data[,'sample'],data3)
colnames(data4)[1] <- 'sample'
#Bind the tumor class
data4 <- inner_join(label,data4)


#### Ordination #### 

## PCA

PCA <- prcomp(data4[,-c(1:2)],center = TRUE, scale. = TRUE)
length(PCA$sdev) # I noticed 
# Select the variables that capture  ~0.95 of the variance
V.explained <- cumsum(PCA$sdev^2/sum(PCA$sdev^2))
PCA_data <- data.frame(PCA$x[,which(V.explained <= 0.95)])
# Add sample data back to the PCA data
PCA_data <- cbind(label,PCA_data)

#convex hull for plot
hull_PCA <- PCA_data %>% group_by(Class) %>% slice(chull(PC1,PC2))
# Create a plot using the first two PCs
# PCA_plot <- autoplot(PCA, data = data4, colour = "Class", frame = TRUE)
PCA_plot <- ggplot(data = PCA_data, aes(x = PC1,y = PC2)) + aes(color = Class, fill = Class)+ geom_point() + geom_polygon(data = hull_PCA, alpha = 0.4)


## UMAP - can change the number of resulting dimensions, the default is 2

UMAP <- umap(data4[,-c(1:2)])
UMAP_data <- data.frame(cbind(label,UMAP$layout))

#convex hull for plot
hull_UMAP <- UMAP_data %>% group_by(Class) %>% slice(chull(X1,X2))
# Create a plot using the first two Axis
UMAP_plot <- ggplot(data = UMAP_data, aes(x = X1,y = X2)) + aes(color = Class, fill = Class)+ geom_point() + geom_polygon(data = hull_UMAP, alpha = 0.4)


## tSNE - can change the number of resulting dimensions, the default is 2

# theta = 0.0 is an exact tSNE
tSNE <- Rtsne(as.matrix(data4[,-c(1:2)]), theta=0.0)
tSNE_data <- cbind(label,tSNE$Y)

#convex hull for plot
hull_tSNE <- tSNE_data %>% group_by(Class) %>% slice(chull(`1`,`2`))
# Create a plot using the first two Axis
tSNE_plot <- ggplot(data = tSNE_data, aes(x = `1`,y = `2`)) + aes(color = Class, fill = Class) + geom_point()+ geom_polygon(data = hull_tSNE, alpha = 0.4)


#### Train Test Splitting ####

# The paper did a (Train:108,Test:77) split (line 90, Normalization.R)

set.seed(134)
# Train test split, current ratio is 9:1 train:test
Train <- createDataPartition(data4$Class, p=0.9, list=FALSE)


#### Model fitting and cross validation within the test set ####

# 10-fold cross validation
cv_data4 <- cv.glmnet(as.matrix(data4[Train,-c(1,2)]), factor(data4[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                      type.multinomial="grouped", nfolds=10)
cv_PCA <- cv.glmnet(as.matrix(PCA_data[Train,-c(1,2)]), factor(PCA_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                    lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                    type.multinomial="grouped", nfolds=10)
cv_UMAP <- cv.glmnet(as.matrix(UMAP_data[Train,-c(1,2)]), factor(UMAP_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                     lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                     type.multinomial="grouped", nfolds=10)
cv_tSNE <- cv.glmnet(as.matrix(tSNE_data[Train,-c(1,2)]), factor(tSNE_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                     lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                     type.multinomial="grouped", nfolds=10)

## Testing



## Evaluations







## Normalization


## Data transformations

# PCA



# tSNE

