#### Load necessary packages#### 
# install.packages('glmnet')
# install.packages("umap") # click no to compiling from source
# install.packages("Rtsne")
# install.packages("caret")
# install.packages("pROC")
# install.packages('ggpubr')
# devtools::install_github("collectivemedia/tictoc")
library(glmnet)
library(dplyr)
library(umap)
library(Rtsne)
library(ggplot2)
library(caret)
library(pROC)
library(ggpubr)
library(tictoc)

## Load data
#data <- read.csv('/Users/lcmas/Desktop/statsproject/data.csv', header = TRUE, stringsAsFactors = FALSE)
data <- read.csv('/homes/anovak9/stats/data.csv', header = TRUE, stringsAsFactors = FALSE)
#data <- read.csv('C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Loyola/Classes/STAT 437/Final-Project/TCGA-PANCAN-HiSeq-801x20531/data.csv', header = TRUE, stringsAsFactors = FALSE)
colnames(data)[1] <- 'sample'
#label <- read.csv('/Users/lcmas/Desktop/statsproject/labels.csv', header = TRUE, stringsAsFactors = FALSE)
label <- read.csv('/homes/anovak9/stats/labels.csv', header = TRUE, stringsAsFactors = FALSE)
#label <- read.csv('C:/Users/Tristan Kosciuch/OneDrive - Loyola University Chicago/Loyola/Classes/STAT 437/Final-Project/TCGA-PANCAN-HiSeq-801x20531/labels.csv') #stringsAsFactors = FALSE
colnames(label)[1] <- 'sample'
# str(data)
# summary(data)


#### Pre-processing/Normalization #### 
tic('Data4 ordination')
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
#Keeping 15398 genes 
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
toc()

#### Ordination #### 

## PCA
tic('PCA ordination')
PCA <- prcomp(data4[,-c(1:2)],center = TRUE, scale. = TRUE)
# Select the variables that capture  ~0.95 of the variance
V.explained <- cumsum(PCA$sdev^2/sum(PCA$sdev^2))
PCA_data <- data.frame(PCA$x[,which(V.explained <= 0.95)])
# Add sample data back to the PCA data
PCA_data <- cbind(label,PCA_data)
toc()

#convex hull for plot
hull_PCA <- PCA_data %>% group_by(Class) %>% slice(chull(PC1,PC2))
# Create a plot using the first two PCs
# PCA_plot <- autoplot(PCA, data = data4, colour = "Class", frame = TRUE)
PCA_plot <- ggplot(data = PCA_data, aes(x = PC1,y = PC2)) + aes(color = Class, fill = Class)+ geom_point() + geom_polygon(data = hull_PCA, alpha = 0.4) 



## UMAP

tic('UMAP ordination, 2 dims')
# Umap with two dimensions
UMAP <- umap(data4[,-c(1:2)])
UMAP_data <- data.frame(cbind(label,UMAP$layout))
toc()

tic('UMAP ordination, 10 dims')
# Umap with 10 dimensions
UMAP_10 <- umap(data4[,-c(1:2)], n_components = 10)
UMAP_10_data <- data.frame(cbind(label,UMAP_10$layout))
toc()

tic('UMAP ordination, PCA dims')
# Umap with the same number of dimensions as PCA_data
UMAP_nPCA <- umap(data4[,-c(1:2)], n_components = length(which(V.explained <= 0.95)))
UMAP_nPCA_data <- data.frame(cbind(label,UMAP_nPCA$layout))
toc()

#convex hull for plot
hull_UMAP <- UMAP_data %>% group_by(Class) %>% slice(chull(X1,X2))
# Create a plot using the first two Axis
UMAP_plot <- ggplot(data = UMAP_data, aes(x = X1,y = X2)) + aes(color = Class, fill = Class)+ geom_point() + geom_polygon(data = hull_UMAP, alpha = 0.4) + labs(x = 'UMAP1', y = 'UMAP2')


## tSNE

tic('tSNE ordination')
# theta = 0.0 is an exact tSNE
tSNE <- Rtsne(as.matrix(data4[,-c(1:2)]), theta=0.0)
tSNE_data <- cbind(label,tSNE$Y)
toc()

#convex hull for plot
hull_tSNE <- tSNE_data %>% group_by(Class) %>% slice(chull(`1`,`2`))
# Create a plot using the first two Axis
tSNE_plot <- ggplot(data = tSNE_data, aes(x = `1`,y = `2`)) + aes(color = Class, fill = Class) + geom_point()+ geom_polygon(data = hull_tSNE, alpha = 0.4) + labs(x = 'tSNE1', y = 'tSNE2')


## Figure of ordination plots
figure <- ggarrange(PCA_plot, UMAP_plot, tSNE_plot, labels = c('A', 'B', 'C'), common.legend = TRUE, legend = 'right')
annotate_figure(figure, top = text_grob('Clustering of data after implementing different ordination methods \n', face = 'bold', size = 18))



#### Train Test Splitting ####

set.seed(134)
# Train test split, current ratio is 9:1 train:test
Train <- createDataPartition(data4$Class, p=0.9, list=FALSE)


#### Model fitting and cross validation within the test set ####

# 10-fold cross validation
tic('CV data4')
cv_data4 <- cv.glmnet(as.matrix(data4[Train,-c(1,2)]), factor(data4[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                      lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                      type.multinomial="grouped", nfolds=10)
toc()
tic('CV PCA')
cv_PCA <- cv.glmnet(as.matrix(PCA_data[Train,-c(1,2)]), factor(PCA_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                    lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                    type.multinomial="grouped", nfolds=10)
toc()
tic('CV UMAP2')
cv_UMAP <- cv.glmnet(as.matrix(UMAP_data[Train,-c(1,2)]), factor(UMAP_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                     lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                     type.multinomial="grouped", nfolds=10)
toc()
tic('CV UMAP10')
cv_UMAP_10 <- cv.glmnet(as.matrix(UMAP_10_data[Train,-c(1,2)]), factor(UMAP_10_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                     lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                     type.multinomial="grouped", nfolds=10)
toc()
tic('CV UMAP_nPCA')
cv_UMAP_nPCA <- cv.glmnet(as.matrix(UMAP_nPCA_data[Train,-c(1,2)]), factor(UMAP_nPCA_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                     lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                     type.multinomial="grouped", nfolds=10)
toc()
tic('CV tSNE')
cv_tSNE <- cv.glmnet(as.matrix(tSNE_data[Train,-c(1,2)]), factor(tSNE_data[Train,c(2)]), family="multinomial", alpha=.93, thresh = 1e-07,
                     lambda = c(0.1,0.05,0.01,0.005,0.001, 0.0005, 0.0001),
                     type.multinomial="grouped", nfolds=10)
toc()

#### Predicting on validation data ####

validate <- function(model,data,train_set){
  
  # use the best lambda from cross validation, lambda.min, during prediction 
  fit <- predict(model, as.matrix(data[-train_set,-c(1,2)]), s = model$lambda.min, type = "class")
  
  fit <- cbind(data[-train_set,c(1,2)],fit)
  
  colnames(fit)[3] <- "Predicted"
  
  return(fit)
  
}

tic('Val Data4')
val_data4 <- validate(cv_data4,data4,Train)
toc()
tic('Val PCA')
val_PCA <- validate(cv_PCA,PCA_data,Train)
toc()
tic('Val UMAP2')
val_UMAP <- validate(cv_UMAP,UMAP_data,Train)
toc()
tic('Val UMAP10')
val_UMAP_10 <- validate(cv_UMAP_10,UMAP_10_data,Train)
toc()
tic('Val UMAP_nPCA')
val_UMAP_nPCA <- validate(cv_UMAP_nPCA,UMAP_nPCA_data,Train)
toc()
tic('Val tSNE')
val_tSNE <- validate(cv_tSNE,tSNE_data,Train)
toc()

#### Evaluations ####

## Function to check how many the model got right, if its 0 all were classified correctly - I know this code could be cleaner...
length(which(val_data4$Class==val_data4$Predicted)) - length(val_data4$Class)
length(which(val_PCA$Class==val_PCA$Predicted)) - length(val_PCA$Class)
length(which(val_UMAP$Class==val_UMAP$Predicted)) - length(val_UMAP$Class)
length(which(val_UMAP_10$Class==val_UMAP_10$Predicted)) - length(val_UMAP_10$Class)
length(which(val_UMAP_nPCA$Class==val_UMAP_nPCA$Predicted)) - length(val_UMAP_nPCA$Class)
length(which(val_tSNE$Class==val_tSNE$Predicted)) - length(val_tSNE$Class)

## Time complexity of different classifiers 
times <- read.csv('/homes/anovak9/stats/Stats/times.csv', header = FALSE, stringsAsFactors = FALSE)
rownames(times) <- c('model', 'ordination', 'cv', 'val', 'tot_runtime')
times <- times[,2:7]
times
#Only use total runtimes for plot
tot_time <- times[c(1,5),]
tot_time <- as.data.frame(t(tot_time))
tot_time
#Plot
ggplot(tot_time, aes(x = model, y = tot_runtime))+
  geom_bar(stat = 'identity', fill = 'deepskyblue2', width = 0.5, col = 'black')+
  labs(title = 'Total run time of each classifier', x = '\n Classifier Type', y = 'Time (sec) \n')+
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(size = 18, face = 'bold', hjust = 0.5), axis.text.y = element_text(size = 13))


#### ROC Curves ####

# Data4 ROC

data4_num <- transform(val_data4, Predicted = factor(Predicted, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)) ,  Class = factor(Class, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)))
data4_roc <- multiclass.roc(as.numeric(data4_num$Class), as.numeric(data4_num$Predicted))
data4_rocs <- data4_roc[['rocs']]
plot.roc(data4_rocs[[1]], col = 2, lwd = 3, print.auc = TRUE, add = FALSE, xlab = "Specificity", ylab = "Sensitivity")


# PCA ROC
PCA_num <- transform(val_PCA, Predicted = factor(Predicted, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)) ,  Class = factor(Class, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)))
PCA_roc <- multiclass.roc(as.numeric(PCA_num$Class), as.numeric(PCA_num$Predicted))
PCA_rocs <- PCA_roc[['rocs']]
plot.roc(PCA_rocs[[1]], col = 3, lwd = 3, print.auc = TRUE, add = FALSE, xlab = "Specificity", ylab = "Sensitivity")


#UMAP ROC

UMAP_num <- transform(val_UMAP, Predicted = factor(Predicted, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)) ,  Class = factor(Class, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)))
UMAP_roc <- multiclass.roc(as.numeric(UMAP_num$Class), as.numeric(UMAP_num$Predicted))
UMAP_rocs <- UMAP_roc[['rocs']]
plot.roc(UMAP_rocs[[1]], col = 4, lwd = 3, print.auc = TRUE, add = FALSE, xlab = "Specificity", ylab = "Sensitivity")


#UMAP10 ROC

UMAP_10_num <- transform(val_UMAP_10, Predicted = factor(Predicted, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)) ,  Class = factor(Class, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)))
UMAP_10_roc <- multiclass.roc(as.numeric(UMAP_10_num$Class), as.numeric(UMAP_10_num$Predicted))
UMAP_10_rocs <- UMAP_10_roc[['rocs']]
plot.roc(UMAP_10_rocs[[1]], col = 5, lwd = 3, print.auc = TRUE, add = FALSE, xlab = "Specificity", ylab = "Sensitivity")


#UMAP_nPCA ROC
UMAP_nPCA_num <- transform(val_UMAP_nPCA, Predicted = factor(Predicted, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)) ,  Class = factor(Class, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)))
UMAP_nPCA_roc <- multiclass.roc(as.numeric(UMAP_nPCA_num$Class), as.numeric(UMAP_nPCA_num$Predicted))
UMAP_nPCA_rocs <- UMAP_nPCA_roc[['rocs']]
plot.roc(UMAP_nPCA_rocs[[1]], col = 6, lwd = 3, print.auc = TRUE, add = FALSE, xlab = "Specificity", ylab = "Sensitivity")

#tSNE ROC

tSNE_num <- transform(val_tSNE, Predicted = factor(Predicted, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)) ,  Class = factor(Class, levels = c("BRCA", "LUAD", "PRAD", "COAD", "KIRC"),labels = c(1, 2, 3, 4, 5)))
tSNE_roc <- multiclass.roc(as.numeric(tSNE_num$Class), as.numeric(tSNE_num$Predicted))
tSNE_rocs <- tSNE_roc[['rocs']]
plot.roc(tSNE_rocs[[1]], col = 7, lwd = 3, print.auc = TRUE, add = FALSE, xlab = "Specificity", ylab = "Sensitivity")



