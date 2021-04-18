#Load necessary packages
suppressMessages(library(TeachingDemos))


#Load Leukemia dataset
leukemia <- read.csv("http://web.stanford.edu/~hastie/CASI_files/DATA/leukemia_small.csv")
#ALL = 0, AML = 1
class <- c(rep(0,27),rep(1,11),rep(0,11),rep(1,5),rep(0,2),rep(1,2),0,rep(1,7),rep(0,6))
#Give patients numeric ids
id <- c(1:72)
#Add above changes to DF leukemia 
colnames(leukemia)<-NULL 
leuk <- rbind(class,leukemia) 
colnames(leuk)<-c(1:72)