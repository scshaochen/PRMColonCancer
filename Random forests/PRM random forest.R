install.packages("randomForest")
library(randomForest)
set.seed(666)

## data import
PRM <- read.csv("PRM_23up del9_diag_RF.csv",header = TRUE,sep = ",",row.names = 1)
PRM$diag_group <- as.factor(PRM$diag_group)

## create the forest
fitTest<-randomForest(diag_group~.,data = PRM,importance=TRUE, proximity=TRUE)
## view the forest results
print(fitTest)


## Look at variable importance based on mean decrease in accuracy index:
ImportanceOrder<-importance(fitTest,type = 1)
ImportanceOrder<-round(ImportanceOrder[order(ImportanceOrder, decreasing = TRUE),],2)
varImpPlot(fitTest)  

###100 random forest replicates of Mean Decrease Accuracy to minimize randomness
result <- as.data.frame(matrix(nc=100,nr=14)) 
  for(i in 1:100){
  fitTest<-randomForest(diag_group~.,data = PRM,importance=TRUE,mtry=5,ntree=1500,na.action = na.omit, proximity=TRUE)
  ImportanceOrder<-importance(fitTest,type = 1)
  ImportanceOrder<-round(ImportanceOrder[order(ImportanceOrder, decreasing = TRUE),],2)
  varImpPlot(fitTest) 
  print(fitTest)
  print(ImportanceOrder)
  result[,i]<-ImportanceOrder[order(names(ImportanceOrder))]
}
rownames(result) <- names(ImportanceOrder)[order(names(ImportanceOrder))]
write.table(result,file="D:Downloads//diag_accuracy 100x.csv",append = FALSE, quote = TRUE, sep = ",",eol = "\n", na ="NA",dec = ".",row.names = TRUE, col.names =TRUE, qmethod = c("escape","double")) 
print(fitTest)


######################################################
## repeat the analysis for metastatic model using the data file "PRM_23up del9_meta_RF.csv"