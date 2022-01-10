library(pROC)
library(caret)

###############################compare DM+LNM vs. NM############################################################################
# Data import and process groups
dat <- read.csv("dot blot samples for metastatic model.csv", stringsAsFactors = F)
dat$group2 <- dat$Group
dat$group2[dat$group2 != "NM"] <- "M"
dat$Group <- dat$group2

# set trainControl mode to leave-one-out method
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, savePredictions = "final", summaryFunction= twoClassSummary)

#list all methods in caret
#names(getModelInfo())

# Train averaged neural network (avNNet) model with four proteins
model <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine, data = dat, 
               method = "avNNet", trControl = ctrl, metric = "ROC")
print(model)

# Train model with four proteins plus CEA
model_plusCEA <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine + CEA..ng.ml., 
                       data = dat, method = "avNNet", trControl = ctrl, metric = "ROC")
print(model_plusCEA)

# Output prediction probabilities
pred_result <- cbind(dat$Group, model$pred$M[order(as.numeric(model$pred$rowIndex))], 
                     model_plusCEA$pred$M[order(as.numeric(model_plusCEA$pred$rowIndex))], dat$CEA..ng.ml.)
rownames(pred_result) <- dat$ID
colnames(pred_result) <- c("group", "panel", "panel+CEA", "CEA")
write.csv(pred_result, "ROC_M_vs_NM_avNNet.csv")


###############################compare DM vs. NM############################################################################
# Data import
dat <- read.csv("dot blot samples for metastatic model_DMvsNM.csv", stringsAsFactors = F)

# set trainControl mode to leave-one-out method
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, savePredictions = "final", summaryFunction= twoClassSummary)

# Train averaged neural network (avNNet) model with four proteins
model <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine, data = dat, 
               method = "avNNet", trControl = ctrl, metric = "ROC")
print(model)

# Train model with four proteins plus CEA
model_plusCEA <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine + CEA..ng.ml., 
                       data = dat, method = "avNNet", trControl = ctrl, metric = "ROC")
print(model_plusCEA)

# Output prediction probabilities
pred_result <- cbind(dat$Group, model$pred$DM[order(as.numeric(model$pred$rowIndex))], 
                     model_plusCEA$pred$DM[order(as.numeric(model_plusCEA$pred$rowIndex))], dat$CEA..ng.ml.)
rownames(pred_result) <- dat$ID
colnames(pred_result) <- c("group", "panel", "panel+CEA", "CEA")
write.csv(pred_result, "ROC_DM_vs_NM_avNNet.csv")


###############################compare LNM vs. NM############################################################################
# Data import
dat <- read.csv("dot blot samples for metastatic model_LNMvsNM.csv", stringsAsFactors = F)

# Set trainControl mode to leave-one-out method
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, savePredictions = "final", summaryFunction= twoClassSummary)

# Train averaged neural network (avNNet) model with four proteins
model <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine, data = dat, 
               method = "avNNet", trControl = ctrl, metric = "ROC")
print(model)

# Train model with four proteins plus CEA
model_plusCEA <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine + CEA..ng.ml., 
                       data = dat, method = "avNNet", trControl = ctrl, metric = "ROC")
print(model_plusCEA)

# Output prediction probabilities
pred_result <- cbind(dat$Group, model$pred$LNM[order(as.numeric(model$pred$rowIndex))], 
                     model_plusCEA$pred$LNM[order(as.numeric(model_plusCEA$pred$rowIndex))], dat$CEA..ng.ml.)
rownames(pred_result) <- dat$ID
colnames(pred_result) <- c("group", "panel", "panel+CEA", "CEA")
write.csv(pred_result, "ROC_LNM_vs_NM_avNNet.csv")

