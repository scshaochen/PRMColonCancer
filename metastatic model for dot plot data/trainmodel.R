library(pROC)
library(caret)

# read data and process groups
dat <- read.csv("D:/work/sunyulin/dot blot samples for metastatic model.csv", stringsAsFactors = F)
#dat <- dat[dat$Group != "DM", ]
dat$group2 <- dat$Group
dat$group2[dat$group2 != "NM"] <- "M"
dat$Group <- dat$group2

# set training mode to leave-one-out method
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, savePredictions = T)

#list all methods in caret
#names(getModelInfo())

#train logistic regression model with four proteins
model <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine, data = dat, method = "glm", trControl = ctrl)
print(model)
model.predict <- predict(model, dat, type="prob")

#train model with four proteins plus CEA
model_plusCEA <- train(factor(Group) ~ CORO1C.creatinine + RAD23B.creatinine + GSPT2.creatinine + NDN.creatinine + CEA..ng.ml., data = dat, method = "glm", trControl = ctrl)
print(model_plusCEA)
#confusionMatrix(model$pred$pred[model$pred$mtry == 2], model$pred$obs[model$pred$mtry == 2])

# output prediction probabilities
pred_result <- cbind(dat$Group, model$pred$DM[model$pred$mtry == 2], model_plusCEA$pred$DM[model_plusCEA$pred$mtry == 5], dat$CEA..ng.ml.)
pred_result <- cbind(dat$Group, model$pred$M, model_plusCEA$pred$M, dat$CEA..ng.ml.)
rownames(pred_result) <- dat$ID
colnames(pred_result) <- c("group", "panel", "panel+CEA", "CEA")

# ROC plot
roc_cea <- roc(dat$Group, as.numeric(pred_result[, "CEA"]), levels = c("NM", "M"))
roc_panel <- roc(dat$Group, as.numeric(pred_result[, "panel"]), levels = c("NM", "M"))
roc_panel_cea <- roc(dat$Group, as.numeric(pred_result[, "panel+CEA"]), levels = c("NM", "M"))

coords(roc_cea, "best", ret = c("specificity", "sensitivity", "accuracy"))

plot(roc_cea)
plot(roc_panel, add = T, col = "red")
plot(roc_panel_cea, add = T, col = "blue")
legend("right", legend = paste(c("CEA", "panel", "panel+CEA"), ", AUC =", round(c(auc(roc_cea), auc(roc_panel), auc(roc_panel_cea)), 2)), fill = c("black", "red", "blue"))

write.csv(pred_result, "D:/work/sunyulin/ROC_M_vs_NM.csv")

# calculate ROC for individual protein
roc_x <- roc(dat$Group, dat$CORO1C.creatinine, levels = c("NM", "M"))
coords(roc_x, "best", ret = c("specificity", "sensitivity", "accuracy"))
auc(roc_x)
cbind(roc_x$sensitivities, roc_x$specificities)
