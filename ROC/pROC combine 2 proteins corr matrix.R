library(pROC)

#### import data
PRM <- read.csv("PRM_23up del9_diag.csv",header = TRUE,row.names = 1)

####calculate AUC values of single protein to generate Fig. 3c
res <- matrix(nrow = ncol(PRM), ncol = 1)
rownames(res) <- colnames(PRM) 
colnames(res) <- c("AUC")
for (i in 1:ncol(PRM)){
  z <- roc(PRM$diag_group, PRM[,i])
  res[i,1] <- z$auc
}

write.table(res, file = "PRM_23up_del9_diag ROC AUC.txt", sep="\t")

##calculate AUC values of the combination of any two proteins
PRM1 <- read.csv("PRM_23up del9_diag_overlap8.csv",header = TRUE,row.names = 1) ##using the overlapping proteins with the top 10 highest mean decrease accuracy values and AUC values
attach(PRM1)
res2 <- matrix(nrow = ncol(PRM1), ncol = ncol(PRM1))
rownames(res2) <- colnames(PRM1)
colnames(res2) <- colnames(PRM1)
for (i in 1:ncol(PRM1)){
  for(j in 1:ncol(PRM1)){
    predictor <- predict(glm(diag_group ~ PRM1[,i] + PRM1[,j],family="binomial"))
    z <- roc(PRM1$diag_group, predictor)
    res2[i,j] <- z$auc
  }
}
res3 <- res2[-1,-1]
write.table(cbind(rownames(res3), res3), "new PRM_23up del9_diag model_ROC AUC_comb 2.csv", sep = ",", row.names = F, col.names = T)
detach(PRM1)

## prepare for drawing corr. matrix
library(ggplot2)
library(reshape2)

cormat<-read.table("new PRM_23up del9_diag model_ROC AUC_comb 2.csv", header=TRUE, sep=",", row.names=1, as.is=T)

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
upper_tri
upper_tri2 <- as.matrix(upper_tri)

# Use correlation between variables as distance
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
upper_tri2 <- as.matrix(upper_tri)
melted_cormat <- melt(upper_tri2, na.rm = TRUE)

# Create a ggheatmap and add correlation coefficients on the heatmap
ggheatmap <- ggplot(data=melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "deepskyblue1", high = "firebrick1", mid = "white", 
                       midpoint = 0.75, limit = c(0.6,0.92), space = "Lab",
                       name="AUC of combined variables") + xlab("") + ylab("")+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   size = 12, hjust = 1), 
        axis.text.y = element_text(size=12))+
  coord_fixed()

print(ggheatmap)

ggheatmap + 
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal")


################################################################################################################################repeat the analysis for metastatic model using the data file "PRM_23up del9_meta.csv.csv" and "PRM_23up del9_meta_overlap8.csv"