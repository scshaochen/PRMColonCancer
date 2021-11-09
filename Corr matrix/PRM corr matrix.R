library(ggplot2)
library(reshape2)


d1<-read.table("PRM_23up for corr.csv", header=TRUE, sep=",", as.is=T)
# Compute the correlation matrix
cormat <- round(cor(d1,method="spearman"),2)
head(cormat)

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri
View(upper_tri)

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Create a ggheatmap and visualization
ggheatmap <- ggplot(data=melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "light gray", high = "red", mid = "white", 
                       midpoint = 0.45, limit = c(-0.15,1.0), space = "Lab", 
                       name="Spearman\nCorrelation") + xlab("") + ylab("")+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   size = 10, hjust = 1))+
  coord_fixed()

print(ggheatmap)

ggheatmap + 
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")

