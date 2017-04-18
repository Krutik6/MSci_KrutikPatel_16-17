library(ggplot2)
library(gplots)
library(vegan)
#get data
read.csv("~/DISK2/glm/FDR<0.1/Together/mydata.csv", row.names = 1) -> my.data
as.matrix(my.data) -> X
####making cluster
hr <- hclust(as.dist(1-cor(t(X), method = "pearson")) ,method= "average")
hc <- hclust(as.dist(1-cor(X, method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.3))
#mycl <- cutree(hr, k=4)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redblue(75))
##heatmap
#x11(width = 9, height = 7, pointsize = 16)
par(oma =c(5,0,0,0))
heatmap.2(X, main="RBP enrichment profiles", Rowv=as.dendrogram(hr), Colv = NULL, dendrogram="col",
          scale="row", col=myheatcol, density.info="none", trace="none",
          RowSideColors= myClusterSideBar) -> HM
###divide clusters
mycl
mycl[hr$order] -> ex
red <- X[mycl == 1,]
green <- X[mycl == 2,]
blue <- X[mycl == 3,]
Y <- cbind(X, clusterID = mycl)
as.matrix(Y) -> Y

write.csv(x = Y, "~/DISK2/polysome profiling/edgeR_clusters.csv")


###########################sub clusters and analysis

hr <- hclust(as.dist(1-cor(t(red), method = "pearson")) ,method= "average")
hc <- hclust(as.dist(1-cor(red, method="pearson")), method="complete")

mycl <- cutree(hr, h=max(hr$height/1.23))
#mycl <- cutree(hr, k=)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redblue(75))
##heatmap
#x11(width = 9, height = 7, pointsize = 16)
par(oma =c(5,0,0,0))
heatmap.2(red, main="RBP enrichment profiles", Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram="both",
          scale="row", col=myheatcol, density.info="none", trace="none",
          RowSideColors= myClusterSideBar) -> HM

mycl
mycl[hr$order] -> ex
Y <- cbind(red, clusterID = mycl)

write.csv(Y, "~/DISK2/polysome profiling/edgeR.cluster1.csv")

as.data.frame(Y) -> Y

Y[which(Y$clusterID == 1),] -> R_2B
Y[which(Y$clusterID == 2),] -> R_2C
Y[which(Y$clusterID == 3),] -> R_3
Y[which(Y$clusterID == 4),] -> R_4
Y[which(Y$clusterID == 5),]-> R_5
Y[which(Y$clusterID == 6),]-> R_6
rbind(R_3, R_4, R_5, R_6) -> R_2A

write.csv(R_2A, "edgeR.RedA.csv")
write.csv(R_2B, "edgeR.RedB.csv")
write.csv(R_2C, "edgeR.RedC.csv")
#############2A

rownames(group1) %in% rownames(R_2A) -> red_2A
cbind(group1, red_2A) -> red_2A
red_2A[which(red_2A$red_2A == T),] -> red_2A

rownames(group2) %in% rownames(R_2A) -> red_2A2
cbind(group2, red_2A2) -> red_2A2
red_2A2[which(red_2A2$red_2A2 == T),] -> red_2A2

rownames(group3A) %in% rownames(R_2A) -> red_2A3
cbind(group3A, red_2A3) -> red_2A3
red_2A3[which(red_2A3$red_2A3 == T),] -> red_2A3

rownames(group3B) %in% rownames(R_2A) -> red_2A4
cbind(group3B, red_2A4) -> red_2A4
red_2A4[which(red_2A4$red_2A4 == T),] -> red_2A4

rownames(group4A) %in% rownames(R_2A) -> red_2A5
cbind(group4A, red_2A5) -> red_2A5
red_2A5[which(red_2A5$red_2A5 == T),] -> red_2A5

rownames(group4B) %in% rownames(R_2A) -> red_2A6
cbind(group4B, red_2A6) -> red_2A6
red_2A6[which(red_2A6$red_2A6 == T),] -> red_2A6

rownames(group4C) %in% rownames(R_2A) -> red_2A7
cbind(group4C, red_2A7) -> red_2A7
red_2A7[which(red_2A7$red_2A7 == T),] -> red_2A7

#############2B

rownames(group1) %in% rownames(R_2B) -> red_2B
cbind(group1, red_2B) -> red_2B
red_2B[which(red_2B$red_2B == T),] -> red_2B

rownames(group2) %in% rownames(R_2B) -> red_2B2
cbind(group2, red_2B2) -> red_2B2
red_2B2[which(red_2B2$red_2B2 == T),] -> red_2B2

rownames(group3A) %in% rownames(R_2B) -> red_2B3
cbind(group3A, red_2B3) -> red_2B3
red_2B3[which(red_2B3$red_2B3 == T),] -> red_2B3

rownames(group3B) %in% rownames(R_2B) -> red_2B4
cbind(group3B, red_2B4) -> red_2B4
red_2B4[which(red_2B4$red_2B4 == T),] -> red_2B4

rownames(group4A) %in% rownames(R_2B) -> red_2B5
cbind(group4A, red_2B5) -> red_2B5
red_2B5[which(red_2B5$red_2B5 == T),] -> red_2B5

rownames(group4B) %in% rownames(R_2B) -> red_2B6
cbind(group4B, red_2B6) -> red_2B6
red_2B6[which(red_2B6$red_2B6 == T),] -> red_2B6

rownames(group4C) %in% rownames(R_2B) -> red_2B7
cbind(group4C, red_2B7) -> red_2B7
red_2B7[which(red_2B7$red_2B7 == T),] -> red_2B7

#############2C

rownames(group1) %in% rownames(R_2C) -> red_2C
cbind(group1, red_2C) -> red_2C
red_2C[which(red_2C$red_2C == T),] -> red_2C

rownames(group2) %in% rownames(R_2C) -> red_2C2
cbind(group2, red_2C2) -> red_2C2
red_2C2[which(red_2C2$red_2C2 == T),] -> red_2C2

rownames(group3A) %in% rownames(R_2C) -> red_2C3
cbind(group3A, red_2C3) -> red_2C3
red_2C3[which(red_2C3$red_2C3 == T),] -> red_2C3

rownames(group3B) %in% rownames(R_2C) -> red_2C4
cbind(group3B, red_2C4) -> red_2C4
red_2C4[which(red_2C4$red_2C4 == T),] -> red_2C4

rownames(group4A) %in% rownames(R_2C) -> red_2C5
cbind(group4A, red_2C5) -> red_2C5
red_2C5[which(red_2C5$red_2C5 == T),] -> red_2C5

rownames(group4B) %in% rownames(R_2C) -> red_2C6
cbind(group4B, red_2C6) -> red_2C6
red_2C6[which(red_2C6$red_2C6 == T),] -> red_2C6

rownames(group4C) %in% rownames(R_2C) -> red_2C7
cbind(group4C, red_2C7) -> red_2C7
red_2C7[which(red_2C7$red_2C7 == T),] -> red_2C7

#######################################~Blue

hr <- hclust(as.dist(1-cor(t(blue), method = "pearson")) ,method= "average")
hc <- hclust(as.dist(1-cor(blue, method="pearson")), method="complete")

mycl <- cutree(hr, h=max(hr$height/1.23))
#mycl <- cutree(hr, k=)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redblue(75))
##heatmap
#x11(width = 9, height = 7, pointsize = 16)
par(oma =c(5,0,0,0))
heatmap.2(blue, main="RBP enrichment profiles", Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram="both",
          scale="row", col=myheatcol, density.info="none", trace="none",
          RowSideColors= myClusterSideBar) -> HM

mycl
mycl[hr$order] -> ex
blu <- cbind(blue, clusterID = mycl)

write.csv(blu, "~/DISK2/polysome profiling/edgeR.cluster3.csv")

as.data.frame(blu) -> Z
Z[which(Z$clusterID == 1),] -> B2
Z[which(Z$clusterID == 2),] -> Bx
Z[which(Z$clusterID == 3),] -> B3
Z[which(Z$clusterID == 4),] -> B4
Z[which(Z$clusterID == 5),] -> B5
rbind(Bx, B4, B5, B3) -> B1

write.csv(B1, "edgeR.blueA.csv")
write.csv(B2, "edgeR.blueB.csv")

#############3A

rownames(group1) %in% rownames(B1) -> B_3A
cbind(group1, B_3A) -> B_3A
B_3A[which(B_3A$B_3A == T),] -> B_3A

rownames(group2) %in% rownames(B1) -> B_3A2
cbind(group2, B_3A2) -> B_3A2
B_3A2[which(B_3A2$B_3A2 == T),] -> B_3A2

rownames(group3A) %in% rownames(B1) -> B_3A3
cbind(group3A, B_3A3) -> B_3A3
B_3A3[which(B_3A3$B_3A3 == T),] -> B_3A3

rownames(group3B) %in% rownames(B1) -> B_3A4
cbind(group3B, B_3A4) -> B_3A4
B_3A4[which(B_3A4$B_3A4 == T),] -> B_3A4

rownames(group4A) %in% rownames(B1) -> B_3A5
cbind(group4A, B_3A5) -> B_3A5
B_3A5[which(B_3A5$B_3A5 == T),] -> B_3A5

rownames(group4B) %in% rownames(B1) -> B_3A6
cbind(group4B, B_3A6) -> B_3A6
B_3A6[which(B_3A6$B_3A6 == T),] -> B_3A6

rownames(group4C) %in% rownames(B1) -> B_3A7
cbind(group4C, B_3A7) -> B_3A7
B_3A7[which(B_3A7$B_3A7 == T),] -> B_3A7

############3B

rownames(group1) %in% rownames(B2) -> B_3B
cbind(group1, B_3B) -> B_3B
B_3B[which(B_3B$B_3B == T),] -> B_3B

rownames(group2) %in% rownames(B2) -> B_3B2
cbind(group2, B_3B2) -> B_3B2
B_3B2[which(B_3B2$B_3B2 == T),] -> B_3B2

rownames(group3A) %in% rownames(B2) -> B_3B3
cbind(group3A, B_3B3) -> B_3B3
B_3B3[which(B_3B3$B_3B3 == T),] -> B_3B3

rownames(group3B) %in% rownames(B2) -> B_3B4
cbind(group3B, B_3B4) -> B_3B4
B_3B4[which(B_3B4$B_3B4 == T),] -> B_3B4

rownames(group4A) %in% rownames(B2) -> B_3B5
cbind(group4A, B_3B5) -> B_3B5
B_3B5[which(B_3B5$B_3B5 == T),] -> B_3B5

rownames(group4B) %in% rownames(B2) -> B_3B6
cbind(group4B, B_3B6) -> B_3B6
B_3B6[which(B_3B6$B_3B6 == T),] -> B_3B6

rownames(group4C) %in% rownames(B2) -> B_3B7
cbind(group4C, B_3B7) -> B_3B7
B_3B7[which(B_3B7$B_3B7 == T),] -> B_3B7





