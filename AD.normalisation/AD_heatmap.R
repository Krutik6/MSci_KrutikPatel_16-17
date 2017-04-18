library(ggplot2)
library(gplots)
#AD
read.csv("~/DISK2/AD/heatmap/My.data.csv", row.names = 1) -> my.data

as.matrix(my.data) -> Y

hr <- hclust(as.dist(1-cor(t(Y), method = "pearson")) ,method= "average")
hc <- hclust(as.dist(1-cor(Y, method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.3))
#mycl <- cutree(hr, k=4)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redblue(75))


heatmap.2(Y, main="AD normalised RBP enrichment profiles", Rowv=as.dendrogram(hr), NULL, dendrogram="none",
          scale="row", col=myheatcol, density.info="none", trace="none",
          RowSideColors= myClusterSideBar) -> HM
par(oma = c(5, 0, 0, 0))

mycl
mycl[hr$order] -> ex
blue <- Y[mycl == 2,]
red <- Y[mycl == 1,]
fooY <- cbind(Y, clusterID = mycl)
write.csv(fooY, "~/DISK2/polysome profiling/AD.clusters.csv")

read.csv("~/DISK2/glm/FDR<0.1/compare with X/group1.csv", row.names = 1) -> group1

#############Finding Group1 genes

rownames(group1) %in% rownames(red) -> G1_red
cbind(group1, G1_red) -> G1_red
G1_red[which(G1_red$G1_red == T),] -> G1_red

rownames(group1) %in% rownames(blue) -> G1_blue
cbind(group1, G1_blue) -> G1_blue
G1_blue[which(G1_blue$G1_blue == T),] -> G1_blue

######### G1 is mostly in blue 74/193, whilst 20/193 in red. 99 not found.

read.csv("~/DISK2/glm/FDR<0.1/compare with X/Group2.csv", row.names = 1) -> group2

rownames(group2) %in% rownames(red) -> G2_red
cbind(group2, G2_red) -> G2_red
G2_red[which(G2_red$G2_red == T),] -> G2_red

rownames(group2) %in% rownames(blue) -> G2_blue
cbind(group2, G2_blue) -> G2_blue
G2_blue[which(G2_blue$G2_blue == T),] -> G2_blue

############most of group2 was in blue 45/195, and 24/195 in red. 126 were not present. 

read.csv("~/DISK2/glm/FDR<0.1/compare with X/group_3A.csv", row.names = 1) -> group3A

rownames(group3A) %in% rownames(red) -> G3A_red
cbind(group3A, G3A_red) -> G3A_red
G3A_red[which(G3A_red$G3A_red == T),] -> G3A_red

rownames(group3A) %in% rownames(blue) -> G3A_blue
cbind(group3A, G3A_blue) -> G3A_blue
G3A_blue[which(G3A_blue$G3A_blue == T),] -> G3A_blue

######################Most of 3A was in blue 176/245 and in red 18/245. 51 not present.

read.csv("~/DISK2/glm/FDR<0.1/compare with X/group3B.csv", row.names = 1) -> group3B

rownames(group3B) %in% rownames(red) -> G3B_red
cbind(group3B, G3B_red) -> G3B_red
G3B_red[which(G3B_red$G3B_red == T),] -> G3B_red

rownames(group3B) %in% rownames(blue) -> G3B_blue
cbind(group3B, G3B_blue) -> G3B_blue
G3B_blue[which(G3B_blue$G3B_blue == T),] -> G3B_blue

######################Most of 3B was in blue 127/148 and in red 1/148. 20 were not found.

read.csv("~/DISK2/glm/FDR<0.1/compare with X/group4A.csv", row.names = 1) -> group4A

rownames(group4A) %in% rownames(red) -> G4A_red
cbind(group4A, G4A_red) -> G4A_red
G4A_red[which(G4A_red$G4A_red == T),] -> G4A_red

rownames(group4A) %in% rownames(blue) -> G4A_blue
cbind(group4A, G4A_blue) -> G4A_blue
G4A_blue[which(G4A_blue$G4A_blue == T),] -> G4A_blue

###################### Most of 4A was in blue 539/629, 25/629 in red. 65  not found.

read.csv("~/DISK2/glm/FDR<0.1/compare with X/group4B.csv", row.names = 1) -> group4B

rownames(group4B) %in% rownames(red) -> G4B_red
cbind(group4B, G4B_red) -> G4B_red
G4B_red[which(G4B_red$G4B_red == T),] -> G4B_red

rownames(group4B) %in% rownames(blue) -> G4B_blue
cbind(group4B, G4B_blue) -> G4B_blue
G4B_blue[which(G4B_blue$G4B_blue == T),] -> G4B_blue

######################## 4B is in red 106/1014, and 480/1014 in blue. 428 not found.

read.csv("~/DISK2/glm/FDR<0.1/compare with X/group4C.csv", row.names = 1) -> group4C

rownames(group4C) %in% rownames(red) -> G4C_red
cbind(group4C, G4C_red) -> G4C_red
G4C_red[which(G4C_red$G4C_red == T),] -> G4C_red

rownames(group4C) %in% rownames(blue) -> G4C_blue
cbind(group4C, G4C_blue) -> G4C_blue
G4C_blue[which(G4C_blue$G4C_blue == T),] -> G4C_blue

######################mostly in blue 103/337 and some in red 192/337. 42 not found.



##############################find sub groups for blue
hr <- hclust(as.dist(1-cor(t(blue), method = "pearson")) ,method= "average")
hc <- hclust(as.dist(1-cor(blue, method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/3.4))
#mycl <- cutree(hr, k=4)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- rev(redblue(75))


heatmap.2(blue, main="AD normalised RBP enrichment profiles", Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram="both",
          scale="row", col=myheatcol, density.info="none", trace="none",
          RowSideColors= myClusterSideBar) -> HM
par(oma = c(5, 0, 0, 0))

mycl
mycl[hr$order] -> ex
blue[mycl == 1,] -> bl1B
blue[mycl == 2,] -> cl2
blue[mycl == 3,] -> cl3
blue[mycl == 4,] -> cl4
blue[mycl == 5,] -> cl5

rbind(cl2, cl5) -> bl1A
rbind(cl4, cl3) -> bl1C
###########################################Blue_cluster 1
rownames(group1) %in% rownames(bl1A) -> b_1a
cbind(group1, b_1a) -> b_1a
b_1a[which(b_1a$b_1a == T),] -> b_1a

rownames(group2) %in% rownames(bl1A) -> b_1a2
cbind(group2, b_1a2) -> b_1a2
b_1a2[which(b_1a2$b_1a2 == T),] -> b_1a2

rownames(group3A) %in% rownames(bl1A) -> b_1a3a
cbind(group3A, b_1a3a) -> b_1a3a
b_1a3a[which(b_1a3a$b_1a3a == T),] -> b_1a3a

rownames(group3B) %in% rownames(bl1A) -> b_1a3b
cbind(group3B, b_1a3b) -> b_1a3b
b_1a3b[which(b_1a3b$b_1a3b == T),] -> b_1a3b

rownames(group4A) %in% rownames(bl1A) -> b_1a4a
cbind(group4A, b_1a4a) -> b_1a4a
b_1a4a[which(b_1a4a$b_1a4a == T),] -> b_1a4a

rownames(group4B) %in% rownames(bl1A) -> b_1a4b
cbind(group4B, b_1a4b) -> b_1a4b
b_1a4b[which(b_1a4b$b_1a4b == T),] -> b_1a4b

rownames(group4C) %in% rownames(bl1A) -> b_1a4C
cbind(group4C, b_1a4C) -> b_1a4C
b_1a4C[which(b_1a4C$b_1a4C == T),] -> b_1a4C

###########################################Blue_cluster 2
rownames(group1) %in% rownames(bl1B) -> b_2a
cbind(group1, b_2a) -> b_2a
b_2a[which(b_2a$b_2a == T),] -> b_2a

rownames(group2) %in% rownames(bl1B) -> b_2a2
cbind(group2, b_2a2) -> b_2a2
b_2a2[which(b_2a2$b_2a2 == T),] -> b_2a2

rownames(group3A) %in% rownames(bl1B) -> b_2a3a
cbind(group3A, b_2a3a) -> b_2a3a
b_2a3a[which(b_2a3a$b_2a3a == T),] -> b_2a3a

rownames(group3B) %in% rownames(bl1B) -> b_2a3b
cbind(group3B, b_2a3b) -> b_2a3b
b_2a3b[which(b_2a3b$b_2a3b == T),] -> b_2a3b

rownames(group4A) %in% rownames(bl1B) -> b_2a4a
cbind(group4A, b_2a4a) -> b_2a4a
b_2a4a[which(b_2a4a$b_2a4a == T),] -> b_2a4a

rownames(group4B) %in% rownames(bl1B) -> b_2a4b
cbind(group4B, b_2a4b) -> b_2a4b
b_2a4b[which(b_2a4b$b_2a4b == T),] -> b_2a4b

rownames(group4C) %in% rownames(bl1B) -> b_2a4C
cbind(group4C, b_2a4C) -> b_2a4C
b_2a4C[which(b_2a4C$b_2a4C == T),] -> b_2a4C

###########################################Blue_cluster 3

rownames(group1) %in% rownames(bl1C) -> b_3a
cbind(group1, b_3a) -> b_3a
b_3a[which(b_3a$b_3a == T),] -> b_3a

rownames(group2) %in% rownames(bl1C) -> b_3a2
cbind(group2, b_3a2) -> b_3a2
b_3a2[which(b_3a2$b_3a2 == T),] -> b_3a2

rownames(group3A) %in% rownames(bl1C) -> b_3a3a
cbind(group3A, b_3a3a) -> b_3a3a
b_3a3a[which(b_3a3a$b_3a3a == T),] -> b_3a3a

rownames(group3B) %in% rownames(bl1C) -> b_3a3b
cbind(group3B, b_3a3b) -> b_3a3b
b_3a3b[which(b_3a3b$b_3a3b == T),] -> b_3a3b

rownames(group4A) %in% rownames(bl1C) -> b_3a4a
cbind(group4A, b_3a4a) -> b_3a4a
b_3a4a[which(b_3a4a$b_3a4a == T),] -> b_3a4a

rownames(group4B) %in% rownames(bl1C) -> b_3a4b
cbind(group4B, b_3a4b) -> b_3a4b
b_3a4b[which(b_3a4b$b_3a4b == T),] -> b_3a4b

rownames(group4C) %in% rownames(bl1C) -> b_3a4C
cbind(group4C, b_3a4C) -> b_3a4C
b_3a4C[which(b_3a4C$b_3a4C == T),] -> b_3a4C


