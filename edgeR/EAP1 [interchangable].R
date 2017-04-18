# Load libraries
library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(data.table)
#load file
setwd("~/DISK2/glm/redoing/EAP1/")
#EAP1.1 <- read.csv("~/DISK2/glm/redoing/EAP1/EAP1.1.csv", row.names=1)
EAP1.2 <- read.csv("~/DISK2/glm/redoing/EAP1/EAP1.2.csv", row.names=1)
#EAP1.3 <- read.csv("~/DISK2/glm/redoing/EAP1/EAP1.3.csv", row.names=1)

#make group
group <- rep(c(0,1))

#make EdgeR object
#y <- DGEList(counts=EAP1.1, group = group)
y <- DGEList(counts=EAP1.2, group = group)
#y <- DGEList(counts=EAP1.3, group = group)
keep <- rowSums(cpm(y) > 20) >= 1
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y,method='TMM')

#design
design <- model.matrix(~group)
rownames(design) <- colnames(y)

#estimate DE
y <- estimateGLMCommonDisp(y)
y <- estimateGLMTrendedDisp(y)
y <- estimateGLMTagwiseDisp(y)

#GLM model
fit <- glmFit(y, design)
lrt <- glmLRT(fit)

#results
topTags(lrt)
#EAP1_allTags.1 <- topTags(lrt,n=length(rownames(lrt$table)))$table
EAP1_allTags.2 <- topTags(lrt,n=length(rownames(lrt$table)))$table
#EAP1_allTags.3 <- topTags(lrt,n=length(rownames(lrt$table)))$table

#write
#write.csv(EAP1_allTags.1, "~/DISK2/glm/redoing/EAP1//Tags.EAP1.1")
write.csv(EAP1_allTags.2, "~/DISK2/glm/redoing/EAP1//Tags.EAP1.2")
#write.csv(EAP1_allTags.3, "~/DISK2/glm/redoing/EAP1/Tags.EAP1.3")