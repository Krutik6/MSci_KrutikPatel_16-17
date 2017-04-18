# Load libraries
library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(data.table)
#load file
setwd("~/DISK2/glm/redoing/EAP1/")
#eIF4E.1 <- read.csv("~/DISK2/glm/redoing/eIF4E/eIF4E.1.csv", row.names=1)
eIF4E.2 <- read.csv("~/DISK2/glm/redoing/eIF4E/eIF4E.2.csv", row.names=1)
#eIF4E.3 <- read.csv("~/DISK2/glm/redoing/eIF4E/eIF4E.3.csv", row.names=1)

#make group
group <- rep(c(0,1))

#make EdgeR object
#y <- DGEList(counts=eIF4E.1, group = group)
y <- DGEList(counts=eIF4E.2, group = group)
#y <- DGEList(counts=eIF4E.3, group = group)
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
#eIF4E_allTags.1 <- topTags(lrt,n=length(rownames(lrt$table)))$table
eIF4E_allTags.2 <- topTags(lrt,n=length(rownames(lrt$table)))$table
#eIF4E_allTags.3 <- topTags(lrt,n=length(rownames(lrt$table)))$table

#write
#write.csv(eIF4E_allTags.1, "~/DISK2/glm/redoing/eIF4E//Tags.eIF4E.1")
write.csv(eIF4E_allTags.2, "~/DISK2/glm/redoing/eIF4E/Tags.eIF4E.2")
#write.csv(eIF4E_allTags.3, "~/DISK2/glm/redoing/eIF4E/Tags.eIF4E.3")