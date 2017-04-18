library(samr)
library(data.table)

#read.csv("~/DISK2/AD/Data/RBPs/AD_norm_data/CAF20.csv", row.names = 1) -> X
#read.csv("~/DISK2/AD/Data/RBPs/AD_norm_data/EAP1.csv", row.names = 1) -> X
#read.csv("~/DISK2/AD/Data/RBPs/AD_norm_data/eIF4E.csv", row.names = 1) -> X
#read.csv("~/DISK2/AD/Data/RBPs/AD_norm_data/eIF4G1.csv", row.names = 1) -> X
#read.csv("~/DISK2/AD/Data/RBPs/AD_norm_data/eIF4G2.csv", row.names = 1) -> X
read.csv("~/DISK2/AD/Data/RBPs/AD_norm_data/Pab1.csv", row.names = 1) -> X

read.csv("~/DISK2/AD/Data/names.csv") -> genes

xcol <- ncol(X)
xrow <- nrow(X)-1
X2 <- X[1:xrow, 1:xcol]
y1 <-c(rep(1,3), rep(2,3))
data=list(x=as.matrix(X2), y=y1,logged2 =TRUE)

samr.obj <- samr(data, resp.type = "Two class unpaired", nperms = 100)
delta.table <- samr.compute.delta.table(samr.obj, min.foldchange = 0.1, nvals = 200)
siggenes.table <- samr.compute.siggenes.table(samr.obj, del = 0, data, delta.table, all.genes = TRUE)

A <- siggenes.table$genes.up 
B <- siggenes.table$genes.lo
c <- rbind(A, B)
lo <- c[as.numeric(c[,8])<1,]

as.data.table(lo) -> fdr

fdr$Row %in% rownames(genes) -> D
cbind(fdr, D) -> E
E[which(E$D == TRUE),] -> F
F$D <- NULL

rownames(genes) %in% fdr$Row -> fdrow
cbind(genes, fdrow) -> comb
comb[which(comb$fdrow == TRUE),] -> yes

yes$fdrow <- NULL
cbind(fdr, yes) -> genes

genes$Row <- genes$`Gene ID` <- genes$`Gene Name` <- genes$`Score(d)` <- genes$`Numerator(r)` <- genes$`Denominator(s+s0)` <- NULL

#write.csv(genes, "~/DISK2/AD/Data/RBPs/1%/CAF20_DE.csv")
#write.csv(genes, "~/DISK2/AD/Data/RBPs/1%/EAP1_DE.csv")
#write.csv(genes, "~/DISK2/AD/Data/RBPs/1%/eIF4E_DE.csv")
#write.csv(genes, "~/DISK2/AD/Data/RBPs/1%/eIF4G1_DE.csv")
#write.csv(genes, "~/DISK2/AD/Data/RBPs/1%/eIF4G2_DE.csv")
write.csv(genes, "~/DISK2/AD/Data/RBPs/1%/Pab1_DE.csv")
