#!/usr/local/bin/Rscript

library(beadarray)

load("results/BSData.quantile.RData")

limma <- read.csv("results/limma_results.csv")

E <- exprs(BSData.quantile)

inds1 <- which(abs(limma[,"logFC"])>1)
inds2 <-which(limma[,"adj.P.Val"]<0.001)

inds <-intersect(inds1, inds2)

ids <- limma[inds,1]

ids <- as.character(ids)

filteredE <- E[ids,]

postscript(file="results/heatmap.ps", horizontal=FALSE)
heatmap(filteredE[,])
dev.off()
