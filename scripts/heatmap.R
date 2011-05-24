#!/usr/local/bin/Rscript

library(beadarray)
library(gplots)

load("results/BSData.quantile.RData")

limma <- read.csv("results/significantuniquegenes.csv")

E <- exprs(BSData.quantile)

inds1 <- which(abs(limma[,"logFC"])>2.763)
inds2 <-which(limma[,"adj.P.Val"]<0.05)

inds <-intersect(inds1, inds2)

ids <- limma[inds,"ID"]

ids <- as.character(ids)

filteredE <- E[ids,]

symbol <- limma[inds,"symbol"]


postscript(file="results/heatmap.ps", horizontal=FALSE)
heatmap.2(filteredE[,],
		Colv=NA,
		col=topo.colors(75), 
		scale="none",
		key=TRUE,
		keysize=0.75,
		symkey=FALSE,
		density.info="none",
		trace="none", 
		labRow=symbol,
		cexRow=0.75,
	)
dev.off()

