#!/usr/local/bin/Rscript


###load limma results file and nearest or overlapping feature data and merge them together

limma.filename <- "/space/matt/mla_neurons_lumixpn/results/limma_results.csv"
NS.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H3K4me3/nearest_peak_to_gene.csv"
NS.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_MLA_H327me3/nearest_peak_to_gene.csv"
N.K4.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_Neuron_H3K4me3/nearest_peak_to_gene.csv"
N.K27.filename <- "/space/matt/MLA_neuron_ChIPSeq/results/SICER_NeuronMLA_H3K27me3/nearest_peak_to_gene.csv"

limma <- read.csv(limma.filename)
NS.K4 <- read.csv(NS.K4.filename)
NS.K27 <- read.csv(NS.K27.filename)
N.K4 <- read.csv(N.K4.filename)
N.K27 <- read.csv(N.K274.filename)

##remove duplicate probe-gene matches - keep highest FC difference
limma.o <- limma[order(abs(limma[,"logFC"]), decreasing=TRUE),]
limma.do <- limma.o[!duplicated(limma.o[,"EnsemblID"]),]

##calculate which genes have a 'real' peak
#for K4 cut off at 1kb and neg10logpVal at 250
#for K27 cut off at 3kb
plot(NS.K4[,"distancetoFeature"],NS.K4[,"neg10log10pVal"], xlim=c(-2000,2000), pch=".")

NS.K4.dist <- NS.K4[abs(NS.K4[,"distancetoFeature"])<=1000,]
NS.K4.sig <- NS.K4.dist[NS.K4.dist[,"neg10log10pVal"]>=250,]

##merge columns by ensembl id
dat <- merge(limma.do, NS.K4.sig, by.x="EnsemblID", by.y="feature", all.x=T)
write.csv(dat, "/space/matt/MLA_neuron_ChIPSeq/results/expression_histones_all.csv")


#remove all the spare columns
rm.cols<-c("ensembl_gene_id", "AveExpr", "t", "P.Value", "B", "space", "names", "start_position.y", "end_postition.y", "peak.id", "nTags_ChIP", "nTags_Cnt", "mgi_symbol")
keep.cols<-!(colnames(dat) %in% rm.cols)

dat.tidy <- dat[,keep.cols]
