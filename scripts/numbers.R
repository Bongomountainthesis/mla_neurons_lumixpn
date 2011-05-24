#!/usr/local/bin/Rscript

neurons.file <- "../mla_neurons_lumixpn/results/limma_results.csv"

neurons<- read.csv(neurons.file)

#remove genes above adj p val of 0.05

neurons.p<-neurons[neurons[,"adj.P.Val"]<=0.05,]

#order by FC and then remove duplicates
neurons.op <- neurons.p[order(abs(neurons.p[,"logFC"]),decreasing=TRUE),]

neurons.dop <- neurons.op[!duplicated(neurons.op[,"symbol"]),]

#write them to a csv file
neurons.dop <- neurons.dop[order(neurons.dop[,"logFC"],decreasing=TRUE),]
write.csv(neurons.dop, "results/significantunqiuegenes.csv")

#find genes that significantly increase
neurons.up <- neurons.dop[neurons.dop[,"logFC"]>=1,] 
neurons.up <- neurons.up[order(neurons.up[,"logFC"],decreasing=TRUE),]
write.csv(neurons.up, "results/significantuniqueUPgenes.csv")

#find genes that significantly decrease
neurons.down <- neurons.dop[neurons.dop[,"logFC"]<=-1,]
neurons.down <- neurons.down[order(neurons.down[,"logFC"], decreasing=FALSE),]
write.csv(neurons.down, "results/significantuniqueDOWNgenes.csv")



