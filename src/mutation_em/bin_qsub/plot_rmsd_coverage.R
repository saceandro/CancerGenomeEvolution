#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly=TRUE)

pdf(paste(argv[1], ".u.pdf", sep=""), height=7, width=7)
d <- read.delim(argv[1])
## par(yaxs="i")
## par(mar = c(5.5, 8.0, 4.1, 2)) #  余白の広さを行数で指定．下，左，上，右の順．
## par(mgp = c(4, 1.2, 0))        #  余白の使い方．説明，ラベル，軸の位置を行で指定．
## par(cex.lab=3)
## par(cex.axis=1.6)
boxplot(d[,3]~d[,1],col="bisque",xlab="Sequencing depth", ylab="RMSD")
dev.off()

pdf(paste(argv[1], ".t.pdf", sep=""), height=7, width=7)
boxplot(d[,4]~d[,1],col="bisque",xlab="Sequencing depth", ylab="RMSD")
dev.off()

pdf(paste(argv[1], ".n.pdf", sep=""), height=7, width=7)
boxplot(d[,5]~d[,1],col="bisque",xlab="Sequencing depth", ylab="RMSD")
dev.off()
