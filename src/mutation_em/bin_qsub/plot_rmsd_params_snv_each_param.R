#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly=TRUE)

d <- read.delim(argv[1])

pdf(paste(argv[1], ".t1.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,3]~d[,2],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Absolute error of ", t[1])))
dev.off()

pdf(paste(argv[1], ".t2.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,4]~d[,2],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Absolute error of ", t[2])))
dev.off()

pdf(paste(argv[1], ".n2.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,6]~d[,2],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Absolute error of ", n[2])))
dev.off()
