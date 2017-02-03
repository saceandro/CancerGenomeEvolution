#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly=TRUE)

d <- read.delim(argv[1])

pdf(paste(argv[1], ".t1.rel.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,2]~d[,1],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Relative error of ", t[1])))
dev.off()

pdf(paste(argv[1], ".t2.rel.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,3]~d[,1],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Relative error of ", t[2])))
dev.off()

pdf(paste(argv[1], ".n1.rel.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,4]~d[,1],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Relative error of ", n[1])))
dev.off()

pdf(paste(argv[1], ".n2.rel.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,5]~d[,1],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Relative error of ", n[2])))
dev.off()

pdf(paste(argv[1], ".t.rel.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,6]~d[,1],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Relative error of ", t)))
dev.off()

pdf(paste(argv[1], ".n.rel.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,7]~d[,1],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Relative error of ", n)))
dev.off()

pdf(paste(argv[1], ".all.rel.pdf", sep=""), height=7, width=7)
par(ps = 24)
par(lwd = 2)
par(mex=1.5)
boxplot(d[,8]~d[,1],col="bisque",xlab=expression(paste("True ", t[2])), ylab=expression(paste("Relative error of ", "all parameters")))
dev.off()
