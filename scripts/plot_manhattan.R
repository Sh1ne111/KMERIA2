#!/usr/bin/env Rscript



# Usage: Rscript plot_manhattan.R sigkmer.plot.txt sigkmer 1e-7


argv<-commandArgs(TRUE)

library('CMplot')

a <- read.table(argv[1], header=F)

## rgb(129,196,240)  skyblue
## rgb(235,70,144  pink

CMplot(a,plot.type = 'm',file='png',file.name=argv[2],col=c(rgb(192,192,192,max=255),rgb(0,0,0,max=255)),threshold = as.numeric(argv[3]),threshold.col='red',threshold.lty = 2,threshold.lwd = 1, amplify = TRUE, signal.cex = c(1,1), signal.pch = c(20,20),signal.col = NULL, dpi=326,cex=0.6,conf.int=T,conf.int.col='grey',box=F,main='',ylab=expression(-log[10](italic(p))), width=12, height = 4)
