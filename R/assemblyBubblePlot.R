#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
option_list <- list(
		make_option(c("-b", "--bam"), default="NULL", type="character", 
					help="File name for the bam file"),
		make_option(c("-f", "--fasta"), default="NULL", type="character", 
					help="File name for the fasta file"),
		make_option(c("-o", "--output"), default="NULL", type="character", 
					help="Output file name for the image")
)
options <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

if( options$bam =="NULL") {
	print_help(OptionParser(option_list = option_list))
	q()
} 
if(options$fasta == "NULL") {
	print_help(OptionParser(option_list = option_list))
	q()
}
if(options$output == "NULL") {
	print_help(OptionParser(option_list = option_list))
	q()
}
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(Rsamtools))
library(plyr)
source("~/scripts/R/averageBamCoverage.R")
source("~/scripts/R/seqStats.R")
f<- read.fasta(options$fasta)

data <- lapply(f, seqStats)
data_as_frame <- ldply(data)

cvg <- lapply(data_as_frame$name,averageBamCoverage, options$bam)
names(cvg) <- data$name
cvg2 <- lapply(cvg, mean)


data2_as_frame <- data.frame(coverage=unlist(cvg2), name=data_as_frame$name)

c <- merge(data_as_frame, data2_as_frame,x.by="name", y.by="name")
svg(options$output)
symbols(c$coverage, c$gc, circles=c$length,fg="white",bg="red",inches=0.7,ylab="GC",xlab="coverage")
dev.off()
