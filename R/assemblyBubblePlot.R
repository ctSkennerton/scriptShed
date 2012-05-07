#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
option_list <- list(
		make_option(c("-b", "--bam"), default="NULL", type="character", 
					help="File name for the bam file"),
		make_option(c("-f", "--fasta"), default="NULL", type="character", 
					help="File name for the fasta file"),
		make_option(c("-o", "--output"), default="NULL", type="character", 
					help="Output file name for the image"),
		make_option(c("-p", "--legendPos"), default="bottomright", type="character", 
					help="position of the legend. choose from: bottomleft, bottomright (default), topleft, topright")			
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
suppressPackageStartupMessages(library(PBSmapping))
source("~/scripts/scriptShed/R/averageBamCoverage.R")
source("~/scripts/scriptShed/R/seqStats.R")
f<- read.fasta(options$fasta)

data <- vapply(f,seqStats, c(name="",length=0,gc=0.0))
data <- t(data)
coverage <- sapply(data[,1],averageBamCoverage, options$bam)

data = cbind(data, coverage)
data <- apply(data,1:2,as.numeric)

svg(options$output)

plot(c(min(data[,4]),max(data[,4])), c(min(data[,3]),max(data[,3])),col="white", xlab="coverage",ylab="gc")
# fix up the names of the columns as the addBubbles call requires them named this way
colnames(data) <- c("EID","Z","Y","X")
addBubbles(data, legend.title="Length", legend.pos=options$legendPos, symbol.bg=rgb(.9,.5,0,.6))
#symbols(data$coverage, data$gc, circles=data$length, fg="white", bg="red", inches=0.7, ylab="GC", xlab="coverage")
dev.off()
