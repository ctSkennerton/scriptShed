#!/usr/bin/env Rscript
###############################################################################
#
#    kmerFrequencies.R
#    
#    Caluclate the abundance of kmers for dna sequences
#
#    Copyright (C) 2012  Connor Skennerton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(seqinr))

option_list <- list(
		make_option(c("-i", "--input"), default="NULL", type="character", 
					help="File name for the bam file"),
		make_option(c("-k", "--kmer"), default="4", type="integer", 
					help="length of the kemrs"),
		make_option(c("-o", "--output"), default="NULL", type="character", 
					help="Output file name for the image")
)
options <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))

if( options$input =="NULL") {
	print_help(OptionParser(option_list = option_list))
	q()

f <- read.fasta(options$input)

#DF <- matrix(nrow=length(f),ncol=4*options$kmer) 
for(i in 1:length(f)){
	# append to the vectors
	count(f[[i]], options$kmer)
}
