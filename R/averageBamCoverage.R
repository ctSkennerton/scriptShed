#averageBamCoverage <- function (bamfile) {
	# read in the bam file
#	bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
	# filter reads without match position
#	ind <- ! is.na(bam$pos)
	## remove non-matches, they are not relevant to us
#	bam <- lapply(bam, function(x) x[ind])
#	ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
	## names of the bam data frame:
	## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
	## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
	## construc: genomic ranges object containing all reads
#	ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
	## returns a coverage for each reference sequence (aka. chromosome) in the bam file
#	aveCov <- mean(coverage(ranges))
#	return(data.frame(name=names(aveCov), coverage=as.real(aveCov)))
#}

averageBamCoverage <- function(seqname, bamFile, ...) {
  	param <- ScanBamParam(what = c("pos", "qwidth"), which = GRanges(seqname, IRanges(1, 1e+07)), flag = scanBamFlag(isUnmappedQuery = FALSE))
 	x <- scanBam(bamFile, ..., param = param)[[1]]
 	mean(coverage(IRanges(x[["pos"]], width = x[["qwidth"]])))
}
