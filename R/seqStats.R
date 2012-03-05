seqStats <- function(fastaObj){
	return(c(attr(fastaObj,"name"),length(fastaObj),GC(fastaObj)))
}