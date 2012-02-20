seqStats <- function(fastaObj){
	return(data.frame(name=attr(fastaObj,"name"),length=length(fastaObj),gc=GC(fastaObj)))
}