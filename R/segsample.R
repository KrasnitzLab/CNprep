#' @title TODO
#' 
#' @description TODO
#' 
#' @param mysegs TODO
#' 
#' @param ratcol TODO
#' 
#' @param startcol TODO
#' 
#' @param endcol TODO
#' 
#' @param blocksize TODO
#' 
#' @param times TODO
#' 
#' @return TODO
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
segsample <-
function(mysegs,ratcol,startcol="StartProbe",endcol="EndProbe",
	blocksize=0,times=0){
	if(blocksize==0&times==0)stop("One of blocksize or times must be set")
	if(blocksize!=0&times!=0)stop("Only one of blocksize or times can be set")
	segtable<-mysegs[,c(startcol,endcol),drop=F]
	## Comment Pascal: at least one result should be different from zero
	if(blocksize!=0)segtable<-
		segtable[rep(1:nrow(segtable),
		times=(segtable[,endcol]-segtable[,startcol]+1)%/%blocksize),]
	if(times!=0)segtable<-segtable[rep(1:nrow(segtable),each=times),]
	return(cbind(segtable, apply(segtable, 1, smedian.sample, v = ratcol)))
}
