#' @title TODO
#' 
#' @description TODO
#' 
#' @param normalmedian TODO
#' 
#' @param normalength TODO
#' 
#' @param tumormedian TODO
#' 
#' @param tumorlength TODO. Default: \code{NULL}.
#' 
#' @param normalmad TODO. Default: \code{NULL}.
#' 
#' @param normalerror TODO. Default: \code{NULL}.
#' 
#' @param tumormad TODO. Default: \code{NULL}.
#' 
#' @param tumorerror TODO. Default: \code{NULL}.
#' 
#' @return TODO
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
normalComparison <- function(normalmedian, normalength, tumormedian, 
                             tumorlength, normalmad=NULL, normalerror=NULL, 
                             tumormad=NULL, tumorerror=NULL){
	nisnull <- c(!(is.null(normalmad)|is.null(tumormad)),
		!(is.null(normalerror)|is.null(tumorerror)))
	nsred <- matrix(ncol=2+sum(nisnull), length(normalmedian),
		data=c(normalength, normalmedian,
		if(nisnull[1])normalmad,if(nisnull[2])normalerror),dimnames=list(NULL,
		c("length","mediandev",if(nisnull[1])"segmad",if(nisnull[2])"segerr")))
	nsred <- nsred[order(nsred[,"mediandev"]),,drop=FALSE]
	lnorm <- sum(nsred[,"length"])
	z <- cbind(c(nsred[,"mediandev"],tumormedian),c(nsred[,"length"],
		rep(0,length(tumormedian))),c(rep(0,nrow(nsred)),1:length(tumormedian)))
	z <- z[order(z[,1]),,drop=FALSE]
	z[,2] <- cumsum(z[,2])/lnorm
	z <- z[z[,3]!=0,,drop=FALSE]
	negtail <- z[order(z[,3]),2]
	if(nisnull[1]){
		z <- cbind(c(nsred[,"mediandev"]/nsred[,"segmad"],
  		tumormedian/tumormad),c(nsred[,"length"],
  		rep(0,length(tumormedian))),
			c(rep(0,nrow(nsred)),1:length(tumormedian)))
		z <- z[order(z[,1]),,drop=FALSE]
		z[,2] <- cumsum(z[,2])/lnorm
		z <- z[z[,3]!=0,,drop=FALSE]
		negtailnormad<-z[order(z[,3]),2]
	}
	if(nisnull[2]){
		z <- cbind(c(nsred[,"mediandev"]/nsred[,"segerr"],tumormedian/tumorerror),
			c(nsred[,"length"],rep(0,length(tumormedian))),
			c(rep(0,nrow(nsred)),1:length(tumormedian)))
		z <- z[order(z[,1]),,drop=FALSE]
		z[,2] <- cumsum(z[,2])/lnorm
		z <- z[z[,3]!=0,,drop=FALSE]
		negtailnormerror <- z[order(z[,3]),2]
	}
	return(matrix(ncol=2+sum(nisnull),data=c(lnorm%/%tumorlength,
		negtail,if(nisnull[1])negtailnormad,if(nisnull[2])negtailnormerror),
		dimnames=list(NULL,c("samplesize", "negtail",if(nisnull[1])"negtailnormad",
		if(nisnull[2])"negtailnormerror"))))
}
