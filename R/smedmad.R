#' @title TODO
#' 
#' @description Compute the median and the median absolute deviation,TODO
#' 
#' @param pos TODO
#' 
#' @param v TODO
#' 
#' @return TODO
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats median mad
#' @keywords internal
smedmad <-
function(pos,v)
	c(median(v[pos[1]:pos[2]],na.rm=T),mad(v[pos[1]:pos[2]],na.rm=T))
