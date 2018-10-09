#' @title TODO
#' 
#' @description Compute the median TODO
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
#' @importFrom stats median
#' @keywords internal
smedian <- function(pos,v) {
    median(v[pos[1]:pos[2]],na.rm=T)
}
