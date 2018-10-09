#' @title TODO
#' 
#' @description Compute the median absolute deviation TODO
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
#' @importFrom stats mad
#' @keywords internal
smad <- function(pos,v) {
    mad(v[pos[1]:pos[2]], na.rm=T)
}
