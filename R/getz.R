#' @title TODO
#' 
#' @description TODO
#' 
#' @param logr TODO
#' 
#' @param emfit An object of class \code{Mclust} providing a 
#' mixture model estimation. 
#' 
#' @param zgroup TODO
#' 
#' @param times A single \code{double} value specifying the number of 
#' time the median of each segment is sampled in order to predict the cluster 
#' assignment for the segment.
#' 
#' @return TODO
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats predict
#' @keywords internal
getz <- function(logr, emfit, zgroup, times) { 
    
    ## Assign gz variable using mixture model when available
    if (is.null(emfit$z)) {
        gz <- matrix(ncol = 1, data = rep(1, length(logr)))
    } else {
        gz <- predict(emfit, newdata=logr)$z
    }
    
    isfin <- matrix(ncol = ncol(gz), nrow = nrow(gz), data = is.finite(gz))
    gz[!isfin] <- 0 #just being honest: we don't know how to assign these
    gz <- matrix(ncol = ncol(gz),
        data = apply(gz, 2, cumsum)[seq(from = times, 
                                            to = nrow(gz), by = times),])
    
    cisfin <- matrix(ncol = ncol(isfin), 
                    data = apply(isfin, 2, cumsum)[seq(from = times, 
                                    to = nrow(isfin), by = times),])
    cisfin <- cisfin -
        rbind(matrix(nrow = 1, data = rep(0, ncol(cisfin))),
                    cisfin[-nrow(gz),,drop = FALSE])
    
    gz <- (gz - rbind(matrix(nrow = 1, data = rep(0, ncol(gz))), 
                        gz[-nrow(gz),,drop = FALSE])) / times
    gz[cisfin < times] <- NA
    
    return(gz %*% t(zgroup))
}
