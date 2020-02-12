#' @title TODO
#' 
#' @description TODO
#' 
#' @param logr A \code{vector} of \code{double} representing the log ratio
#' values.
#' 
#' @param emfit An object of class \code{Mclust} providing a 
#' mixture model estimation. 
#' 
#' @param zgroup A \code{matrix} of \code{double}, used as 
#' \code{integer}, with the number of rows corresponding to the current 
#' number of clusters while the number of columns is 
#' corresponding to the initial number of clusters. The presence of \code{1} 
#' in position \emph{[i,k]} indicates that the initial \emph{i}th cluster 
#' is now part of the new \emph{k}th cluster.
#' 
#' @param times A single \code{double} value specifying the number of 
#' time the median of each segment is sampled in order to predict the cluster 
#' assignment for the segment.
#' 
#' @return A \code{matrix} of \code{double} containing the probability for each
#' observation to belong to each class of the mixture model with the updated 
#' number of clusters.
#' 
#' @examples
#' 
#' # TODO
#' 
#' ## Load 
#' data(EMexample)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats predict
#' @keywords internal
getz <- function(logr, emfit, zgroup, times) { 
    
    ## Assign gz variable using mixture model when available
    ## When available, the matrix containing the probability for each
    ## observation to belong to each class is used (emfit$z)
    if (is.null(emfit$z)) {
        gz <- matrix(ncol = 1, data = rep(1, length(logr)))
    } else {
        gz <- predict(emfit, newdata=logr)$z
    }
    
    
    ## Identify positions where probability is not finite and assign zero
    isfin <- matrix(ncol = ncol(gz), nrow = nrow(gz), data = is.finite(gz))
    gz[!isfin] <- 0 #just being honest: we don't know how to assign these
    
    ## Sum the probability for each time the model is run
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
