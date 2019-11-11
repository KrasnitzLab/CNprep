#' @title TODO
#' 
#' @description TODO
#' 
#' @param logr TODO a \code{vector} of \code{numeric}
#' 
#' @param emfit An object of class \code{Mclust} providing a 
#' mixture model estimation. 
#' 
#' @param times A single \code{double} value specifying the number of 
#' time the median of each segment is sampled in order to predict the cluster 
#' assignment for the segment.
#' 
#' @param center TODO
#' 
#' @return TODO a \code{vector} of \code{numeric}
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats pnorm
#' @keywords internal
centerprob <- function(logr, emfit, zgroup, times, center) { 
    
    ## Create new matrix using input data
    gz <- matrix(nrow = length(logr), ncol = length(emfit$parameters$mean), 
                    data=logr)
    
    #gz<-t(emfit$pro*exp(-0.5*(t(gz)-emfit$mu)^2/emfit$sigma)/
    #sqrt(2*pi*emfit$sigma))
    
    ## Assign correct value to epro object according to mean length
    if (length(emfit$parameters$mean) == 1) { 
        epro <- as.vector(1)
    }
    
    if (length(emfit$parameters$mean) > 1) {
        epro <- emfit$parameters$pro        
    }
    
    gz <- t(epro * pnorm(-abs(t(gz) - emfit$parameters$mean) /
                sqrt(emfit$parameters$variance$sigmasq)))
    
    ## Combine columns of z table using indicator matrix zgroup
    gz <- gz %*% t(zgroup)  
    
    gz <- matrix(ncol = ncol(gz),
            data = apply(gz, 2, cumsum)[seq(from = times, 
                                    to = nrow(gz), by = times),] / times)
    
    ## Mean value within each segment
    gz <- (gz[,center] - c(0, gz[-nrow(gz), center])) / 
                            sum(epro[zgroup[center,] == 1])
    
    ## Ensure that the return value is always inferior to 0.5
    return(ifelse(gz < 0.5, gz, 1 - gz))
}
