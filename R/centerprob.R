#' @title TODO
#' 
#' @description TODO
#' 
#' @param logr a \code{vector} of \code{double} TODO
#' 
#' @param emfit an object of class \code{Mclust} providing a 
#' mixture model estimation for the clusters. 
#' 
#' @param times a \code{integer} TODO
#' 
#' @param center a \code{integer} TODO
#' 
#' @return TODO
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats pnorm
#' @keywords internal
centerprob <- function(logr, emfit, zgroup, times, center) { 
    ## Create a matrix with one column per cluster
    gz <- matrix(nrow=length(logr), ncol=length(emfit$parameters$mean), 
                    data=logr)
    #gz<-t(emfit$pro*exp(-0.5*(t(gz)-emfit$mu)^2/emfit$sigma)/
    #sqrt(2*pi*emfit$sigma))
    if (length(emfit$parameters$mean) == 1) {
        epro <- as.vector(1)
    }
    if (length(emfit$parameters$mean) > 1) {
        epro <- emfit$parameters$pro        
    }
    
    gz <- t(epro * pnorm(-abs(t(gz) - emfit$parameters$mean)/
                        sqrt(emfit$parameters$variance$sigmasq)))
    gz <- gz%*%t(zgroup) # combine columns of z table using indicator matrix zgroup 
    gz <- matrix(ncol=ncol(gz),
        data=apply(gz,2,cumsum)[seq(from=times,to=nrow(gz),by=times),]/times)
    #mean value within each segment
    gz <- (gz[,center] - c(0, gz[-nrow(gz), center])) / 
                                sum(epro[zgroup[center,] == 1]) 
    
    ## Ensure that returned values that are inferior to 0.5
    return(ifelse(gz < 0.5, gz, 1 - gz))
}
