#' @title TODO
#' 
#' @description TODO
#' 
#' @param emfit An object of class \code{Mclust} providing a 
#' mixture model estimation. 
#' 
#' @param mincenter TODO
#' 
#' @return TODO A \code{list} containing:
#' \itemize{
#' \item \code{mu} TODO}
#' \itme \code{pro} A vector whose \emph{k}th component is the mixing 
#' proportion for the \emph{k}th component of the mixture model. If missing, 
#' equal proportions are assumed.
#' \item \code{z} A matrix whose \emph{[i,k]}th entry is the probability that 
#' observation \emph{i} in the test data belongs to the \emph{k}th class.
#' \item \code{groups} TODO
#' \item \code{ngroups} TODO
#' \item \code{sigmasq} TODO
#' \item \code{center} TODO
#' }
#' 
#' @examples
#'
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
get.center <- function(emfit, mincenter) 
{
    #emfit must have come out of consolidate! 
    newem<-list(mu=emfit$mu, pro=emfit$pro, z=emfit$z, 
                    groups=emfit$groups, ngroups=emfit$ngroups, 
                    sigmasq=emfit$sigmasq, center=which.min(abs(emfit$mu)))
    
    while (newem$ngroups > 1) {
        omu <- order(abs(newem$mu))
        if (newem$pro[omu[1]] < mincenter) {
            gl <- min(omu[seq_len(2)])
            gr <- max(omu[seq_len(2)])
            newem$z[,gl] <- newem$z[, gl] + newem$z[, gr]
            newem$z <- newem$z[, -gr, drop=FALSE]
            numu <- (newem$mu[gl] * newem$pro[gl] + 
                            newem$mu[gr] * newem$pro[gr])/
                            (newem$pro[gl] + newem$pro[gr])
            newem$sigmasq[gl] <- (newem$pro[gl] * (newem$sigmasq[gl] + 
                    newem$mu[gl]^2) + newem$pro[gr] * (newem$sigmasq[gr] + 
                                                        newem$mu[gr]^2))/
                    (newem$pro[gl]+newem$pro[gr])-numu^2
            newem$mu[gl] <- numu
            newem$mu <- newem$mu[-gr]
            newem$sigmasq <- newem$sigmasq[-gr]
            newem$pro[gl] <- newem$pro[gl] + newem$pro[gr]
            newem$pro <- newem$pro[-gr]
            newem$groups[gl,] <- newem$groups[gl,] + newem$groups[gr,]
            newem$groups <- newem$groups[-gr, , drop=FALSE]
            newem$ngroups <- newem$ngroups - 1
            newem$center <- gl
        }
        else break
    }
    return(newem)
}
