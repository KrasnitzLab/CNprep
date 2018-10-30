#' @title TODO
#' 
#' @description TODO
#' 
#' @param emfit a \code{list} containing information about the updated mixture
#' model for the clusters as given by an object of class 
#' \code{Mclust}.  
#' 
#' @param mincenter TODO
#' 
#' @return a \code{list} containing information about the updated mixture
#' model for the clusters: 
#' \itemize{
#' \item{ \code{mu}: a \code{vector} of \code{double} which are the mean 
#' of each cluster in the model.}
#' \item{ \code{pro}: a \code{vector} of \code{double} which are the mixing 
#' proportion fo each cluster in the model.}
#' \item{ \code{z}: a \code{matrix} of \code{double} whose \emph{[i,k]}th 
#' entry is
#' the probability that observation \emph{i} in the dataset belongs to the 
#' \emph{k}th cluster.}
#' \item{ \code{groups}: a \code{matrix}  }
#' \item{ \code{ngroups}: a \code{integer} specifying the number of clusters
#' in the model. }
#' \item{ \code{sigmasq}: a \code{vector} of \code{double} whose \emph{k}th 
#' component is 
#' the variance for the \emph{k}th cluster in the mixture model.} 
#' 
#' @details The \code{emfit} parameter correspond to an object of class 
#' \code{Mclust} that has been modified by the \code{CNprep:::consolidate}
#' function.
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
    newem <- list(mu = emfit$mu, pro = emfit$pro, z = emfit$z, 
                    groups = emfit$groups, ngroups = emfit$ngroups, 
                    sigmasq = emfit$sigmasq, center = which.min(abs(emfit$mu)))
    
    ## Only process when the number of cluster is superior to one
    while(newem$ngroups > 1) {
        omu <- order(abs(newem$mu))
        ## Group clusters when the mixing proportion is inferior to the
        ## specified minimal share of the central cluster
        if (newem$pro[omu[1]] < mincenter) {
            gl <- min(omu[seq_len(2)])
            gr <- max(omu[seq_len(2)])
            newem$z[, gl] <- newem$z[, gl] + newem$z[, gr]
            newem$z <- newem$z[, -gr, drop=FALSE]
            numu <- (newem$mu[gl]*newem$pro[gl]+newem$mu[gr]*newem$pro[gr])/
                            (newem$pro[gl]+newem$pro[gr])
            newem$sigmasq[gl] <- (newem$pro[gl]*(newem$sigmasq[gl] + 
                    newem$mu[gl]^2) + newem$pro[gr]*(newem$sigmasq[gr] + 
                                                        newem$mu[gr]^2))/
                    (newem$pro[gl]+newem$pro[gr])-numu^2
            newem$mu[gl] <- numu
            newem$mu <- newem$mu[-gr]
            newem$sigmasq <- newem$sigmasq[-gr]
            newem$pro[gl] <- newem$pro[gl] + newem$pro[gr]
            newem$pro <- newem$pro[-gr]
            newem$groups[gl,] <- newem$groups[gl,] + newem$groups[gr,]
            newem$groups <- newem$groups[-gr,,drop=FALSE]
            newem$ngroups <- newem$ngroups - 1
            newem$center <- gl
        }
        else break
    }
    return(newem)
}
