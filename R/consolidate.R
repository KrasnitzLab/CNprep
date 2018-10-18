#' @title Extract, modify and return information from a provided mixture model 
#' estimation for the clusters.
#' 
#' @description Extract information from the provided mixture model 
#' estimation for the clusters and merge clusters that exceed a specified 
#' overlap value. The udpated information about the model is returned.
#' 
#' @param emfit An object of class \code{Mclust} providing a 
#' mixture model estimation for the clusters. 
#' 
#' @param minover A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the degree of overlap above which two clusters will be joined 
#' into one.
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
#' }
#'
#' @examples
#'
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom mclust Mclust
#' @keywords internal
consolidate <- function(emfit, minover) { 
    newem <- list(mu=emfit$parameters$mean, pro=emfit$parameters$pro,
                z=emfit$z, groups=matrix(nrow=length(emfit$parameters$mean),
                ncol=length(emfit$parameters$mean), data=0),
                ngroups=length(emfit$parameters$mean),
                sigmasq=emfit$parameters$variance$sigmasq)
    if (is.null(emfit$z)) {
        newem$z<-matrix(ncol=1, data=rep(1, emfit$n))
    }
    
    if (length(newem$sigmasq) == 1) {
        newem$sigmasq <- rep(newem$sigmasq, length(emfit$parameters$mean))
    }
    
    diag(newem$groups) <- 1
    
    while(newem$ngroups > 1) {
        #note the asymmetry; a fraction of me in you != a fraction of you in me
        avz <- (t(newem$z)/colSums(newem$z))%*%newem$z 
        diag(avz) <- 0
        if (max(avz) > minover) {
            g1 <- col(avz)[which.max(avz)]
            g2 <- row(avz)[which.max(avz)]
            gl <- min(g1, g2)
            gr <- max(g1, g2)
            newem$z[,gl] <- newem$z[,gl] + newem$z[,gr]
            newem$z <- matrix(ncol=ncol(newem$z)-1, nrow=nrow(newem$z),
                                data=newem$z[,-gr])
            numu <- (newem$mu[gl]*newem$pro[gl] + newem$mu[gr]*newem$pro[gr])/
                            (newem$pro[gl]+newem$pro[gr])
            newem$sigmasq[gl] <- (newem$pro[gl]*(newem$sigmasq[gl] + 
                                                        newem$mu[gl]^2) +
                    newem$pro[gr]*(newem$sigmasq[gr] + newem$mu[gr]^2))/
                    (newem$pro[gl] + newem$pro[gr])-numu^2
            newem$mu[gl] <- numu
            newem$mu <- newem$mu[-gr]
            newem$sigmasq <- newem$sigmasq[-gr]
            newem$pro[gl] <- newem$pro[gl] + newem$pro[gr]
            newem$pro <- newem$pro[-gr]
            newem$groups[gl,] <- newem$groups[gl,] + newem$groups[gr,]
            newem$groups <- matrix(ncol=ncol(newem$groups),
                                    nrow=newem$ngroups-1,
                                    data=newem$groups[-gr,])
            newem$ngroups <- newem$ngroups - 1
        }
        else break
    }
    return(newem)
}
