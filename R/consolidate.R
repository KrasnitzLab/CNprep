#' @title Join clusters until the minimum degree of overlap is reached
#' 
#' @description The function join clusters that have the minimum ratio of
#' overlap as specified by user.
#' 
#' @param emfit an object of class \code{Mclust} providing a 
#' mixture model estimation. 
#' 
#' @param minover a single \code{numeric} value between \code{0} and \code{1} 
#' specifying the degree of overlap above which two clusters will be joined 
#' into one.
#' 
#' @return a \code{list} containing information about the updated
#' clusters obtained from a mixture model estimation:
#' \itemize{
#' \item \code{mu} a \code{numeric} \code{vector} representing the mean 
#'     for each component. If there is more than one component, the 
#'     \emph{k}th element is the mean of the \emph{k}th component of the 
#'     mixture model.
#' \item \code{pro} a \code{vector} whose \emph{k}th component is the mixing 
#'     proportion for the \emph{k}th component of the mixture model. If 
#'     missing, equal proportions are assumed.
#' \item \code{z} a \code{numeric} \code{matrix} whose \emph{[i,k]}th entry 
#'     is the probability that observation \emph{i} in the test data belongs 
#'     to the \emph{k}th class.
#' \item \code{groups} a \code{matrix} of \code{double}, used as 
#'     \code{integer}, with the number of rows corresponding
#'     to the current number of clusters while the number of columns is 
#'     corresponding to the initial number of clusters. The presence of 
#'     \code{1} in position \emph{[i,k]} indicates that the initial 
#'     \emph{i}th cluster is now part of the new \emph{k}th cluster.
#' \item \code{ngroups} a \code{numeric}, used as an integer, giving the final
#'     number of clusters.
#' \item \code{sigmasq}  a \code{numeric} \code{vector} giving the common 
#'     variance for each component in the mixture model "E".
#' }
#' 
#' @examples
#'
#' ## Load Mclust object
#' data(EMexample)
#' 
#' ## Group clusters that have at least 0.4% of overlap
#' ## The inital object has 5 clusters while the return object has only 
#' ## 4 clusters
#' CNprep:::consolidate(EMexample, minover=0.004)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom mclust Mclust
#' @keywords internal
consolidate <- function(emfit, minover) { 
    
    ## Create new object using input parameter
    newem <- list(mu=emfit$parameters$mean, 
                    pro=emfit$parameters$pro,
                    z=emfit$z,
                    groups=matrix(nrow=length(emfit$parameters$mean),
                                ncol=length(emfit$parameters$mean), data=0),
                    ngroups=length(emfit$parameters$mean),
                    sigmasq=emfit$parameters$variance$sigmasq)
    
    ## Ensure that z entry is a matrix in the newem object
    if (is.null(emfit$z)) {
        newem$z <- matrix(ncol=1, data=rep(1, emfit$n))
    }
    
    ## Ensure that sigmasq entry has the same length than mu in newem object
    if (length(newem$sigmasq) == 1) {
        newem$sigmasq <- rep(newem$sigmasq, length(emfit$parameters$mean))
    }
    
    diag(newem$groups) <- 1
    
    while(newem$ngroups > 1) {
        #note the asymmetry; a fraction of me in you != a fraction of you in me
        avz <- (t(newem$z)/colSums(newem$z)) %*% newem$z 
        diag(avz) <- 0
        if (max(avz) > minover) {
            g1 <- col(avz)[which.max(avz)]
            g2 <- row(avz)[which.max(avz)]
            gl <- min(g1, g2)
            gr <- max(g1, g2)
            newem$z[,gl] <- newem$z[,gl] + newem$z[,gr]
            newem$z <- matrix(ncol=ncol(newem$z) - 1, 
                                nrow=nrow(newem$z),
                                data=newem$z[, -gr])
            
            numu <- (newem$mu[gl]*newem$pro[gl]+newem$mu[gr]*newem$pro[gr])/
                            (newem$pro[gl]+newem$pro[gr])
            
            newem$sigmasq[gl] <- (newem$pro[gl] * 
                                    (newem$sigmasq[gl] + newem$mu[gl]^2) +
                                    newem$pro[gr] * 
                                    (newem$sigmasq[gr] + newem$mu[gr]^2))/
                                    (newem$pro[gl] + newem$pro[gr]) - numu^2
            
            newem$mu[gl] <- numu
            newem$mu <- newem$mu[-gr]
            newem$sigmasq <- newem$sigmasq[-gr]
            newem$pro[gl] <- newem$pro[gl] + newem$pro[gr]
            newem$pro <- newem$pro[-gr]
            newem$groups[gl,] <- newem$groups[gl,] + newem$groups[gr,]
            newem$groups <- matrix(ncol = ncol(newem$groups),
                                        nrow = newem$ngroups - 1,
                                        data = newem$groups[-gr,])
            newem$ngroups <- newem$ngroups - 1
        }
        else break
    }
    
    return(newem)
}
