#' @title TODO
#' 
#' @description TODO
#' 
#' @param segrat a \code{list}
#' 
#' @param blsize A single \code{integer} specifying the bootstrap 
#' sampling rate of segment medians to generate input for model-based 
#' clustering. The number of times a segment is sampled is then given by the 
#' (integer) division of the segment length in internal units by \code{blsize}.
#' 
#' @param minjoin A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the degree of overlap above which two clusters will be joined 
#' into one.
#' 
#' @param ntrial A single \code{integer} specifying the number of times 
#' a model-based clustering is attempted for each profile in order to 
#' achieve the highest Bayesian information criterion (BIC).
#' 
#' @param bestbic A single \code{numeric} value for initalizing BIC 
#' maximization. A large negative value is recommended. The default 
#' is \code{-1e7}.
#' 
#' @param modelNames A \code{vector} of \code{character} strings specifying 
#' the names of models to be used in model-based clustering (see package 
#' \code{mclust} for further details). The default is \code{"E"}.
#' 
#' @param cweight A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the minimal share of the central cluster in each profile.
#' 
#' @param bstimes A single \code{double} value specifying the number of 
#' time the median of each segment is sampled in order to predict the cluster 
#' assignment for the segment.
#' 
#' @param chromrange A \code{integer} vector enumerating chromosomes from 
#' which segments are to be used for initial model-based clustering.
#' 
#' @return TODO
#'
#' @examples
#'
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats median var
#' @importFrom mclust Mclust
#' @keywords internal
CNclusterNcenter <- function(segrat, blsize, minjoin, ntrial, bestbic,
    modelNames, cweight, bstimes, chromrange) {

    ## Create and assign seeds
    # .lec.CreateStream(segrat$stream)
    # .lec.SetSeed(segrat$stream,seedme)
    # for( j in seq_len(segrat$sub)) .lec.ResetNextSubstream(segrat$stream)
    #     .lec.CurrentStream(segrat$stream)

    startcol <- "StartProbe"
    endcol   <- "EndProbe"
    chromcol <- "chrom"
    medcol   <- "segmedian"
    madcol   <- "segmad"

    segrat$seg<-cbind(segrat$seg,
                        t(apply(segrat$seg[,c(startcol, endcol, chromcol),
                        drop=FALSE], 1, smedmad, v=segrat$rat)))

    dimnames(segrat$seg)[[2]] <- c(startcol, endcol, chromcol, medcol, madcol)
    
    seguse <- segrat$seg[segrat$seg[, chromcol] %in% chromrange,, drop=FALSE]
    
    aux <- rep(0, length(segrat$rat))
    aux[seguse[, startcol]] <- 1
    aux[seguse[, endcol]]   <- (-1)
    aux <- cumsum(aux)
    aux[seguse[, endcol]]   <- 1
    
    ratuse <- segrat$rat[aux == 1]

    for(j in seq_len(ntrial)) {
        aaa <- segsample(seguse, ratuse, blocksize=blsize)
        if (all(unique(aaa[,3]) == 0)) { 
            aaa[,3] <- 1e-10 
        }
        emfit <- Mclust(aaa[,3], maxG=15, modelNames=modelNames)
        if (emfit$bic >= bestbic) {
            bestaaa <- aaa
            bestem  <- emfit
            bestbic <- emfit$bic
        }
    }

    newem <- consolidate(bestem, minjoin)
    newem <- get.center(newem, cweight)
    
    if (length(bestem$parameters$mean) == 1) { 
        profcenter <- median(bestaaa[,3]) 
    } else { 
        profcenter <- weighted.median(bestaaa[,3], newem$z[, newem$center]) 
    }

    mediandev <- segrat$seg[, medcol] - profcenter
    segs <- segsample(segrat$seg, segrat$rat, times = bstimes)

    if (all(unique(aaa[,3]) == 1e-10)) { 
        segs[segs[,3] == 0, 3] <- 1e-10 
    }

    segzall <- getz(segs[,3], bestem, newem$groups, times=bstimes)
    centerz <- segzall[, newem$center]
    maxz <- segzall[nrow(segzall) * (max.col(segzall) - 1) + 
                            (seq_len(nrow(segzall)))]
    
    maxzcol <- max.col(segzall)
    maxzmean <- newem$mu[maxzcol] - newem$mu[newem$center]
    maxzsigma <- sqrt(newem$sigmasq[maxzcol])
    cpb <- centerprob(segs[,3], bestem, newem$groups, times=bstimes, 
                newem$center)
    w <- t(matrix(nrow=bstimes, data=segs[,3]))
    segerr <- sqrt(apply(w, 1, var, na.rm=TRUE))
    
    # .lec.CurrentStreamEnd()
    # .lec.DeleteStream(segrat$stream)
        
    return(cbind(segrat$seg[, medcol], segrat$seg[,madcol], mediandev, segerr, 
                    centerz, cpb, maxz, maxzmean, maxzsigma))
}
