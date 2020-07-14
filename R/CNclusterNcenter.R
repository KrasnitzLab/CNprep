#' @title TODO
#' 
#' @description TODO
#' 
#' @param segrat a \code{list} containing information for one specific profile:
#' \itemize{
#' \item{\code{seg}}{ a \code{data.frame} with 3 columns:} \itemize{
#'     \item{\code{StartProbe}}{ a \code{numeric} that tabulates the (integer) 
#'     start position of each segment in internal units such as probe numbers.}
#'     \item{\code{EndProbe}}{ a \code{numeric} that tabulates the (integer) 
#'     end position of each segment in internal units such as probe numbers.}
#'     \item{\code{chrom}}{ a \code{numeric} representing the chromosome.}}
#' \item{\code{rat}}{ a \code{numeric} \code{vector} of elements that are 
#'     functions of copy number, most often log ratios of copy number to 
#'     the expected standard value, such as 2 in diploid genomes for a specific
#'     profile. }
#' \item{\code{stream}}{ a \code{character} \code{string} representing the 
#'     profile ID.}
#' \item{\code{sub}}{ a \code{numeric} representing the position of the current
#'     profile ID in a \code{vector} of profiles.}
#' \item{\code{weight}}{ a \code{numeric} \code{vector} of elements that are 
#' weight the bin} 
#' }
#' 
#' @param blsize a single \code{integer} specifying the bootstrap 
#' sampling rate of segment medians to generate input for model-based 
#' clustering. The number of times a segment is sampled is then given by the 
#' (integer) division of the segment length in internal units by \code{blsize}.
#' 
#' @param minJoin a single \code{numeric} value between \code{0} and \code{1} 
#' specifying the degree of overlap above which two clusters will be joined 
#' into one.
#' 
#' @param nTrial a single \code{integer} specifying the number of times 
#' a model-based clustering is attempted for each profile in order to 
#' achieve the highest Bayesian information criterion (BIC).
#' 
#' @param bestBIC a single \code{numeric} value for initializing BIC 
#' maximization. A large negative value is recommended. 
#' 
#' @param modelNames a \code{vector} of \code{character} strings specifying 
#' the names of models to be used in model-based clustering (see package 
#' \code{mclust} for further details).
#' 
#' @param cweight a single \code{numeric} value between \code{0} and \code{1} 
#' specifying the minimal share of the central cluster in each profile.
#' 
#' @param bstimes a single \code{double} value specifying the number of 
#' time the median of each segment is sampled in order to predict the cluster 
#' assignment for the segment.
#' 
#' @param keepClust a single \code{logical} if the mClust object is keep 
#' 
#' @param chromRange a \code{integer} vector enumerating chromosomes from 
#' which segments are to be used for initial model-based clustering.
#' 
#' @return a \code{matrix} with 9 columns:
#' \itemize{
#' \item{\code{untitled}}{ a \code{numeric}, the median value of the copy 
#'     number elements forming each segment}
#' \item{\code{untitled}}{ TODO}
#' \item{\code{mediandev}}{ a \code{numeric}, the median function of copy 
#' number relative to its central value for each segment}
#' \item{\code{segerr}}{ a \code{numeric}, the error estimate for the 
#' function of copy number for each segment}
#' \item{\code{centerz}}{ a \code{numeric} between \code{0} and \code{1}, the 
#' probability that the segment is in the central cluster}
#' \item{\code{cpb}}{ TODO}
#' \item{\code{maxz}}{ TODO}
#' \item{\code{maxzmean}}{ TODO}
#' \item{\code{maxzsigma}}{ TODO}
#' }
#'
#' @examples
#'
#' # TODO
#' 
#' ## Create a list that contain information about the segments and bins 
#' ## related to one 
#' 
#' 
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats median var
#' @importFrom mclust Mclust mclustBIC
#' @keywords internal
CNclusterNcenter <- function(segrat, blsize, minJoin, nTrial, bestBIC,
    modelNames, cweight, bstimes, chromRange, keepClust=FALSE) {


    startcol <- "StartProbe"
    endcol   <- "EndProbe"
    chromcol <- "chrom"
    medcol   <- "segmedian"
    madcol   <- "segmad"

    ## Add to the current data.frame a column with the median and a column
    ## with the median absolute deviation for each group of bins forming 
    ## a specific segment
    # Modified for weigh
    segrat$seg <- cbind(segrat$seg,
                        t(apply(segrat$seg[,c(startcol, endcol, chromcol),
                        drop=FALSE], 1, smedmad, v=segrat$rat,
                        w=segrat$weight)))
    
    ## Assigning the column names of the updated data.frame
    dimnames(segrat$seg)[[2]] <- c(startcol, endcol, chromcol, medcol, madcol)
    
    ## Only retain the segments for the selected chromosomes
    seguse <- segrat$seg[segrat$seg[, chromcol] %in% chromRange,, drop=FALSE]
    
    ## Identify bins that are associated to the retained segments
    ## Only those bins will be kept for the analysis
    # aux <- rep(0, length(segrat$rat))
    # aux[seguse[, startcol]] <- 1  ## start position == 1
    # aux[seguse[, endcol]]   <- (-1) ## end position == -1
    # aux <- cumsum(aux) ## bins to retained are tagged 1 (except end position)
    # aux[seguse[, endcol]]   <- 1 ## end position == 1
    # 
    # ratuse <- segrat$rat[aux == 1]
    # 
    # # Modified for weight
    # if (length(segrat$weight) > 0) {
    #     weightuse <- segrat$weight[aux == 1]
    # } else{
    #     weightuse <- NULL
    # }
    
    ## Run trials and keep best cluster from all trial (best BIC value)
    for(j in seq_len(nTrial)) {
        # Modified for weight
        aaa <- segsample(seguse, segrat$rat, blocksize=blsize, weightcol=segrat$weight)
        if (all(unique(aaa[,3]) == 0)) { 
            aaa[,3] <- 1e-10 
        }
        emfit <- Mclust(aaa[,3], maxG=15, modelNames=modelNames)
        if (emfit$bic >= bestBIC) {
            bestaaa <- aaa
            bestem  <- emfit
            bestBIC <- emfit$bic
        }
    }
    
    ## Join clusters with minimum overlap
    ## The returned object is no more a Mclust object
    newem <- consolidate(bestem, minJoin)
    
    ## Join clusters until the main cluster contain the minimum required 
    ## ratio of data
    newem <- get.center(newem, cweight)
    clustRes <- NULL
    if(keepClust){
        clustRes <- list(bestaaa=bestaaa,
                         bestem=bestem,
                         newem=newem)
    } 
    ## Get the median of the central cluster. The central cluster is
    ## the one with the mean closer to zero
    if (length(bestem$parameters$mean) == 1) { 
        profcenter <- median(bestaaa[,3]) 
    } else { 
        profcenter <- weighted.median(bestaaa[,3], newem$z[, newem$center]) 
    }

    ## Calculate the median deviation to the central cluster
    mediandev <- segrat$seg[, medcol] - profcenter
    
    ## Bin sampling
    # Modified for weight
    segs <- segsample(segrat$seg, segrat$rat, times=bstimes, 
                      weightcol=segrat$weight)

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
    res <- list(seg=cbind(segrat$seg[, medcol], segrat$seg[,madcol], mediandev, segerr, 
                          centerz, cpb, maxz, maxzmean, maxzsigma),
                clustRes=clustRes)
    return(res)
}
