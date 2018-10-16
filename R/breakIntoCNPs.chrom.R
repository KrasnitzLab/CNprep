#' @title TODO
#' 
#' @description TODO
#' 
#' @param segtable a \code{matrix} or a \code{data.frame} with columns 
#' named or enumerated by the values of 
#' \code{chrom, startPos, endPos, startProbe, endProbe, eventIndex}. 
#' 
#' @param chrom a \code{character} string specifying the name for the column in 
#' \code{segtable} tabulating the (integer) chromosome number for each segment.
#' 
#' @param startPos a \code{character} string or integer specifying the 
#' name or number of columns in \code{segtable} that tabulates the (integer) 
#' genomic start coordinate of each segment.
#' 
#' @param endPos a \code{character} string or integer specifying the 
#' name or number of columns in \code{segtable} that tabulates the (integer) 
#' genomic end coordinate of each segment.
#' 
#' @param startProbe a \code{character} string specifying the names of 
#' columns in \code{segtable} that tabulates the (integer) start postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param endProbe a \code{character} string specifying the names of 
#' columns in \code{segtable} that tabulates the (integer) end postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param eventIndex a \code{character} string giving the name of a column in 
#' \code{segtable} where copy number variation status of the segments is 
#' tabulated. 
#' 
#' @param cnptable a \code{matrix} or a \code{data.frame} with columns named or 
#' enumerated as given by \code{cnpchrom, cnpstart, cnpend, cnpindex} and 
#' with rows corresponding to genomic intervals that comprise the mask.
#' 
#' @param cnpchrom a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{cnptable} that tabulates 
#' the chromosome number of the intervals comprising the mask. 
#' 
#' @param cnpstart a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{masktable} that tabulates 
#' the genomic start coordinates of the intervals comprising the mask.
#' 
#' @param cnpend a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{masktable} that tabulates 
#' the genomic end coordinates of the intervals comprising the mask. 
#' 
#' @param cnpindex a \code{numeric} \code{vector} corresponding to 
#' \code{eventIndex}, specifying copy number events status for measuring units.
#' 
#' @param mincover A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the degree of overlap above which two clusters will be joined 
#' into one.
#' 
#' @param indexvals a \code{numeric} \code{vector} of length 2 
#' specifying the two values in \code{maskindex} to be matched 
#' with values in \code{eventIndex} to 
#' determine the events that are to be masked.
#' 
#' @return TODO
#'
#' @examples
#'
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom mclust Mclust
#' @keywords internal
breakIntoCNPs.chrom <- function(segtable, chrom, startPos, endPos, startProbe,
    endProbe, eventIndex, cnptable, cnpchrom, cnpstart, cnpend, cnpindex, 
    mincover, indexvals)
{
    toremove <- rep(0,nrow(segtable))
    segtable <- cbind(segtable,toremove)
    chr <- segtable[1,chrom]
    if(sum(segtable[,eventIndex]!=0)==0|sum(cnptable[,cnpchrom]==chr)==0)
        return(as.matrix(segtable[,c(startProbe,endProbe,"toremove")]))
    cnpsinchr <- cnptable[cnptable[,cnpchrom]==chr,,drop=FALSE]
    for(i in indexvals)
        if(sum(segtable[,eventIndex]==i)>0&sum(cnpsinchr[,cnpindex]==i)>0) {
            acnpinchr <- cnpsinchr[cnpsinchr[,cnpindex]==i,,drop=FALSE]
            amps <- which(segtable[,eventIndex]==i)
            segstartmat <- matrix(ncol=nrow(acnpinchr),
                            data=rep(segtable[amps,startPos],nrow(acnpinchr)))
            segendmat <- matrix(ncol=nrow(acnpinchr),
                            data=rep(segtable[amps,endPos],nrow(acnpinchr)))
            cnpstartmat <- t(matrix(ncol=length(amps),
                            data=rep(acnpinchr[,cnpstart],length(amps))))
            cnpendmat <- t(matrix(ncol=length(amps),
                            data=rep(acnpinchr[,cnpend],length(amps))))
            cnpcover <- rowSums(pmax(matrix(nrow=nrow(cnpendmat),
                            ncol=ncol(cnpendmat),data=0),
                            (pmin(segendmat,cnpendmat)-
                            pmax(segstartmat,cnpstartmat)+1)))/
                            (segtable[amps,endPos]-segtable[amps,startPos]+1)
            toremove[amps[cnpcover>mincover]] <- 1
        }
    segtable[,"toremove"] <- toremove
    if(sum(toremove)>0) segtable[,c(startProbe,endProbe)] <-
        breakIntoGaps(segtable,"toremove",startProbe,endProbe)
    return(as.matrix(segtable[,c(startProbe,endProbe,"toremove")]))
}
