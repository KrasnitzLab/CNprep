#' @title Apply a mask to a table of copy number events for a specified
#' chromosome.
#' 
#' @description A mask is applied to amplified or deleted segments specific
#' to one chromosome as 
#' tabulated in \code{segtable}. A decision whether to mask a segment 
#' is taken based on what portion of the segment is covered by the mask. A 
#' position is chosen at random within a segment to be masked, the flanking 
#' segments are extended to that position and the segment to be masked is 
#' indicated as such in the value returned.
#' 
#' @param segtable a \code{matrix} or a \code{data.frame} with columns 
#' named or enumerated by the values of 
#' \code{chrom, startPos, endPos, startProbe, endProbe, eventIndex}. The
#' \{matrix} must contain information for only one chromosome.
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
#' specifying the two values in \code{cnpindex} to be matched 
#' with values in \code{eventIndex} to 
#' determine the events that are to be masked.
#' 
#' @return a \code{matrix} with same number of observations/rows as 
#' \code{segtable} and with following three columns:
#' \itemize{
#' \item{startProbe,endProbe}{ An integer vector for the start and end 
#' positions of the segments after masking. }
#' \item{toremove}{ An integer vector whose values are 1 if the segment 
#' is masked and 0 otherwise. }
#' }
#' 
#' @details Masking is performed separately for each value in 
#' \code{indexvals}. Segments (rows of \code{segtable}) with that 
#' value of \code{eventIndex} are examined for coverage by mask intervals 
#' with that value of \code{cnpindex} in \code{cnptable}. If the coverage 
#' is at least \code{mincover}, the segment is slated for masking, while its 
#' flanking segments are extended to a random point within the segment 
#' being masked.
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
    ## Add a column named "toremove" and filled with zero to segtable
    toremove <- rep(0, nrow(segtable))
    segtable <- cbind(segtable, toremove)
    
    ## Get the current analysed chromosome
    chr <- segtable[1, chrom]
    
    ## AD: should be a long form ||
    ## When no segment are marked as containing a potential event (with
    ## non zero value in eventIndex column) or when no genomic interval related
    ## to mask are found for the chromosome, a matrix with unmodified results
    ## is returned
    if (sum(segtable[, eventIndex] != 0) == 0 | 
            sum(cnptable[, cnpchrom] == chr) == 0) {
        return(as.matrix(segtable[, c(startProbe, endProbe, "toremove")]))
    }
    
    ## Retain only genomic intervals related to current analysed chromosome
    cnpsinchr <- cnptable[cnptable[, cnpchrom] == chr, , drop=FALSE]
    
    ## Treat each event to be masked (as identified by an unique number) 
    ## separatly
    for (i in indexvals) {
        ## AD: should be a long form &&
        if (sum(segtable[, eventIndex] == i) > 0 & 
                sum(cnpsinchr[, cnpindex] == i) > 0) {
            acnpinchr <- cnpsinchr[cnpsinchr[, cnpindex] == i, , drop=FALSE]
            amps <- which(segtable[, eventIndex] == i)
            segstartmat <- matrix(ncol=nrow(acnpinchr),
                            data=rep(segtable[amps, startPos], nrow(acnpinchr)))
            segendmat <- matrix(ncol=nrow(acnpinchr),
                            data=rep(segtable[amps, endPos], nrow(acnpinchr)))
            cnpstartmat <- t(matrix(ncol=length(amps),
                            data=rep(acnpinchr[, cnpstart], length(amps))))
            cnpendmat <- t(matrix(ncol=length(amps),
                            data=rep(acnpinchr[, cnpend], length(amps))))
            cnpcover <- rowSums(pmax(matrix(nrow=nrow(cnpendmat),
                            ncol=ncol(cnpendmat), data=0),
                            (pmin(segendmat, cnpendmat) -
                            pmax(segstartmat,cnpstartmat) + 1)))/
                            (segtable[amps, endPos] - 
                            segtable[amps, startPos] + 1)
            toremove[amps[cnpcover > mincover]] <- 1
        }
    }
    
    segtable[, "toremove"] <- toremove
    
    if (sum(toremove) > 0) { 
        segtable[, c(startProbe, endProbe)] <- breakIntoGaps(segtable, 
                                                            "toremove",
                                                            startProbe,
                                                            endProbe)
    }
    
    return(as.matrix(segtable[, c(startProbe, endProbe, "toremove")]))
}
