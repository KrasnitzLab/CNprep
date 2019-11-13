#' @title Redistribute bins in removed segments into the adjacent segments
#' 
#' @description This function redistributes the bins from a removed segment 
#' into the adjacent segments for one specific chromosome at the time.
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
#' # Table containing information about segments in the chromosome
#' segTable <- data.frame(ID = c(rep("WZ1", 5)), 
#'     start = c(1, 16, 23, 31, 38),
#'     end = c(15, 22, 30, 37, 50),
#'     seg.median = c(0.03779, -0.51546, 0.2431, -0.2259, 0.0372),
#'     chrom = c(rep(1, 5)),
#'     chrom.pos.start = c(932544, 16004440, 38093655, 78729960, 103416416),
#'     chrom.pos.end = c(15844870, 37974708, 78619856, 103394039, 142176090),
#'     eventIndex = c(0, 0, 1, 0, -1))
#'     
#' # Table with copy number information
#' cnptable <- matrix(c(rep(1, 3), c(932544, 38093688, 123416446), 
#'     c(11844870, 48619856, 182176090), rep(1,3)), ncol = 4, byrow = FALSE,
#'     dimnames=list(NULL, c("chrom", "start", "end", "cnpindex")))
#'     
#' # This function redistributes the bins from a removed segment into the 
#' # adjacent segments  
#' CNprep:::breakIntoCNPs.chrom(segtable = segTable, chrom = "chrom", 
#'     startPos = "chrom.pos.start", endPos = "chrom.pos.end", 
#'     startProbe = "start", endProbe = "end", 
#'     eventIndex = "eventIndex", 
#'     cnptable = cnptable, cnpchrom = "chrom",  
#'     cnpstart = "start", cnpend = "end", cnpindex = "cnpindex", 
#'     mincover = 0.005, indexvals = c(-1, 1))
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
breakIntoCNPs.chrom <- function(segtable, chrom, startPos, endPos, startProbe,
    endProbe, eventIndex, cnptable, cnpchrom, cnpstart, cnpend, cnpindex, 
    mincover, indexvals)
{
    toremove <- rep(0, nrow(segtable))
    segtable <- cbind(segtable, toremove)
    chr <- segtable[1, chrom]
    if (sum(segtable[,eventIndex] != 0) == 0 | 
                sum(cnptable[, cnpchrom] == chr) == 0) {
        return(as.matrix(segtable[,c(startProbe, endProbe, "toremove")]))
    }
    
    cnpsinchr <- cnptable[cnptable[, cnpchrom] == chr,, drop = FALSE]
    
    for (i in indexvals) {
        if (sum(segtable[,eventIndex] == i) > 0 & 
                    sum(cnpsinchr[, cnpindex] == i) > 0) {
            acnpinchr <- cnpsinchr[cnpsinchr[, cnpindex] == i,, drop = FALSE]
            amps <- which(segtable[, eventIndex] == i)
            segstartmat <- matrix(ncol = nrow(acnpinchr),
                                    data = rep(segtable[amps, startPos], 
                                    nrow(acnpinchr)))
            segendmat <- matrix(ncol = nrow(acnpinchr),
                            data = rep(segtable[amps, endPos], nrow(acnpinchr)))
            cnpstartmat <- t(matrix(ncol = length(amps),
                            data = rep(acnpinchr[, cnpstart], length(amps))))
            cnpendmat <- t(matrix(ncol = length(amps),
                            data = rep(acnpinchr[, cnpend], length(amps))))
            cnpcover <- rowSums(pmax(matrix(nrow = nrow(cnpendmat),
                            ncol = ncol(cnpendmat), data = 0), 
                            (pmin(segendmat, cnpendmat) -
                            pmax(segstartmat, cnpstartmat) + 1))) /
                            (segtable[amps, endPos] - 
                            segtable[amps, startPos] + 1)
            
            toremove[amps[cnpcover > mincover]] <- 1
        }
    }
    
    segtable[, "toremove"] <- toremove
    
    if (sum(toremove) > 0) {
        segtable[, c(startProbe, endProbe)] <-
            breakIntoGaps(segtable, "toremove", startProbe, endProbe)
    }
    
    return(as.matrix(segtable[, c(startProbe, endProbe, "toremove")]))
}
