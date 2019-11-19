#' @title Determine which segments have to be removed and redistribute 
#' probes/bins from removed segments
#' 
#' @description This function determines which segments have to be removed and
#' redistributes the probes/bins from a removed segment 
#' into the adjacent segments for one specific chromosome at the time.
#' 
#' @param segTable a \code{matrix} or a \code{data.frame} with columns 
#' named or enumerated by the values of 
#' \code{chrom, startPos, endPos, startProbe, endProbe, eventIndex}. 
#' 
#' @param chrom a \code{character} string specifying the name for the column in 
#' \code{segTable} tabulating the (integer) chromosome number for each segment.
#' 
#' @param startPos a \code{character} string or integer specifying the 
#' name or number of columns in \code{segTable} that tabulates the (integer) 
#' genomic start coordinate of each segment.
#' 
#' @param endPos a \code{character} string or integer specifying the 
#' name or number of columns in \code{segTable} that tabulates the (integer) 
#' genomic end coordinate of each segment.
#' 
#' @param startProbe a \code{character} string specifying the name of 
#' column in \code{segTable} that tabulates the (integer) start position 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param endProbe a \code{character} string specifying the name of the
#' column in \code{segTable} that tabulates the (integer) end position 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param eventIndex a \code{character} string giving the name of the column in 
#' \code{segTable} where copy number variation status of the segments is 
#' tabulated. 
#' 
#' @param cnpTable a \code{matrix} or a \code{data.frame} with columns named or 
#' enumerated as given by \code{cnpChrom, cnpStart, cnpEnd, cnpIndex} and 
#' with rows corresponding to genomic intervals that comprise the mask.
#' 
#' @param cnpChrom a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{cnpTable} that tabulates 
#' the chromosome number of the intervals comprising the mask. 
#' 
#' @param cnpStart a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{masktable} that tabulates 
#' the genomic start coordinates of the intervals comprising the mask.
#' 
#' @param cnpEnd a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{masktable} that tabulates 
#' the genomic end coordinates of the intervals comprising the mask. 
#' 
#' @param cnpIndex a \code{numeric} \code{vector} corresponding to 
#' \code{eventIndex}, specifying copy number events status for measuring units.
#' 
#' @param minCover A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the degree of overlap above which two clusters will be joined 
#' into one.
#' 
#' @param indexVals a \code{numeric} \code{vector} of length 2 
#' specifying the two values in \code{maskindex} to be matched 
#' with values in \code{eventIndex} to 
#' determine the events that are to be masked.
#' 
#' @return a \code{matrix} of \code{numeric} with 3 columns: 
#' \enumerate{
#'     \item the updated start position (integer) of each segment in internal 
#'     units. This column has the same name than the \code{startProbe} 
#'     parameter.
#'     \item the updated end position (integer) of each segment in internal 
#'     units. This column has the same name than the \code{endProbe} parameter.
#'     \item the masked result: \code{0} when the segment is retained, 
#'     \code{1} when the segment
#'     is removed. This column is called "toremove".
#' }
#'
#' @examples
#'
#' # Table containing information about segments in the chromosome
#' segTable <- data.frame(ID=c(rep("WZ1", 5)), 
#'     start=c(1, 16, 23, 31, 38),
#'     end=c(15, 22, 30, 37, 50),
#'     seg.median=c(0.03779, -0.51546, 0.2431, -0.2259, 0.0372),
#'     chrom=c(rep(1, 5)),
#'     chrom.pos.start=c(932544, 16004440, 38093655, 78729960, 103416416),
#'     chrom.pos.end=c(15844870, 37974708, 78619856, 103394039, 142176090),
#'     eventIndex=c(0, 0, 1, 0, -1))
#'     
#' # Table with copy number information
#' cnptable <- matrix(c(rep(1, 3), c(932544, 38093688, 123416446), 
#'     c(11844870, 48619856, 182176090), rep(1,3)), ncol=4, byrow=FALSE,
#'     dimnames=list(NULL, c("chrom", "start", "end", "cnpindex")))
#'     
#' # This function redistributes the bins/probes from a removed segment into 
#' # the adjacent segments  
#' CNprep:::breakIntoCNPs.chrom(segTable=segTable, chrom="chrom", 
#'     startPos="chrom.pos.start", endPos="chrom.pos.end", 
#'     startProbe="start", endProbe="end", 
#'     eventIndex="eventIndex", 
#'     cnpTable=cnptable, cnpChrom="chrom",  
#'     cnpStart="start", cnpEnd="end", cnpIndex="cnpindex", 
#'     minCover=0.005, indexVals=c(-1, 1))
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
breakIntoCNPs.chrom <- function(segTable, chrom, startPos, endPos, startProbe,
    endProbe, eventIndex, cnpTable, cnpChrom, cnpStart, cnpEnd, cnpIndex, 
    minCover, indexVals)
{
    ## List of segments that will be masked
    ## Initialy, none are masked (all set to zero)
    toremove <- rep(0, nrow(segTable))
    segTable <- cbind(segTable, toremove)
    
    ## Get information about the currently analysed chromosome
    chr <- segTable[1, chrom]
    
    ## When there is not segment with an event or when the copy number
    ## table have not entry for the specific chromosome, 
    ## an empty matrix is returned
    if (sum(segTable[,eventIndex] != 0) == 0 | 
                sum(cnpTable[, cnpChrom] == chr) == 0) {
        return(as.matrix(segTable[,c(startProbe, endProbe, "toremove")]))
    }
    
    ## Only retain copy number in the specified chromosome
    cnpsinchr <- cnpTable[cnpTable[, cnpChrom] == chr,, drop=FALSE]
    
    ## Select segments that should be removed
    ## Each type of event (amplification/deletion) is treated separately 
    for (i in indexVals) {
        ## At least one event must be present in the segment table and in 
        ## the copy number table to run the analysis
        if (sum(segTable[,eventIndex] == i) > 0 & 
                    sum(cnpsinchr[, cnpIndex] == i) > 0) {
            
            ## Only retain copy numbers specific to the analysed event
            acnpinchr <- cnpsinchr[cnpsinchr[, cnpIndex] == i,, drop=FALSE]
            
            ## Create segment matrix containing only the analysed event
            amps <- which(segTable[, eventIndex] == i)
            segstartmat <- matrix(ncol = nrow(acnpinchr),
                                    data = rep(segTable[amps, startPos], 
                                    nrow(acnpinchr)))
            segendmat <- matrix(ncol = nrow(acnpinchr),
                            data = rep(segTable[amps, endPos], nrow(acnpinchr)))
            cnpstartmat <- t(matrix(ncol = length(amps),
                            data = rep(acnpinchr[, cnpStart], length(amps))))
            cnpendmat <- t(matrix(ncol = length(amps),
                            data = rep(acnpinchr[, cnpEnd], length(amps))))
            
            ## For each segment with a specific event, calculate the 
            ## ratio of coverage by all copy number event associated to 
            ## the same event
            cnpcover <- rowSums(pmax(matrix(nrow = nrow(cnpendmat),
                            ncol = ncol(cnpendmat), data = 0), 
                            (pmin(segendmat, cnpendmat) -
                            pmax(segstartmat, cnpstartmat) + 1))) /
                            (segTable[amps, endPos] - 
                            segTable[amps, startPos] + 1)
            
            ## The segments with a coverage superior to the minimum coverage
            ## specified are marked to remove
            toremove[amps[cnpcover > minCover]] <- 1
        }
    }
    
    segTable[, "toremove"] <- toremove
    
    ## When at least one segment is marked to be removed, 
    ## the probes/bins in the removed segments are assigned to adjacent segments
    if (sum(toremove) > 0) {
        segTable[, c(startProbe, endProbe)] <-
            breakIntoGaps(segTable, "toremove", startProbe, endProbe)
    }
    
    return(as.matrix(segTable[, c(startProbe, endProbe, "toremove")]))
}
