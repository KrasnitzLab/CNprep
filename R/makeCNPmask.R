#' @title Given a set of copy-number events, create a DNA copy-number mask.
#' 
#' @description The function takes as an input a set of intervals with 
#' integer-valued boundary positions. It then finds interval regions where 
#' the event count is above each of two thresholds, upper and lower, and 
#' returns those interval regions with the count above the lower 
#' threshold that contain interval regions with the count above the upper 
#' threshold.
#' 
#' @param imat A \code{matrix} or a \code{data.frame} tabulating the chromosome 
#' numbers and startpoint/endpoint positions of the interval events. 
#' 
#' @param chromCol A \code{character} string or \code{numeric} representing an
#' integer specifying the column of \code{imat} containing the chromosome 
#' numbers of the interval events.
#' 
#' @param startCol A \code{character} string or \code{numeric} representing an
#' integer   
#' specifying the column of \code{imat} containing the left (right) endpoint of
#' the interval events.
#' 
#' @param endCol A \code{character} string or \code{numeric} representing an
#' integer specifying the column of \code{imat} containing the right endpoint 
#' of the interval events.
#' 
#' @param nProf A positive \code{numeric} acting as an integer specifying 
#' the number of copy number profiles from which the events originated. 
#' Default: \code{1}.
#' 
#' @param uThresh A \code{numeric} specifying the upper threshold for 
#' the event frequency or (if \code{nProf = 1}) for the event count.
#' 
#' @param dThresh A \code{numeric} specifying the upper and lower thresholds 
#' for the event frequency or (if \code{nProf = 1}) for the event count.
#' 
#' @return A \code{matrix} of \code{numeric} (used as integer) 
#' with three columns, called 
#' "chrom","start" and "end", specifying the chromosome number and 
#' boundary positions of the mask.
#' 
#' @details Masking is performed separately for each value in 
#' \code{indexvals}. Segments (rows of \code{segtable}) with that 
#' value of \code{eventIndex} are examined for coverage by mask intervals 
#' with that value of \code{maskindex} in \code{masktable}. If the coverage 
#' is at least \code{mincover}, the segment is slated for masking, while its 
#' flanking segments are extended to a random point within the segment 
#' being masked.
#' 
#' @examples
#' 
#' ## Load a table of copy number events collected from 1203 profiles.
#' data(cnpexample)
#' 
#' ## Create a table of gain (amplification) events only.
#' amps <- cnpexample[cnpexample[,"copy.num"] == "amp",]
#' 
#' # Create a mask using this table.
#' ampCNPmask <- makeCNPmask(imat=amps, chromCol="chrom",
#'     startCol="chrom.start", endCol="chrom.end", nProf=1203,
#'     uThresh=0.02, dThresh=0.008)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @export
makeCNPmask <- function(imat, chromCol=1, startCol=2, endCol=3, nProf=1,
                        uThresh, dThresh)
{
    ## Call makeCNPmask.chrom for each chromosome separately
    CNPmask <- by(imat, INDICES=as.factor(imat[, chromCol]), 
                    FUN=makeCNPmask.chrom,
                    startcol=startCol,
                    endcol=endCol,
                    nprof=nProf,
                    uthresh=uThresh,
                    dthresh=dThresh, simplify=TRUE)
    
    ## Create a matrix containing all results
    myCNPmask <- matrix(ncol=2, byrow=TRUE, 
                        data=unlist(lapply(CNPmask, t)))
    
    myCNPmask <- cbind(unlist(lapply(seq_len(length(unique(imat[, chromCol]))),
                    FUN=function(x) rep(as.numeric(names(CNPmask)[x]), 
                                nrow(CNPmask[[x]])))), myCNPmask)
    
    dimnames(myCNPmask)[[2]] <- c("chrom", "start", "end")
    
    return(myCNPmask)
}
