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
#' numbers and endpoint positions of the interval events. 
#' 
#' @param chromcol A \code{character} string or \code{numeric} representing an
#' integer  
#' specifying the column of \code{imat} containing the chromosome numbers 
#' of the interval events.
#' 
#' @param startcol A \code{character} string or \code{numeric} representing an
#' integer   
#' specifying the column of \code{imat} containing the left (right) endpoint of
#' the interval events.
#' 
#' @param endcol A \code{character} string or \code{numeric} representing an
#' integer 
#' specifying the column of \code{imat} containing the right endpoint 
#' of the interval events.
#' 
#' @param nprof An \code{numeric} acting as an integer specifying the number 
#' of copy number profiles from which the events originated. Default: \code{1}.
#' 
#' @param uthresh A \code{numeric} specifying the upper threshold for 
#' the event frequency or (if \code{nprof = 1}) for the event count.
#' 
#' @param dthresh A \code{numeric} specifying the upper and lower thresholds for 
#' the event frequency or (if \code{nprof = 1}) for the event count.
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
#' ampCNPmask <- makeCNPmask(imat=amps, chromcol="chrom",
#'     startcol="chrom.start", endcol="chrom.end", nprof=1203,
#'     uthresh=0.02, dthresh=0.008)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @export
makeCNPmask <- function(imat, chromcol=1, startcol=2, endcol=3, nprof=1,
                        uthresh, dthresh)
{
    ## Call makeCNPmask.chrom for each chromosome separately
    CNPmask <- by(imat, INDICES = as.factor(imat[, chromcol]), 
                    FUN=makeCNPmask.chrom,
                    startcol=startcol,
                    endcol=endcol,
                    nprof=nprof,
                    uthresh=uthresh,
                    dthresh=dthresh, simplify=TRUE)
    
    ## Create a matrix containing all results
    myCNPmask <- matrix(ncol=2, byrow=TRUE, 
                        data=unlist(lapply(CNPmask, t)))
    
    myCNPmask <- cbind(unlist(lapply(seq_len(length(unique(imat[, chromcol]))),
                    FUN=function(x) rep(as.numeric(names(CNPmask)[x]), 
                                nrow(CNPmask[[x]])))), myCNPmask)
    
    dimnames(myCNPmask)[[2]] <- c("chrom", "start", "end")
    
    return(myCNPmask)
}
