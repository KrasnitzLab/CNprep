#' @title Given a set of copy-number events, create a DNA copy-number mask for
#' one specific chromosome.
#' 
#' @description The function takes as an input a set of intervals with 
#' integer-valued boundary positions. It then finds interval regions where 
#' the event count is above each of two thresholds, upper and lower, and 
#' returns those interval regions with the count above the lower 
#' threshold that contain interval regions with the count above the upper 
#' threshold. This is done for one specific chromosome.
#' 
#' @param imat A \code{matrix} or a \code{data.frame} tabulating the chromosome 
#' numbers and endpoint positions of the interval events.  
#' 
#' @param startcol A \code{character} string or \code{numeric} (representing an
#' integer) specifying the column of \code{imat} containing the left 
#' endpoint of the interval events. Default: \code{1}.
#' 
#' @param endcol A \code{character} string or \code{numeric} (representing an
#' integer) specifying the column of \code{imat} containing the right endpoint 
#' of the interval events. Default: \code{2}.
#' 
#' @param nprof A \code{numeric} acting as an integer specifying the number 
#' of copy number profiles from which the events originated. Default: \code{1}.
#' 
#' @param uthresh A \code{numeric} specifying the upper threshold for 
#' the event frequency or (if \code{nprof = 1}) for the event count.
#' 
#' @param dthresh A \code{numeric} specifying the upper and lower thresholds 
#' for the event frequency or (if \code{nprof = 1}) for the event count.
#' 
#' @return A \code{matrix} of \code{numeric} (used as integer) 
#' with three columns that represent the masked regions:
#' \itemize{
#' \item "chrom", the chromosome number
#' \item "start", the starting position of a masked region
#' \item  "end", the ending position of a masked region
#' }
#'
#' @examples
#'
#' ## Load a table of copy number events collected from 1203 profiles.
#' data(cnpexample)
#' 
#' ## Create a table with chromosome 2 only.
#' del <- cnpexample[cnpexample$chrom == 2,]
#' 
#' ## Create a table of deletion events only.
#' del <- del[del[,"copy.num"] == "del",]
#' 
#' # Create a mask for chromosome 2 using this table.
#' delCNPmask <- CNprep:::makeCNPmask.chrom(imat=del,
#'     startcol="chrom.start", endcol="chrom.end", nprof=1203,
#'     uthresh=0.02, dthresh=0.008)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
makeCNPmask.chrom <- function(imat, startcol=1, endcol=2, nprof=1, uthresh,
                                dthresh)
{
    astart <- imat[, startcol]
    aend <- imat[, endcol]
    z <- cbind(c(astart, aend, aend + 1),
        c(rep(1, length(astart)), rep(0,length(aend)), rep(-1, length(aend))))
    z <- z[order(z[, 1]), ]
    z[, 2] <- cumsum(z[, 2])
    z <- z[nrow(z) - rev(match(rev(unique(z[, 1])), rev(z[, 1]))) + 1, ]
    #z[,1] gives unique start and end positions; z[,2] gives event counts there
    
    ## Mark positions w/counts above upper thresh
    z <- cbind(z,z[,2] >= (uthresh*nprof)) 
    zsteps <- z[,3] - c(0, z[-nrow(z), 3])
    ustart <- z[zsteps == 1, 1]
    zsteps <- z[,3] - c(z[-1, 3], 0)
    
    ## Mark starts and ends of intervals w/count above upper thresh
    uend <- z[zsteps == 1, 1] 
    z[,3] <- z[,2] >= (dthresh*nprof)
    zsteps <- z[,3] - c(0, z[-nrow(z), 3])
    dstart <- z[zsteps == 1, 1]
    zsteps <- z[,3] - c(z[-1, 3], 0)
    dend <- z[zsteps == 1, 1] #likewise for the lower thresh
    if (length(ustart) > 0) {
        ci <- containment.indicator(ustart, uend, dstart, dend)
        return(matrix(ncol=2, data=c(dstart[ci[, 2] >= ci[, 1]],
                dend[ci[, 2] >= ci[, 1]]), 
                dimnames=list(NULL, c("start", "end"))))
    } #ie intervals above lower thresh with counts above upper thresh inside
    else {
        ## Return an empty matrix
        return(matrix(ncol=3, nrow=0, dimnames=list(NULL, 
                                            c("chrom", "start", "end"))))
    }
}
