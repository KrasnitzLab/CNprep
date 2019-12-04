#' @title Parameters validation for the \code{\link{makeCNPmask}} function
#' 
#' @description Validation of all parameters needed by the public
#' \code{\link{makeCNPmask}} function.
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
#' 
#' @param uThresh A \code{numeric} specifying the upper ratio threshold for 
#' the event frequency or (if \code{nProf = 1}) for the event count.
#' The upper ratio threshold must be equal or superior to the lower 
#' threshold (\code{dThresh}). The upper ratio threshold must also be between
#' 0 and 1.
#' 
#' @param dThresh A \code{numeric} specifying the lower ratio threshold 
#' for the event frequency or (if \code{nProf = 1}) for the event count.
#' The lower ratio threshold must be equal or inferior to the upper  
#' threshold (\code{dThresh}). 
#' The lower ratio threshold must also be between 0 and 1.
#' 
#' @return \code{0}. 
#' 
#' @examples
#'
#' data(segexample)
#' 
#' ## Return zero as all parameters are valid
#' CNprep:::makeCNPmask(imat=segexample, chromCol="chrom",
#'     startCol="start", endCol="end", nProf=1203, uThresh=0.20, dThresh=0.10)
#' 
#' @author Astrid Deschênes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateMakeCNPmask <- function(imat, chromCol, startCol, endCol, nProf,
                                uThresh, dThresh) 
{
    
    ## Validate that nProf is an positive integer
    if (!(isSingleInteger(nProf) || isSingleNumber(nProf)) ||
        as.integer(nProf) < 1) {
        stop("nProf must be a positive integer")
    }
    
    ## Validate that dThresh is a number
    if (!(isSingleNumber(dThresh))) {
        stop("dThresh must be a ratio between 0 and 1")
    }
    
    ## Validate that uThresh is a number
    if (!(isSingleNumber(uThresh))) {
        stop("uThresh must be a ratio between 0 and 1")
    }
    
    ## Validate that dThresh is betweeen 0 and 1
    if (!((dThresh >= 0) && (dThresh <= 1))) {
        stop("dThresh must be between 0 and 1")
    }
    
    ## Validate that uThresh is betweeen 0 and 1
    if (!((uThresh >= 0) && (uThresh <= 1))) {
        stop("uThresh must be between 0 and 1")
    }
    
    ## Validate that uThresh is superior or equal to dThresh
    if (!(uThresh >= dThresh)) {
        stop("uThresh must be equal or superior to dThresh")
    }
    
    return(0L)
}