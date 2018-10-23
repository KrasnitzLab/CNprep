#' @title Update the start and end positions of each segment in function of 
#' the segments to be removed.
#' 
#' @description Dispatch the bins associated to a segment to be removed to
#' the surronding segments.
#' 
#' @param segtable a \code{matrix} or a \code{data.frame} for segmented copy 
#' number profiles that contains a column that identifies the segments to 
#' remove.  
#' 
#' @param gapind a \code{character} string specifying the name of the 
#' column in \code{segtable} that identify the segments to remove.
#' 
#' @param StartProbe a \code{character} string specifying the name of 
#' column in \code{segtable} that tabulates the (integer) start position 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param EndProbe a \code{character} string specifying the name of 
#' column in \code{segtable} that tabulates the (integer) end position 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @return a \code{matrix} containing 2 columns that represent the updated
#' start and end positions of each segment.
#' 
#' @examples
#' 
#' ## Load dataset
#' data(segexample)
#' 
#' ## Use only a subsection of the dataset
#' mysegs <- segexample[1:10, ]
#' 
#' ## Add colum identifying segment to remove
#' mysegs$toremove <- rep(c(0, 1), 5)
#' 
#' ## Dispatch the bins from the segments marked to remove (value of 1 in 
#' ## the toremove column) to the adjacent segments
#' CNprep:::breakIntoGaps(segtable = mysegs, gapind="toremove", 
#'     StartProbe="start", EndProbe="end")
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats runif
#' @keywords internal
breakIntoGaps <- function(segtable, gapind, StartProbe, EndProbe) {
    ## When no segment marked to remove, simply return the current
    ## start and end values
    if (sum(segtable[, gapind]) == 0) {
        return(as.matrix(segtable[, c(StartProbe,EndProbe)]))
    }
    
    gapstep <- segtable[, gapind] - c(0, segtable[-nrow(segtable), gapind])
    gapstart <- which(gapstep == 1)
    gapend <- which(gapstep == -1) - 1
    
    if (length(gapend) < length(gapstart)) {
        gapend <- c(gapend,nrow(segtable))
    }
    
    ranfrac <- runif(n = length(gapend))
    ranfrac[gapstart == 1] <- 1
    ranfrac[gapend == nrow(segtable)] <- 0
    
    midpoint <- round(ranfrac*segtable[gapstart, StartProbe] +
            (1-ranfrac) * segtable[gapend, EndProbe])
    
    segtable[(gapend + 1)[gapend != nrow(segtable)], StartProbe] <-
                        midpoint[gapend != nrow(segtable)]
    
    segtable[(gapstart - 1)[gapstart != 1], EndProbe] <-
                    ifelse(gapstart[gapstart != 1] != nrow(segtable),
                        midpoint[gapstart != 1] - 1, midpoint[gapstart != 1])
    
    return(as.matrix(segtable[, c(StartProbe, EndProbe)]))
}
