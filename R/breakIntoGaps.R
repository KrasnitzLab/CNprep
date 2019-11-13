#' @title Redistribute probes/bins in removed segments into the adjacent 
#' segments.
#' 
#' @description This function redistributes the probes/bins from a removed 
#' segment into the adjacent segments for one specific chromosome at the time.
#' 
#' @param segtable a \code{matrix} or a \code{data.frame} with 3 columns 
#' named or enumerated by the values of 
#' \code{StartProbe, EndProbe, gapind}. 
#' 
#' @param gapind a \code{character} string specifying the name of the column
#' in \code{segtable} that tabulates the segments to remove. The value \code{1}
#' identifies segments to remove while \code{0} identifies segments to keep.
#' 
#' @param StartProbe a \code{character} string specifying the name of the 
#' column in \code{segtable} that tabulates the (integer) start postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param EndProbe a \code{character} string specifying the name of the
#' column in \code{segtable} that tabulates the (integer) end postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @return a \code{matrix} of \code{numeric} with 2 columns: 
#' \enumerate{
#'     \item the updated start position (integer) of each segment in internal 
#'     units. This column has the same name than the \code{StartProbe} 
#'     parameter.
#'     \item the updated end position (integer) of each segment in internal 
#'     units. This column has the same name than the \code{EndProbe} parameter.
#' }
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats runif
#' @keywords internal
breakIntoGaps <- function(segtable, gapind, StartProbe, EndProbe) {
    
    ## There is no probe/bin to remove, the positions can be returned without 
    ## modification
    if (sum(segtable[, gapind]) == 0) {
        return(as.matrix(segtable[, c(StartProbe, EndProbe)]))
    }
    
    gapstep  <- segtable[, gapind] - c(0, segtable[-nrow(segtable), gapind])
    gapstart <- which(gapstep == 1)
    gapend   <- which(gapstep == -1) - 1
    
    if(length(gapend) < length(gapstart)) {
        gapend <- c(gapend, nrow(segtable))
    }
    
    ranfrac <- runif(n = length(gapend))
    ranfrac[gapstart == 1] <- 1
    ranfrac[gapend == nrow(segtable)] <- 0
    midpoint <- round(ranfrac*segtable[gapstart, StartProbe] +
                (1 - ranfrac) * segtable[gapend, EndProbe])
    segtable[(gapend + 1)[gapend != nrow(segtable)], StartProbe] <-
                midpoint[gapend != nrow(segtable)]
    segtable[(gapstart - 1)[gapstart != 1], EndProbe] <-
                ifelse(gapstart[gapstart != 1] != nrow(segtable),
                midpoint[gapstart != 1] - 1, midpoint[gapstart != 1])
    
    return(as.matrix(segtable[, c(StartProbe, EndProbe)]))
}
