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
#' # Table containing information about segments in the chromosome
#' # The table must have a column with the information about segments
#' # to be removed
#' segTable <- data.frame(ID = c(rep("WZ1", 5)), 
#'     start = c(1, 16, 23, 31, 38),
#'     end = c(15, 22, 30, 37, 50),
#'     seg.median = c(0.03779, -0.51546, 0.2431, -0.2259, 0.0372),
#'     chrom = c(rep(1, 5)),
#'     chrom.pos.start = c(932544, 16004440, 38093655, 78729960, 103416416),
#'     chrom.pos.end = c(15844870, 37974708, 78619856, 103394039, 142176090),
#'     eventIndex = c(0, 0, 1, 0, -1),
#'     toremove = c(0, 1, 0, 1, 0))
#'     
#' # This function redistributes the bins/probes from a removed segment into 
#' # the adjacent segments  
#' CNprep:::breakIntoGaps(segtable = segTable, gapind = "toremove",
#'     StartProbe = "start", EndProbe = "end")
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
    
    if (length(gapend) < length(gapstart)) {
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
