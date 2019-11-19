#' @title Redistribute probes/bins in removed segments into the adjacent 
#' segments.
#' 
#' @description This function redistributes the probes/bins from a removed 
#' segment into the adjacent segments for one specific chromosome at the time.
#' The function uses an uniform distribution to select the proportion of
#' probes/bins each adjacent segments will receive.
#' 
#' @param segTable a \code{matrix} or a \code{data.frame} with 3 columns 
#' named or enumerated by the values of 
#' \code{startProbe, endProbe, gapInd}. 
#' 
#' @param gapInd a \code{character} string specifying the name of the column
#' in \code{segTable} that tabulates the segments to remove. The value \code{1}
#' identifies segments to remove while \code{0} identifies segments to keep.
#' 
#' @param startProbe a \code{character} string specifying the name of the 
#' column in \code{segTable} that tabulates the (integer) start postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param endProbe a \code{character} string specifying the name of the
#' column in \code{segTable} that tabulates the (integer) end postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @return a \code{matrix} of \code{numeric} with 2 columns: 
#' \enumerate{
#'     \item the updated start position (integer) of each segment in internal 
#'     units. This column has the same name than the \code{startProbe} 
#'     parameter.
#'     \item the updated end position (integer) of each segment in internal 
#'     units. This column has the same name than the \code{endProbe} parameter.
#' }
#' 
#' @examples
#' 
#' ## Table containing information about segments in the chromosome
#' ## The table must have a column with the information about segments
#' ## to be removed (in the example, the column "toremove")
#' ## The table must also have start and end columns
#' segTable <- data.frame(ID=c(rep("WZ1", 5)), 
#'     start=c(1, 16, 23, 31, 38),
#'     end=c(15, 22, 30, 37, 50),
#'     seg.median=c(0.03779, -0.51546, 0.2431, -0.2259, 0.0372),
#'     chrom=c(rep(1, 5)),
#'     chrom.pos.start = c(932544, 16004440, 38093655, 78729960, 103416416),
#'     chrom.pos.end = c(15844870, 37974708, 78619856, 103394039, 142176090),
#'     eventIndex = c(0, 0, 1, 0, -1),
#'     toremove = c(0, 1, 0, 1, 0))
#'     
#' ## This function redistributes the bins/probes from a removed segment into 
#' ## the adjacent segments  
#' CNprep:::breakIntoGaps(segTable = segTable, gapInd = "toremove",
#'     startProbe = "start", endProbe = "end")
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats runif
#' @keywords internal
breakIntoGaps <- function(segTable, gapInd, startProbe, endProbe) {
    
    ## When there is no probe/bin to remove, the positions can be returned 
    ## without modification
    if (sum(segTable[, gapInd]) == 0) {
        return(as.matrix(segTable[, c(startProbe, endProbe)]))
    }
    
    gapstep  <- segTable[, gapInd] - c(0, segTable[-nrow(segTable), gapInd])
    gapstart <- which(gapstep == 1)
    gapend   <- which(gapstep == -1) - 1
    
    if (length(gapend) < length(gapstart)) {
        gapend <- c(gapend, nrow(segTable))
    }
    
    ranfrac <- runif(n = length(gapend))
    ranfrac[gapstart == 1] <- 1
    ranfrac[gapend == nrow(segTable)] <- 0
    midpoint <- round(ranfrac*segTable[gapstart, startProbe] +
                (1 - ranfrac) * segTable[gapend, endProbe])
    segTable[(gapend + 1)[gapend != nrow(segTable)], startProbe] <-
                midpoint[gapend != nrow(segTable)]
    segTable[(gapstart - 1)[gapstart != 1], endProbe] <-
                ifelse(gapstart[gapstart != 1] != nrow(segTable),
                midpoint[gapstart != 1] - 1, midpoint[gapstart != 1])
    
    return(as.matrix(segTable[, c(startProbe, endProbe)]))
}
