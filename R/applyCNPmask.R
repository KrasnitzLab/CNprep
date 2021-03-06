#' @title Apply a mask to a table of copy number events.
#' 
#' @description A mask is applied to amplified or deleted segments as 
#' tabulated in \code{segTable}. A decision whether to mask a segment 
#' is taken based on what portion of the segment is covered by the mask. A 
#' position is chosen at random within a segment to be masked, the flanking 
#' segments are extended to that position and the segment to be masked is 
#' indicated as such in the value returned.
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
#' @param startProbe a \code{character} string specifying the names of 
#' columns in \code{segTable} that tabulates the (integer) start postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param endProbe a \code{character} string specifying the names of 
#' columns in \code{segTable} that tabulates the (integer) end postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param eventIndex a \code{character} string giving the name of a column in 
#' \code{segTable} where copy number variation status of the segments is 
#' tabulated. 
#' 
#' @param maskTable a \code{matrix} or a \code{data.frame} with columns named 
#' or enumerated as given by \code{maskChrom, maskStart, maskEnd, maskIndex} 
#' and with rows corresponding to genomic intervals that comprise the mask.
#' 
#' @param maskChrom a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{maskTable} that tabulates 
#' the chromosome number of the intervals comprising the mask. 
#' 
#' @param maskStart a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{maskTable} that tabulates 
#' the genomic start coordinates of the intervals comprising the mask. 
#' 
#' @param maskEnd a \code{character} string or \code{integer} 
#' specifying the name or number of columns in \code{maskTable} that tabulates 
#' the genomic end coordinates of the intervals comprising the mask. 
#' 
#' @param maskIndex a \code{numeric} \code{vector} corresponding to 
#' \code{eventIndex}, specifying copy number events status for measuring units.
#' 
#' @param minCover a \code{numeric} value specifying the minimal portion of 
#' the segment that must be covered by the mask in order to trigger masking.
#' Default: \code{1}.
#' 
#' @param indexVals a \code{numeric} \code{vector} of length 2 specifying 
#' the two values in \code{maskIndex} to be matched with values in 
#' \code{eventIndex} to determine the events that are to be masked.
#' Default: \code{c(-1, 1)}.
#' 
#' @return a \code{matrix} with same number of observations/rows as 
#' \code{segTable} and with following three columns:
#' \itemize{
#' \item{StartProbe}{ an \code{numeric}, used as \code{integer},
#' for the start position of the segments after masking. }
#' \item{EndProbe}{ an \code{numeric}, used as integer
#' for the end position of the segments after masking. }
#' \item{toremove}{ an \code{numeric} \code{vector} used as integer 
#' whose values are 1 if the segment is masked and 0 otherwise. }
#' }
#' 
#' @details Masking is performed separately for each value in 
#' \code{indexVals}. Segments (rows of \code{segTable}) with that 
#' value of \code{eventIndex} are examined for coverage by mask intervals 
#' with that same value of \code{maskIndex} in \code{maskTable}. If the 
#' coverage is at least \code{minCover}, the segment is slated for masking, 
#' while its flanking segments are extended to a random point within the 
#' segment being masked.
#' 
#' @examples
#' 
#' ## Load datasets
#' data(segexample)
#' data(ratexample)
#' data(normsegs)
#' data(cnpexample)
#' 
#' ## Create a table with segment information (table of copy number events)
#' segtable <- CNpreprocessing(segall = segexample[segexample[,"ID"] == "WZ1",],
#'     ratall = ratexample, idCol = "ID", startCol = "start", endCol = "end",
#'     chromCol = "chrom", bpStartCol = "chrom.pos.start", 
#'     bpEndCol = "chrom.pos.end", blsize = 50, minJoin = 0.25, cWeight = 0.4,
#'     bsTimes = 50, chromRange = 1:22, 
#'     modelNames = "E", normalLength = normsegs[,1], 
#'     normalMedian = normsegs[,2])
#' 
#' ## Add an eventIndex column to segtable that identifies the 
#' ## amplication (marked as 1) and deletion (marked as -1) events
#' eventIndex <- rep(0, nrow(segtable))
#' eventIndex[segtable[,"marginalprob"] < 1e-4 & segtable[,"negtail"] > 0.999 & 
#'     segtable[,"mediandev"] < 0] <- -1
#' eventIndex[segtable[,"marginalprob"] < 1e-4 & segtable[,"negtail"] > 0.999 &
#'     segtable[,"mediandev"] > 0] <- 1
#' segtable <- cbind(segtable, eventIndex)
#' 
#' ## Create a mask table using amplification and deletion regions as input
#' namps17 <- cnpexample[cnpexample[,"copy.num"] == "amp",]
#' aCNPmask <- makeCNPmask(imat=namps17, chromCol=2, startCol=3, 
#'     endCol=4, nProf=1203, uThresh=0.02, dThresh=0.008)
#' ndels17 <- cnpexample[cnpexample[,"copy.num"] == "del",]
#' dCNPmask <- makeCNPmask(imat=ndels17, chromCol=2, startCol=3, 
#'     endCol=4, nProf=1203, uThresh=0.02, dThresh=0.008)
#' maskTable <- rbind(cbind(aCNPmask, cnpindex=1), 
#'     cbind(dCNPmask, cnpindex=-1))
#' 
#' ## Apply a mask to a table of copy number events
#' myCNPtable <- applyCNPmask(segTable=segtable, chrom="chrom",
#'     startPos="chrom.pos.start", endPos="chrom.pos.end", 
#'     startProbe="start", endProbe="end", eventIndex="eventIndex",
#'     maskTable=maskTable, maskChrom="chrom", maskStart="start", 
#'     maskEnd="end", maskIndex="cnpindex", minCover=0.005,
#'     indexVals=c(-1, 1))
#' 
#' ## Show some results
#' tail(myCNPtable)
#' 
#' 
#' @author Alexander Krasnitz
#' @export
applyCNPmask <- function(segTable, chrom, startPos, endPos, startProbe,
                    endProbe, eventIndex, maskTable, maskChrom, maskStart,
                    maskEnd, maskIndex, minCover=1, indexVals=c(-1, 1)) 
{
    breakCNPs <- by(segTable, INDICES=as.factor(segTable[,chrom]),
        FUN=breakIntoCNPs.chrom, chrom=chrom, startPos=startPos, 
        endPos=endPos, startProbe=startProbe, endProbe=endProbe, 
        eventIndex=eventIndex, cnpTable=maskTable, cnpChrom=maskChrom, 
        cnpStart=maskStart, cnpEnd=maskEnd, cnpIndex=maskIndex, 
        minCover=minCover, indexVals=indexVals, simplify=TRUE)

    myCNPs <- matrix(ncol=3, byrow=TRUE, 
                        data=unlist(lapply(breakCNPs, t)))
    
    dimnames(myCNPs)[[2]] <- c("StartProbe", "EndProbe", "toremove")
    
    return(as.matrix(myCNPs))
}
