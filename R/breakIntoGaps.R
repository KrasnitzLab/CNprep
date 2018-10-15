#' @title TODO
#' 
#' @description TODO
#' 
#' @param segtable a \code{matrix} or a \code{data.frame} TODO 
#' 
#' @param gapind a \code{character} string TODO
#' 
#' @param StartProbe a \code{character} string specifying the names of 
#' columns in \code{segtable} that tabulates the (integer) start postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @param EndProbe a \code{character} string specifying the names of 
#' columns in \code{segtable} that tabulates the (integer) end postion 
#' of each segment in internal units such as probe numbers for 
#' data of CGH microarray origin.
#' 
#' @return TODO
#' 
#' @examples
#' 
#' # TODO
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats runif
#' @keywords internal
breakIntoGaps <- function(segtable, gapind, StartProbe, EndProbe) {
    if(sum(segtable[,gapind])==0)
        return(as.matrix(segtable[,c(StartProbe,EndProbe)]))
    gapstep<-segtable[,gapind]-c(0,segtable[-nrow(segtable),gapind])
    gapstart<-which(gapstep==1)
    gapend<-which(gapstep==-1)-1
    if(length(gapend)<length(gapstart))gapend<-c(gapend,nrow(segtable))
    ranfrac<-runif(n=length(gapend))
    ranfrac[gapstart==1]<-1
    ranfrac[gapend==nrow(segtable)]<-0
    midpoint<-round(ranfrac*segtable[gapstart,StartProbe]+
        (1-ranfrac)*segtable[gapend,EndProbe])
    segtable[(gapend+1)[gapend!=nrow(segtable)],StartProbe]<-
        midpoint[gapend!=nrow(segtable)]
    segtable[(gapstart-1)[gapstart!=1],EndProbe]<-
        ifelse(gapstart[gapstart!=1]!=nrow(segtable),
        midpoint[gapstart!=1]-1,midpoint[gapstart!=1])
    return(as.matrix(segtable[,c(StartProbe,EndProbe)]))
}
