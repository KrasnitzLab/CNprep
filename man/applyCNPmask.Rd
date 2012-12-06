\name{applyCNPmask}
\alias{applyCNPmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply a mask to a table of copy number events.
}
\description{
A mask is applied to amplified or deleted segments as tabulated in 
\code{segtable}. A decision whether to mask a segment is taken based on what 
portion of the segment is covered by the mask. A position is chosen at random 
within a segment to be masked, the flanking segments are extended to that 
position and the segment to be masked is indicated as such in the value 
returned.
}
\usage{
applyCNPmask(segtable,chrom,startPos,endPos,startProbe,endProbe,
	eventIndex,masktable,maskchrom,maskstart,maskend,maskindex,mincover=1,
	indexvals=c(-1,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segtable}{
A matrix or a data frame with columns named or enumerated by the values of 
\code{chrom, startPos, endPos, startProbe, endProbe, eventIndex}. 
}
  \item{chrom}{
A character string specifying the name for the column in \code{segtable}
tabulating the (integer) chromosome number for each segment.
}
  \item{startPos, endPos}{
Character strings or integers specifying the names or numbers of columns in 
\code{segtable} that tabulate the (integer) genomic start and end coordinates 
of each segment.
}
  \item{startProbe, endProbe}{
Character strings specifying the names of columns in \code{segtable} that
tabulate the (integer) start and end postions of each segment in internal
units such as probe numbers for data of CGH microarray origin.
}
  \item{eventIndex}{
A character string giving the name of a column in \code{segtable} where copy number variation status of the segments is tabulated. 
}
  \item{masktable}{
A matrix or a data frame with columns named or enumerated as given by
\code{maskchrom, maskstart, maskend, maskindex} and with rows corresponding to
genomic intervals that comprise the mask.
}
  \item{maskchrom,maskstart,maskend}{
Character strings or integers specifying the names or numbers of columns 
in \code{masktable} that tabulate the chromosome number and (integer) genomic 
start and end coordinates of the intervals comprising the mask. 
}
  \item{maskindex}{
A numeric vector corresponding to \code{eventIndex},specifying copy number 
events status for measuring units.
}
  \item{mincover}{
A numeric value specifying the minimal portion of the segment that must be 
covered by the mask in order to trigger masking.
}
  \item{indexvals}{
A numeric vector of length 2 specifying the two values in \code{maskindex} to 
be matched with values in \code{eventIndex} to determine the events that are to 
be masked.
}
}
\details{
Masking is performed separately for each value in \code{indexvals}. Segments 
(rows of \code{segtable}) with that value of \code{eventIndex} are examined
for coverage by mask intervals with that value of \code{maskindex} in 
\code{masktable}. If the coverage is at least \code{mincover}, the segment is
slated for masking, while its flanking segments are extended to a random 
point within the segment being masked.
}
\value{
A matrix with same number of observations/rows as \code{segtable} and with 
following three columns:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{startProbe,endProbe}{
An integer vector for the start and end positions of the segments after masking.
}
\item{toremove}{
An integer vector whose values are 1 if the segment is masked and 0 otherwise.
}
}
\author{
Alex Krasnitz
}
\examples{
\dontrun{
data(segexample)
data(ratexample)
data(normsegs)
data(cnpexample)
segtable<-CNpreprocessing(segall=segexample[segexample[,"ID"]=="WZ1",],
ratall=ratexample,"ID","start","end",chromcol="chrom",bpstartcol="chrom.pos.start",
bpendcol="chrom.pos.end",blsize=50,minjoin=0.25,cweight=0.4,bstimes=50,
chromrange=1:22,distrib="Rparallel",njobs=2,modelNames="E",normalength=normsegs[,1],
normalmedian=normsegs[,2])
#form a eventIndex vector
eventIndex<-rep(0,nrow(segtable))
eventIndex[segtable[,"marginalprob"]<1e-4&segtable[,"negtail"]>
0.999&segtable[,"mediandev"]<0] <- -1
eventIndex[segtable[,"marginalprob"]<1e-4&segtable[,"negtail"]>
0.999&segtable[,"mediandev"]>0] <- 1
segtable<-cbind(segtable,eventIndex)
#form a cnpindex vector
namps17<-cnpexample[cnpexample[,"copy.num"]=="amp",]
aCNPmask<-makeCNPmask(imat=namps17,chromcol=2,startcol=3,endcol=4,
nprof=1203,uthresh=0.02,dthresh=0.008)
ndels17<-cnpexample[cnpexample[,"copy.num"]=="del",]
dCNPmask<-makeCNPmask(imat=ndels17,chromcol=2,startcol=3,endcol=4,
nprof=1203,uthresh=0.02,dthresh=0.008)
cnptable<-rbind(cbind(aCNPmask,cnpindex=1),cbind(dCNPmask,cnpindex=-1))
#run the CNP test
myCNPtable<-applyCNPmask(segtable,"chrom",startPos="chrom.pos.start",
endPos="chrom.pos.end","start","end","eventIndex",masktable=cnptable,"chrom",
maskstart="start",maskend="end",maskindex="cnpindex",mincover=0.005,indexvals=c(-1,1))
}
}
