\name{makeCNPmask}
\alias{makeCNPmask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Given a set of copy-number events, create a DNA copy-number mask 
}
\description{
The function takes as an input a set of intervals with integer-valued boundary
positions. It then finds interval regions where the event count is above each 
of two thresholds, upper and lower, and returns those interval regions with
the count above the lower threshold that contain interval regions with the 
count above the upper threshold.
}
\usage{
makeCNPmask(imat,chromcol=1,startcol=2,endcol=3,nprof=1,uthresh,dthresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{imat}{
A matrix or a data frame tabulating the chromosome number and 
endpoint positions of the interval events. 
}
  \item{chromcol,startcol, endcol}{
Character strings or integers specifying the columns of 
\code{imat} containing the chromosome number and the left (right) endpoint 
positions of the interval events.
}
  \item{nprof}{
An integer specifying the number of copy number profiles from which the events
originate.
}
  \item{uthresh, dthresh}{
Numeric, specifying the upper and lower thresholds for the event frequency or
(if \code{nprof =1}) for the event count.
}
}
\value{
An integer matrix with three columns, called "chrom","start" and "end", specifying 
the chromosome number and boundary positions of the mask.
}
\author{
Alex Krasnitz, Guoli Sun
}
\examples{
#load a table of copy number events collected from 1203 profiles.
data(cnpexample)
#Create a table of gain (amplification) events only.
amps<-cnpexample[cnpexample[,"copy.num"]=="amp",]
#Create a mask using this table.
ampCNPmask<-makeCNPmask(imat=amps,chromcol="chrom",startcol="chrom.start",
endcol="chrom.end",nprof=1203,uthresh=0.02,dthresh=0.008)
}
