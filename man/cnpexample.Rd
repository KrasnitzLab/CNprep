\name{cnpexample}
\alias{cnpexample}
\docType{data}
\title{
Example of a boundary positions table.
}
\description{
A table of genomic positions for DNA copy-number changing events, collected 
from genomes of 1203 individuals 
using Representational Oligonucleotide Microarray Analysis (ROMA) 
platform.
}
\usage{data(cnpexample)}
\format{
  A data frame with 19188 rows and 4 columns.
  \describe{
    \item{\code{copy.num}}{a character vector indicating whether an event is
a gain ("amp") or a loss ("del"). }
    \item{\code{chrom}}{a numeric vector indicating which chromosome the 
event is in.}
    \item{\code{chrom.start}}{a numeric vector of event start positions.}
    \item{\code{chrom.end}}{a numeric vector of event start positions.}
  }
}
\source{
Strong association of de novo copy number mutations with autism.
Sebat J, Lakshmi B, Malhotra D, Troge J, Lese-Martin C, Walsh T, Yamrom B, 
Yoon S, Krasnitz A, Kendall J, Leotta A, Pai D, Zhang R, Lee YH, Hicks J, 
Spence SJ, Lee AT, Puura K, Lehtimaki T, Ledbetter D, Gregersen PK, Bregman J, 
Sutcliffe JS, Jobanputra V, Chung W, Warburton D, King MC, Skuse D, 
Geschwind DH, Gilliam TC, Ye K, Wigler M.
Science. 2007 Apr 20;316(5823):445-9.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
data(cnpexample)
}
\keyword{datasets}
