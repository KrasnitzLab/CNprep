\name{normsegs}
\alias{normsegs}
\docType{data}
\title{
A reference set of segments
}
\description{
A table of segment lengths and log copy number ratios for a large set of
human diploid genomes.
}
\usage{data(normsegs)}
\format{
 A data matrix with 43497 rows/segments and 2 columns/variables.
 \describe{
    \item{\code{length}}{a numeric vector of segment genomic length}
    \item{\code{segmedian}}{a numeric vector of segment median computed 
    from log copy number ratio}
}
}
\details{
The table originates in a set of copy number profiles of over a 1000 
individuals, obtained using  Representational Oligonucleotide Microarray 
Analysis (ROMA) technology. To ensure ploidy of 2 segments from X and Y 
chromosomes and segemnts shorter than 5Mb were excluded.
}
\source{
Science. 2007 Apr 20;316(5823):445-9. Epub 2007 Mar 15.

Strong association of de novo copy number mutations with autism.

Sebat J, Lakshmi B, Malhotra D, Troge J, Lese-Martin C, Walsh T, Yamrom B, 
Yoon S, Krasnitz A, Kendall J, Leotta A, Pai D, Zhang R, Lee YH, Hicks J, 
Spence SJ, Lee AT, Puura K, Lehtimaki T, Ledbetter D, Gregersen PK, Bregman J, 
Sutcliffe JS, Jobanputra V, Chung W, Warburton D, King MC, Skuse D, 
Geschwind DH, Gilliam TC, Ye K, Wigler M.
}
\examples{
data(normsegs)
}
\keyword{datasets}
