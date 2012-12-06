\name{segexample}
\alias{segexample}
\docType{data}
\title{
Example of a segmented copy number table.
}
\description{
Segmented log ratio data for 5 breast cancer genomes, derived using Representational
Oligonucleotide Microarray Analysis (ROMA) platform. ROMA detects genomic 
amplifications and deletions with boundaries defined at a resolution of 50 kb.
In this segmented table, each row represnts a segment.
}
\usage{data(segexample)}
\format{
  A data frame with 479 rows/segments and 12 columns/variables.
  \describe{
    \item{\code{ID}}{a character vector of profile IDs}
    \item{\code{start}}{a numeric vector (segment start probe number)}
    \item{\code{end}}{a numeric vector (segment end probe number)}
    \item{\code{num.probes}}{a numeric vector (number of probes in the segment)}
    \item{\code{seg.median}}{a numeric vector (median log ratio)}
    \item{\code{chrom}}{a numeric vector (chromosome number)}
    \item{\code{chrom.pos.start}}{a numeric vector (genomic start)}
    \item{\code{chrom.pos.end}}{a numeric vector (genomic end)}
    \item{\code{cytoband.start}}{a character vector (cytogenetic band start)}
    \item{\code{cytoband.end}}{a character vector (cytogenetic band end)}
    \item{\code{abs.pos.start}}{a numeric vector (genomic start, absolute)}
    \item{\code{abs.pos.end}}{a numeric vector (genomic end, absolute)}
  }
}
\details{
Segment medians are computed from log copy number ratio. The corresponding raw data
table is \code{ratexample} in this package.
}
\source{
Genome Res. 2006 Dec;16(12):1465-79.

Novel patterns of genome rearrangement and their association with survival in
breast cancer.

Hicks J, Krasnitz A, Lakshmi B, Navin NE, Riggs M, Leibu E, Esposito D, 
Alexander J, Troge J, Grubor V, Yoon S, Wigler M, Ye K, Borresen-Dale AL, 
Naume B, Schlicting E, Norton L, Hagerstrom T, Skoog L, Auer G, Maner S,
Lundin P, Zetterberg A.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
data(segexample)
}
\keyword{datasets}
