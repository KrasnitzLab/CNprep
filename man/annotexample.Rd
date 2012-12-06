\name{annotexample}
\alias{annotexample}
\docType{data}
\title{
Annotation table for ROMA CGH platform and human genome version 17.
}

\description{
Whole genome annotation table using Representational Oligonucleotide Microarray
Analysis (ROMA) CGH platform, human genome version 17.
}
\usage{data(annotexample)}
\format{
  A data frame with 83055 observations on the following 3 variables.
  \describe{
    \item{\code{PROBEID}}{a character vector}
    \item{\code{CHROM}}{a numeric vector}
    \item{\code{CHROM.POS}}{a numeric vector}
  }
}
\details{
The values in the chromosome column are all integer, with 23 corresponding to X,
24 to Y and 25 to a set of non-human test probes.
}
\source{
GEO accession GPL9775, http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL9775
}
\examples{
data(annotexample)
}
\keyword{datasets}
