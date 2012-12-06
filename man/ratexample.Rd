\name{ratexample}
\alias{ratexample}
\docType{data}
\title{
Example of copy number log ratio data
}
\description{
Log ratio data for 5 breast cancer genomes, derived using Representational 
Oligonucleotide Microarray Analysis (ROMA), an array-based hybridization 
method that uses genomic complexity reduction based on representations.
}
\usage{data(ratexample)}
\format{
A log ratio data matrix with rows of 83055 oligonucleotide probes, and columns of 5
breast tumors.
}
\details{
The values are natural log copy number ratios, consistent with data in
\code{segexample} (segmented data for these tumors) and \code{normsegs}.
These copy number ratios are normalized using an intensity-based lowess 
curve fitting algorithm.
}
\source{
Genome Res. 2006 Dec;16(12):1465-79.

Novel patterns of genome rearrangement and their association with survival in
breast cancer.

Hicks J, Krasnitz A, Lakshmi B, Navin NE, Riggs M, Leibu E, Esposito D, 
Alexander J, Troge J, Grubor V, Yoon S, Wigler M, Ye K, Borresen-Dale AL, 
Naume B, Schlicting E, Norton L, Hagerstrom T, Skoog L, Auer G, Maner S,
Lundin P, Zetterberg A.
}
\examples{
data(ratexample)
#Plot the whole genome log ratio data for the first profile 
#Note X and Y chromosomes at the far  right of the plot
plot(ratexample[,1])
}
\keyword{datasets}
