% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/consolidate.R
\name{consolidate}
\alias{consolidate}
\title{Join clusters until the minimum degree of overlap is reached}
\usage{
consolidate(emfit, minover)
}
\arguments{
\item{emfit}{an object of class \code{Mclust} providing a 
mixture model estimation.}

\item{minover}{a single \code{numeric} value between \code{0} and \code{1} 
specifying the degree of overlap above which two clusters will be joined 
into one.}
}
\value{
a \code{list} containing information about the updated
clusters obtained from a mixture model estimation:
\itemize{
\item \code{mu} a \code{numeric} \code{vector} representing the mean 
    for each component. If there is more than one component, the 
    \emph{k}th element is the mean of the \emph{k}th component of the 
    mixture model.
\item \code{pro} a \code{vector} whose \emph{k}th component is the mixing 
    proportion for the \emph{k}th component of the mixture model. If 
    missing, equal proportions are assumed.
\item \code{z} a \code{numeric} \code{matrix} whose \emph{[i,k]}th entry 
    is the probability that observation \emph{i} in the test data belongs 
    to the \emph{k}th class.
\item \code{groups} a \code{matrix} of \code{double}, used as 
    \code{integer}, with the number of rows corresponding
    to the current number of clusters while the number of columns is 
    corresponding to the initial number of clusters. The presence of 
    \code{1} in position \emph{[i,k]} indicates that the initial 
    \emph{i}th cluster is now part of the new \emph{k}th cluster.
\item \code{ngroups} a \code{numeric}, used as an integer, giving the final
    number of clusters.
\item \code{sigmasq}  a \code{numeric} \code{vector} giving the common 
    variance for each component in the mixture model "E".
}
}
\description{
The function join clusters that have the minimum ratio of
overlap as specified by user.
}
\examples{

## Load Mclust object
data(EMexample)

## Group clusters that have at least 0.4\% of overlap
## The inital object has 5 clusters while the return object has only 
## 4 clusters
CNprep:::consolidate(EMexample, minover=0.004)

}
\author{
Alexander Krasnitz, Guoli Sun
}
\keyword{internal}
