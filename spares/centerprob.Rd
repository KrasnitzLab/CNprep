\name{centerprob}
\alias{centerprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
centerprob(logr, emfit, zgroup, times, center)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logr}{
%%     ~~Describe \code{logr} here~~
}
  \item{emfit}{
%%     ~~Describe \code{emfit} here~~
}
  \item{zgroup}{
%%     ~~Describe \code{zgroup} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (logr, emfit, zgroup, times, center) 
{
    gz <- matrix(nrow = length(logr), ncol = length(emfit$parameters$mean), 
        data = logr)
    if (length(emfit$parameters$mean) == 1) 
        epro <- as.vector(1)
    if (length(emfit$parameters$mean) > 1) 
        epro <- emfit$parameters$pro
    gz <- t(epro * pnorm(-abs(t(gz) - emfit$parameters$mean)/sqrt(emfit$parameters$variance$sigmasq)))
    gz <- gz \%*\% t(zgroup)
    gz <- matrix(ncol = ncol(gz), data = apply(gz, 2, cumsum)[seq(from = times, 
        to = nrow(gz), by = times), ]/times)
    gz <- (gz[, center] - c(0, gz[-nrow(gz), center]))/sum(epro[zgroup[center, 
        ] == 1])
    return(ifelse(gz < 0.5, gz, 1 - gz))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
