\name{normalComparison}
\alias{normalComparison}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
normalComparison(normalmedian, normalength, tumormedian, tumorlength, normalmad = NULL, normalerror = NULL, tumormad = NULL, tumorerror = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{normalmedian}{
%%     ~~Describe \code{normalmedian} here~~
}
  \item{normalength}{
%%     ~~Describe \code{normalength} here~~
}
  \item{tumormedian}{
%%     ~~Describe \code{tumormedian} here~~
}
  \item{tumorlength}{
%%     ~~Describe \code{tumorlength} here~~
}
  \item{normalmad}{
%%     ~~Describe \code{normalmad} here~~
}
  \item{normalerror}{
%%     ~~Describe \code{normalerror} here~~
}
  \item{tumormad}{
%%     ~~Describe \code{tumormad} here~~
}
  \item{tumorerror}{
%%     ~~Describe \code{tumorerror} here~~
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
function (normalmedian, normalength, tumormedian, tumorlength, 
    normalmad = NULL, normalerror = NULL, tumormad = NULL, tumorerror = NULL) 
{
    nisnull <- c(!(is.null(normalmad) | is.null(tumormad)), !(is.null(normalerror) | 
        is.null(tumorerror)))
    nsred <- matrix(ncol = 2 + sum(nisnull), length(normalmedian), 
        data = c(normalength, normalmedian, if (nisnull[1]) normalmad, 
            if (nisnull[2]) normalerror), dimnames = list(NULL, 
            c("length", "mediandev", if (nisnull[1]) "segmad", 
                if (nisnull[2]) "segerr")))
    nsred <- nsred[order(nsred[, "mediandev"]), , drop = F]
    lnorm <- sum(nsred[, "length"])
    z <- cbind(c(nsred[, "mediandev"], tumormedian), c(nsred[, 
        "length"], rep(0, length(tumormedian))), c(rep(0, nrow(nsred)), 
        1:length(tumormedian)))
    z <- z[order(z[, 1]), , drop = F]
    z[, 2] <- cumsum(z[, 2])/lnorm
    z <- z[z[, 3] != 0, , drop = F]
    negtail <- z[order(z[, 3]), 2]
    if (nisnull[1]) {
        z <- cbind(c(nsred[, "mediandev"]/nsred[, "segmad"], 
            tumormedian/tumormad), c(nsred[, "length"], rep(0, 
            length(tumormedian))), c(rep(0, nrow(nsred)), 1:length(tumormedian)))
        z <- z[order(z[, 1]), , drop = F]
        z[, 2] <- cumsum(z[, 2])/lnorm
        z <- z[z[, 3] != 0, , drop = F]
        negtailnormad <- z[order(z[, 3]), 2]
    }
    if (nisnull[2]) {
        z <- cbind(c(nsred[, "mediandev"]/nsred[, "segerr"], 
            tumormedian/tumorerror), c(nsred[, "length"], rep(0, 
            length(tumormedian))), c(rep(0, nrow(nsred)), 1:length(tumormedian)))
        z <- z[order(z[, 1]), , drop = F]
        z[, 2] <- cumsum(z[, 2])/lnorm
        z <- z[z[, 3] != 0, , drop = F]
        negtailnormerror <- z[order(z[, 3]), 2]
    }
    return(matrix(ncol = 2 + sum(nisnull), data = c(lnorm\%/\%tumorlength, 
        negtail, if (nisnull[1]) negtailnormad, if (nisnull[2]) negtailnormerror), 
        dimnames = list(NULL, c("samplesize", "negtail", if (nisnull[1]) "negtailnormad", 
            if (nisnull[2]) "negtailnormerror"))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
