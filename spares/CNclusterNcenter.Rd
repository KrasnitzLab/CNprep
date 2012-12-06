\name{CNclusterNcenter}
\alias{CNclusterNcenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CNclusterNcenter(segrat, blsize, minjoin, ntrial, bestbic, modelNames, cweight, bstimes, chromrange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segrat}{
%%     ~~Describe \code{segrat} here~~
}
  \item{blsize}{
%%     ~~Describe \code{blsize} here~~
}
  \item{minjoin}{
%%     ~~Describe \code{minjoin} here~~
}
  \item{ntrial}{
%%     ~~Describe \code{ntrial} here~~
}
  \item{bestbic}{
%%     ~~Describe \code{bestbic} here~~
}
  \item{modelNames}{
%%     ~~Describe \code{modelNames} here~~
}
  \item{cweight}{
%%     ~~Describe \code{cweight} here~~
}
  \item{bstimes}{
%%     ~~Describe \code{bstimes} here~~
}
  \item{chromrange}{
%%     ~~Describe \code{chromrange} here~~
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
function (segrat, blsize, minjoin, ntrial, bestbic, modelNames, 
    cweight, bstimes, chromrange) 
{
    RNGkind("L'Ecuyer-CMRG")
    .Random.seed <- segrat$rngseed
    startcol <- "StartProbe"
    endcol <- "EndProbe"
    chromcol <- "chrom"
    medcol <- "segmedian"
    madcol <- "segmad"
    segrat$seg <- cbind(segrat$seg, t(apply(segrat$seg[, c(startcol, 
        endcol, chromcol), drop = F], 1, smedmad, v = segrat$rat)))
    dimnames(segrat$seg)[[2]] <- c(startcol, endcol, chromcol, 
        medcol, madcol)
    seguse <- segrat$seg[segrat$seg[, chromcol] \%in\% chromrange, 
        , drop = F]
    aux <- rep(0, length(segrat$rat))
    aux[seguse[, startcol]] <- 1
    aux[seguse[, endcol]] <- (-1)
    aux <- cumsum(aux)
    aux[seguse[, endcol]] <- 1
    ratuse <- segrat$rat[aux == 1]
    for (j in 1:ntrial) {
        aaa <- segsample(seguse, ratuse, blocksize = blsize)
        emfit <- Mclust(aaa[, 3], maxG = 15, modelNames = modelNames)
        if (emfit$bic >= bestbic) {
            bestaaa <- aaa
            bestem <- emfit
            bestbic <- emfit$bic
        }
    }
    newem <- consolidate(bestem, minjoin)
    newem <- get.center(newem, cweight)
    profcenter <- weighted.median(bestaaa[, 3], newem$z[, newem$center])
    mediandev <- segrat$seg[, medcol] - profcenter
    segs <- segsample(segrat$seg, segrat$rat, times = bstimes)
    segz <- getz(segs[, 3], bestem, newem$groups, times = bstimes)[, 
        newem$center]
    cpb <- centerprob(segs[, 3], bestem, newem$groups, times = bstimes, 
        newem$center)
    w <- t(matrix(nrow = bstimes, data = segs[, 3]))
    segerr <- sqrt(apply(w, 1, var, na.rm = T))
    return(cbind(segrat$seg[, medcol], segrat$seg[, madcol], 
        mediandev, segerr, segz, cpb))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
