\name{segsample}
\alias{segsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segsample(mysegs, ratcol, startcol = "StartProbe", endcol = "EndProbe", blocksize = 0, times = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mysegs}{
%%     ~~Describe \code{mysegs} here~~
}
  \item{ratcol}{
%%     ~~Describe \code{ratcol} here~~
}
  \item{startcol}{
%%     ~~Describe \code{startcol} here~~
}
  \item{endcol}{
%%     ~~Describe \code{endcol} here~~
}
  \item{blocksize}{
%%     ~~Describe \code{blocksize} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
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
function (mysegs, ratcol, startcol = "StartProbe", endcol = "EndProbe", 
    blocksize = 0, times = 0) 
{
    if (blocksize == 0 & times == 0) 
        stop("One of blocksize or times must be set")
    if (blocksize != 0 & times != 0) 
        stop("Only one of blocksize or times can be set")
    segtable <- mysegs[, c(startcol, endcol), drop = F]
    if (blocksize != 0) 
        segtable <- segtable[rep(1:nrow(segtable), times = (segtable[, 
            endcol] - segtable[, startcol] + 1)\%/\%blocksize), 
            ]
    if (times != 0) 
        segtable <- segtable[rep(1:nrow(segtable), each = times), 
            ]
    return(cbind(segtable, apply(segtable, 1, smedian.sample, 
        v = ratcol)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
