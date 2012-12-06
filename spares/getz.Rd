\name{getz}
\alias{getz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getz(logr, emfit, zgroup, times)
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
function (logr, emfit, zgroup, times) 
{
    gz <- predict(emfit, newdata = logr)$z
    gz[is.nan(gz)] <- 0
    gz <- gz \%*\% t(zgroup)
    gz <- matrix(ncol = ncol(gz), data = apply(gz, 2, cumsum)[seq(from = times, 
        to = nrow(gz), by = times), ]/times)
    return(gz - rbind(matrix(nrow = 1, data = rep(0, ncol(gz))), 
        matrix(ncol = ncol(gz), data = gz[-nrow(gz), ])))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
