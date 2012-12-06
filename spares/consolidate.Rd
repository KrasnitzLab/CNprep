\name{consolidate}
\alias{consolidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
consolidate(emfit, minover)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{emfit}{
%%     ~~Describe \code{emfit} here~~
}
  \item{minover}{
%%     ~~Describe \code{minover} here~~
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
function (emfit, minover) 
{
    newem <- list(mu = emfit$parameters$mean, pro = emfit$parameters$pro, 
        z = emfit$z, groups = matrix(nrow = length(emfit$parameters$mean), 
            ncol = length(emfit$parameters$mean), data = 0), 
        ngroups = length(emfit$parameters$mean))
    diag(newem$groups) <- 1
    while (newem$ngroups > 1) {
        avz <- (t(newem$z)/colSums(newem$z)) \%*\% newem$z
        diag(avz) <- 0
        if (max(avz) > minover) {
            g1 <- col(avz)[rev(order(avz))][1]
            g2 <- row(avz)[rev(order(avz))][1]
            gl <- min(g1, g2)
            gr <- max(g1, g2)
            newem$z[, gl] <- newem$z[, gl] + newem$z[, gr]
            newem$z <- matrix(ncol = ncol(newem$z) - 1, nrow = nrow(newem$z), 
                data = newem$z[, -gr])
            newem$mu[gl] <- (newem$mu[gl] * newem$pro[gl] + newem$mu[gr] * 
                newem$pro[gr])/(newem$pro[gl] + newem$pro[gr])
            newem$mu <- newem$mu[-gr]
            newem$pro[gl] <- newem$pro[gl] + newem$pro[gr]
            newem$pro <- newem$pro[-gr]
            newem$groups[gl, ] <- newem$groups[gl, ] + newem$groups[gr, 
                ]
            newem$groups <- matrix(ncol = ncol(newem$groups), 
                nrow = newem$ngroups - 1, data = newem$groups[-gr, 
                  ])
            newem$ngroups <- newem$ngroups - 1
        }
        else break
    }
    return(newem)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
