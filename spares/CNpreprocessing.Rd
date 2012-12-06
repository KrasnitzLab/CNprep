\name{CNpreprocessing}
\alias{CNpreprocessing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pre-process DNA copy number (CN) data for detection of CN events.
}
\description{
Description: The package evaluates DNA copy number data, using both their 
intitial form (copy number as a noisy function of genomic position) and their
approximation by a piecewise-constant function (segmentation), for the purpose
of identifying genomic regions where the copy number differs from the norm.
}
\usage{
CNpreprocessing(segall, ratall = NULL, idcol = NULL, startcol = NULL, endcol = NULL, medcol = NULL, madcol = NULL, errorcol = NULL, chromcol = NULL, bpstartcol = NULL, bpendcol = NULL, annot = NULL, annotstartcol = NULL, annotendcol = NULL, annotchromcol = NULL, useend = F, blsize = NULL, minjoin = NULL, ntrial = 10, bestbic = -1e+07, modelNames = "E", cweight = NULL, bstimes = NULL, chromrange = NULL, myseed = NULL, distrib = c("vanilla", "Rparallel", "sge"), njobs = 1, normalength = NULL, normalmedian = NULL, normalmad = NULL, normalerror = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segall}{
A matrix or a data frame for segmented copy number profiles. It may have a 
character column, with a name specified by \code{idcol}, and/or numeric 
columns with names specified by \code{startcol, endcol, medcol, madcol,
errorcol, chromcol, bpstartcol, bpendcol}. Each row of \code{segall} 
corresponds to a segment belonging to one of the profiles to be pre-processed.
}
  \item{ratall}{
A matrix whose rows correspond to genomic positions and columns to copy number
profiles. Its matrix elements are functions of copy number, most often
log ratios of copy number to the expected standard value, such as 2 in diploid 
genomes.
}
  \item{idcol}{
A character string specifying the name for the column in \code{segall} 
tabulating the profile IDs.
}
  \item{startcol, endcol}{
Character strings specifying the names of columns in \code{segall} that 
tabulate the (integer) start and end postions of each segment in internal 
units such as probe numbers for data of CGH microarray origin.
}
  \item{medcol, madcol, errorcol}{
Character strings specifying the names of columns in \code{segall} that,
for the function of copy number used in the study (typically log ratios),
tabulate the (numeric) values for the function (\code{medcol}), a measure
of its spread (\code{madcol}) and its error (\code{errorcol}) for the segment.
}
  \item{chromcol}{
A character string specifying the name for the column in \code{segall} 
tabulating the (integer) chromosome number for each segment.
}
  \item{bpstartcol, bpendcol}{
Character strings specifying the names of columns in \code{segall} that 
tabulate the (integer) genomic start and end coordinates of each segment.
}
  \item{annot}{
A matrix or a data frame that contains the annotation for the copy number 
measurement platform in the study. It is generally expected to contain columns
with names specified by \code{annotstartcol, annotendcol, annotchromcol}.
}
  \item{annotstartcol, annotendcol, annotchromcol}{
Character strings specifying the names of columns in \code{annot} that 
tabulate the (integer) genomic start and end coordinates and the chromosome
number for each copy number measuring unit, such as a probe in case of CGH
microarrays.
}
  \item{useend}{
A single logical value specifying whether the segment end positions as given by
the \code{bpendcol} of \code{segall} are to be looked up in the 
\code{annotendcol} column of \code{annot} (if \code{useend=TRUE}) or in the
\code{annotstartcol} column (default).
}
  \item{blsize}{
A single integer specifying the bootstrap sampling rate of segment medians to 
generate input for model-based clustering. The number of times a segment is 
sampled is then given by the (integer) division of the segment length in 
internal units by \code{blsize}.
}
  \item{minjoin}{
A single numeric value between 0 and 1 specifying the degree of overlap above 
which two clusters will be joined into one. 
}
  \item{ntrial}{
A single integer specifying the number of times a model-based clustering is
attempted for each profile in order to achieve the highest Bayesian information
criterion (BIC).
}
  \item{bestbic}{
A single numeric value for initalizing BIC maximization. A large negative value
is recommended. The default is \code{-1e7}.
}
  \item{modelNames}{
A vector of character strings specifying the names of models to be used in
model-based clustering (see package \code{mclust} for further details). The
default is \code{"E"}.
}
  \item{cweight}{
A single numeric value between 0 and 1 specifying the minimal share of the 
central cluster in each profile.
}
  \item{bstimes}{
A single integer value specifying the number of time the median of each segment
is sampled in order to predict the cluster assignment for the segment.
}
  \item{chromrange}{
A numeric vector enumerating chromosomes from which segments are to be used for
initial model-based clustering.
}
  \item{myseed}{
A single integer value to seed the random number generator.
}
  \item{distrib}{
One of \code{"vanilla", "Rparallel"} or {"sge"} to specify the distributed
computing option for the cluster assignment step. For \code{"vanilla"} (default)
no distributed computing is performed. For \code{"Rparallel"} the 
\code{parallel} package of \code{R} core is used for multi-core processing. 
For \code{"sge"} SGE is invoked using CRAN contributed package \code{Rsge}.
}
  \item{njobs}{
A single integer specifying the number of worker jobs to create in case of 
distributed computation.
}
  \item{normalength}{
An integer vector specifying the genomic lengths of segments in the normal 
reference data.
}
  \item{normalmedian, normalmad, normalerror}{
Numeric vectors, of the same length as \code{normalength}, specifying the 
segment values, value spreads and errors of the normal reference segments.
}
}
\details{
Depending on the availability of input, the function will perform the following
operations for each copy number profile.

If raw data are available in addition to segment start and end positions,
median and MAD of each segment will be computed. For each profile, bootstrap 
sampling of the segment median values will be performed, and the sample will be
used to estimate the error in the median for each segment. Model-dependent
clustering (fitting to a gaussian mixture) of the sample will be performed. 
The central cluster (the one nearest the expected unaltered value) will be 
identified and, if necessary, merged with adjacent clusters in order to
comprise the minimal required fraction of the data. Deviation of each segment
from the center, its probabilty to belong to the central cluster and its 
marginal probability in the central cluster will be computed.

If segment medians or median deviations are available or have been computed, 
and, in addition, genomic lengths and average values are given for a 
collection of segments with unaltered copy number, additional estimates will be
performed. If median values are available for the unaltered segments, the
marginal probability of the observed median or median deviation in the unaltered
set will be computed for each segment. Likewise, marginal probabilities for
median/MAD and/or median/error will be computed if these statistics are 
available. 
}
\value{
The input \code{segall} data frame to which some or all of the following columns
may be bound, depending on the availability of input:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{segmedian}{Median function of copy number}
\item{segmad}{MAD for the function of copy number}
\item{mediandev}{median function of copy number relative to its central value}
\item{segerr}{error estimate for the function of copy number}
\item{segz}{the probability that the segment is in the central cluster}
\item{marginalprob}{marginal probability for the segment in the central cluster}
\item{negtail}{the probability of finding the deviation as observed or larger 
in a collection of central segments}
\item{negtailnormad}{the probability of finding the deviation/MAD as observed 
or larger in a collection of central segments}
\item{negtailnormerror}{the probability of finding the deviation/error as 
observed or larger in a collection of central segments}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Alex Krasnitz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (segall, ratall = NULL, idcol = NULL, startcol = NULL, 
    endcol = NULL, medcol = NULL, madcol = NULL, errorcol = NULL, 
    chromcol = NULL, bpstartcol = NULL, bpendcol = NULL, annot = NULL, 
    annotstartcol = NULL, annotendcol = NULL, annotchromcol = NULL, 
    useend = F, blsize = NULL, minjoin = NULL, ntrial = 10, bestbic = -1e+07, 
    modelNames = "E", cweight = NULL, bstimes = NULL, chromrange = NULL, 
    myseed = NULL, distrib = c("vanilla", "Rparallel", "sge"), 
    njobs = 1, normalength = NULL, normalmedian = NULL, normalmad = NULL, 
    normalerror = NULL) 
{
    if (is.null(idcol)) {
        cat("Found a single segmented profile with no ID", "\n")
        if (!is.null(ratall)) {
            if (sum(apply(ratall, 2, data.class) == "numeric") > 
                1) 
                stop("Ambiguity: more than 1 numeric column in raw data table\n")
            else {
                idrat <- which(apply(ratall, 2, data.class) == 
                  "numeric")
                segall <- data.frame(rep(as.character(idrat), 
                  nrow(segall)), segall)
                idcol <- "ID"
                dimnames(segall)[[2]][1] <- idcol
            }
        }
    }
    if (is.null(ratall)) 
        cat("No raw table, proceeding to comparison\n")
    else {
        profnames <- unique(segall[, idcol])
        if (!all(profnames \%in\% dimnames(ratall)[[2]])) 
            stop("Found unmatched segmented profile IDs\n")
        if (is.null(startcol) | is.null(endcol)) {
            if (is.null(bpstartcol) | is.null(bpendcol) | is.null(chromcol)) 
                stop("Unable to proceed: incomplete segment annotation\n")
            if (is.null(chromrange)) 
                chromrange <- sort(unique(segall[, chromcol]))
            if (is.null(annot)) 
                stop("No annotation table; unable to determine boundary probes/bins\n")
            if (is.null(annotstartcol) | is.null(annotchromcol)) 
                stop("No start and chrom column names provided for annotation table\n")
            if (useend & is.null(annotendcol)) 
                stop("End column name required but not provided in annotation table\n")
            maxbpstart <- max(c(segall[, bpstartcol], annot[, 
                annotstartcol])) + 1
            maxbpend <- ifelse(useend, max(c(segall[, bpendcol], 
                annot[, annotendcol])), max(c(segall[, bpendcol], 
                annot[, annotstartcol]))) + 1
            startprobe <- match((segall[, chromcol] - 1) * maxbpstart + 
                segall[, bpstartcol], (annot[, annotchromcol] - 
                1) * maxbpstart + annot[, annotstartcol])
            endprobe <- ifelse(useend, match((segall[, chromcol] - 
                1) * maxbpend + segall[, bpendcol], (annot[, 
                annotchromcol] - 1) * maxbpend + annot[, annotendcol]), 
                match((segall[, chromcol] - 1) * maxbpend + segall[, 
                  bpendcol], (annot[, annotchromcol] - 1) * maxbpend + 
                  annot[, annotstartcol]))
            if (!all(!is.na(startprobe) & !is.na(endprobe))) 
                stop("Incomplete start and end annotation of segments\n")
            segall <- data.frame(segall, startprobe, endprobe)
            startcol <- "StartProbe"
            endcol <- "EndProbe"
        }
        profpack <- vector(mode = "list", length = length(profnames))
        names(profpack) <- profnames
        RNGkind("L'Ecuyer-CMRG")
        if (!is.null(myseed)) 
            set.seed(myseed)
        s <- .Random.seed
        distrib <- match.arg(distrib)
        for (pn in profnames) {
            profpack[[pn]] <- vector(mode = "list", length = 3)
            names(profpack[[pn]]) <- c("seg", "rat", "rngseed")
            profpack[[pn]]$seg <- segall[segall[, idcol] == pn, 
                c(startcol, endcol, chromcol), drop = F]
            dimnames(profpack[[pn]]$seg)[[2]] <- c("StartProbe", 
                "EndProbe", "chrom")
            profpack[[pn]]$rat <- ratall[, pn]
            s <- nextRNGStream(s)
            profpack[[pn]]$rngseed <- s
        }
        distrib <- match.arg(distrib)
        if (distrib == "Rparallel") {
            ncores <- min(njobs, length(profnames), detectCores())
            cl <- parallel::makeCluster(getOption("cl.cores", 
                ncores))
            parallel::clusterEvalQ(cl = cl, expr = library(parallel))
            parallel::clusterEvalQ(cl = cl, expr = library(mclust))
        }
        if (distrib == "sge") 
            ncores <- min(njobs, length(profnames))
        processed <- switch(distrib, vanilla = lapply(X = profpack, 
            FUN = CNclusterNcenter, blsize = blsize, minjoin = minjoin, 
            ntrial = ntrial, bestbic = bestbic, modelNames = modelNames, 
            cweight = cweight, bstimes = bstimes, chromrange = chromrange), 
            Rparallel = parLapply(cl, X = profpack, fun = CNclusterNcenter, 
                blsize = blsize, minjoin = minjoin, ntrial = ntrial, 
                bestbic = bestbic, modelNames = modelNames, cweight = cweight, 
                bstimes = bstimes, chromrange = chromrange), 
            sge = sge.parLapply(cl, X = profpack, FUN = CNclusterNcenter, 
                blsize = blsize, minjoin = minjoin, ntrial = ntrial, 
                bestbic = bestbic, modelNames = modelNames, cweight = cweight, 
                bstimes = bstimes, chromrange = chromrange, packages = c("parallel", 
                  "mclust")))
        segall <- cbind(segall, t(matrix(nrow = 6, data = unlist(lapply(processed, 
            t)))))
        dimnames(segall)[[2]][(ncol(segall) - 5):ncol(segall)] <- c("segmedian", 
            "segmad", "mediandev", "segerr", "segz", "marginalprob")
        medcol <- "mediandev"
        madcol <- "segmad"
        errorcol <- "segerr"
    }
    if (!(is.null(normalength) | is.null(normalmedian) | is.null(medcol))) {
        if (is.null(bpstartcol) | is.null(bpendcol)) {
            if (is.null(startcol) | is.null(endcol) | is.null(annot) | 
                is.null(annotstartcol) | is.null(annotendcol)) 
                stop("Insufficient annotation for comparison")
            tumorlength <- annot[segall[, endcol], annotendcol] - 
                annot[segall[, startcol], annotstartcol] + 1
        }
        else tumorlength <- segall[, bpendcol] - segall[, bpstartcol] + 
            1
        tumormedian <- segall[, medcol]
        if (!is.null(madcol)) 
            tumormad <- segall[, madcol]
        if (!is.null(errorcol)) 
            tumorerror <- segall[, errorcol]
        segall <- cbind(segall, normalComparison(normalmedian, 
            normalength, tumormedian, tumorlength, normalmad, 
            normalerror, tumormad, tumorerror))
    }
    return(segall)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
