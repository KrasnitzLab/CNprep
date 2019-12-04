#' CNprep: Pre-process DNA Copy Number (CN) Data for Detection of CN Events
#'
#' This package evaluates DNA copy number data, using both their initial 
#' form (copy number as a noisy function of genomic position) and their 
#' approximation by a piecewise-constant function (segmentation), for the 
#' purpose of identifying genomic regions where the copy number differs 
#' from the norm.
#'
#' @docType package
#'
#' @name CNprep-package
#'
#' @aliases CNprep-package CNprep
#'
#' @author Alexander Krasnitz and Guoli Sun
#'
#' Maintainer:
#' Guoli Sun <guolisun87@gmail.com>
#'
#' @seealso
#' \itemize{
#' \item \code{\link{CNpreprocessing}} { for pre-processing DNA copy number 
#' (CN) data for detection of CN events.}
#' \item \code{\link{makeCNPmask}} { for creating a mask given a set of
#' copy number events.}
#' \item \code{\link{applyCNPmask}} { for applying a mask to a set of
#' copy number events.}
#' }
#'
#' @keywords package
NULL


#' @name annotexample
#'
#' @docType data
#'
#' @aliases annotexample
#' 
#' @title Annotation table for ROMA CGH platform and human genome 
#' version 17.
#' 
#' @description Whole genome annotation table using Representational 
#' Oligonucleotide Microarray Analysis (ROMA) CGH platform, 
#' human genome version 17.
#'
#' @format A data frame with 83055 observations on the following 3 variables.
#' \describe{
#' \item{\code{PROBEID}}{a character vector}
#' \item{\code{CHROM}}{a numeric vector}
#' \item{\code{CHROM.POS}}{a numeric vector}
#' }
#'
#' @return A data frame with 83055 observations on the following 3 variables.
#' \describe{
#' \item{\code{PROBEID}}{a character vector}
#' \item{\code{CHROM}}{a numeric vector}
#' \item{\code{CHROM.POS}}{a numeric vector}
#' }
#'
#' @details The values in the chromosome column are all integer, with 
#' 23 corresponding to X, 24 to Y and 25 to a set of non-human test probes.
#' 
#' @seealso
#' \itemize{
#' \item \code{\link{CNpreprocessing}} {for pre-process DNA copy number (CN) 
#' data for detection of CN events.}
#' }
#'
#' @usage data(annotexample)
#'
#' @keywords datasets
#'
#' @source GEO accession GPL9775, 
#' http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL9775
#' 
#' @examples
#'
#' ## Loading annotation table dataset
#' data(annotexample)
#' 
#' ## Loading other datasets
#' data(ratexample)
#' data(segexample)
#' data(normsegs)
#' 
#' ## How to use annotexample, when segment table does not have columns 
#' ## of integer postions in terms of measuring units(probes), such 
#' ## as "mysegs" below 
#' 
#' mysegs <- segexample[,c(1,5:12)]
#' 
#' ## Analysis limited to chromosomes 1 and 2
#' ## The bsTimes variable should be higher for a real analysis
#' segtable <- CNpreprocessing(segall=mysegs, ratall=ratexample, "ID",
#'     chromCol="chrom", bpStartCol="chrom.pos.start", bpEndCol="chrom.pos.end",
#'     annot=annotexample, annotStartCol="CHROM.POS", annotEndCol="CHROM.POS",
#'     annotChromCol="CHROM", blsize=50, minJoin=0.25, cWeight=0.4, bsTimes=5,
#'     chromRange=1:2, nJobs=1, modelNames="E", 
#'     normalLength=normsegs[,1], normalMedian=normsegs[,2])
#' 
NULL


#' @title Example of a boundary positions table.
#' 
#' @description A table of genomic positions for DNA copy-number changing 
#' events, collected from genomes of 1203 individuals using Representational 
#' Oligonucleotide Microarray Analysis (ROMA) platform.
#' 
#' @name cnpexample
#'
#' @docType data
#'
#' @aliases cnpexample
#'
#' @format A data frame with 19188 rows and 4 columns. 
#' \describe{
#' \item{\code{copy.num}}{a character vector indicating whether an event is
#' a gain ("amp") or a loss ("del"). }
#' \item{\code{chrom}}{a numeric vector indicating which chromosome the 
#' event is in.}
#' \item{\code{chrom.start}}{a numeric vector of event start positions.}
#' \item{\code{chrom.end}}{a numeric vector of event start positions.}
#' }
#'
#' @return A data frame with 19188 rows and 4 columns. 
#' \describe{
#' \item{\code{copy.num}}{a character vector indicating whether an event is
#' a gain ("amp") or a loss ("del"). }
#' \item{\code{chrom}}{a numeric vector indicating which chromosome the 
#' event is in.}
#' \item{\code{chrom.start}}{a numeric vector of event start positions.}
#' \item{\code{chrom.end}}{a numeric vector of event start positions.}
#' }
#' 
#' @seealso
#' \itemize{
#' \item \code{\link{applyCNPmask}} {for applying a mask to a table of 
#' copy number events.}
#' }
#'
#' @usage data(cnpexample)
#'
#' @keywords datasets
#'
#' @source Strong association of de novo copy number mutations with autism.
#' Sebat J, Lakshmi B, Malhotra D, Troge J, Lese-Martin C, Walsh T, Yamrom B, 
#' Yoon S, Krasnitz A, Kendall J, Leotta A, Pai D, Zhang R, Lee YH, Hicks J, 
#' Bregman J, Sutcliffe JS, Jobanputra V, Chung W, Warburton D, King MC, 
#' Skuse D, Geschwind DH, Gilliam TC, Ye K, Wigler M. 
#' Science. 2007 Apr 20;316(5823):445-9.
#' 
#' @examples
#'
#' ## Loading dataset
#' data(cnpexample)
#'
#' ## Create masked table usign cnpindex vector 
#' ##namps17 <- cnpexample[cnpexample[,"copy.num"]=="amp",] 
#' ##aCNPmask <- makeCNPmask(imat=namps17, chromCol=2, startCol=3, endCol=4, 
#' ##    nProf=1203, uThresh=0.02, dThresh=0.008) 
#' ##ndels17 <- cnpexample[cnpexample[,"copy.num"]=="del",] 
#' ##dCNPmask <- makeCNPmask(imat=ndels17, chromCol=2, startCol=3, endCol=4, 
#' ##    nProf=1203, uThresh=0.02, dThresh=0.008) 
#' ##cnptable <- rbind(cbind(aCNPmask,cnpindex=1),cbind(dCNPmask,cnpindex=-1))
#' 
#' # TODO
#' ## Run the CNP test using masked table
#' #myCNPtable <- applyCNPmask(segtable,"chrom",startPos="chrom.pos.start", 
#' #    endPos="chrom.pos.end","start","end","eventIndex",masktable=cnptable,
#' #    "chrom", maskstart="start",maskend="end",maskindex="cnpindex",
#' #    mincover=0.005,indexvals=c(-1,1))
#' 
NULL

#' @title Annotation table for ROMA CGH platform and human genome version 17.
#' 
#' @description Whole genome annotation table using Representational 
#' Oligonucleotide Microarray Analysis (ROMA) CGH platform, 
#' human genome version 17.
#' 
#' @name annotexample
#'
#' @docType data
#'
#' @aliases annotexample
#'
#' @format A \code{data.frame} with 83055 observations on the following 3 
#' variables.
#' \describe{
#' \item{\code{PROBEID}}{a \code{character} \code{vector} of probe names.}
#' \item{\code{CHROM}}{a \code{numeric} \code{vector} of chromosome positions.}
#' \item{\code{CHROM.POS}}{a \code{numeric} \code{vector} of genomic 
#' positions.}
#' }
#'
#' @return a \code{data.frame} with 83055 observations on the following 3 
#' variables. 
#' \describe{
#' \item{\code{PROBEID}}{a \code{character} \code{vector} of probe names.}
#' \item{\code{CHROM}}{a \code{numeric} \code{vector} of chromosome positions.}
#' \item{\code{CHROM.POS}}{a \code{numeric} \code{vector} of genomic 
#' positions.}
#' }
#' 
#' @details The values in the chromosome column are all integer, with 
#' 23 corresponding to X, 24 to Y and 25 to a set of non-human test probes.
#' 
#' @seealso
#' \itemize{
#' \item \code{\link{CNpreprocessing}} {for pre-process DNA copy number (CN) 
#' data for detection of CN events.}
#' }
#'
#' @usage data(annotexample)
#'
#' @keywords datasets
#'
#' @source GEO accession GPL9775, 
#' http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL9775
#' 
#' @examples
#'
#' ## Loading annotation table dataset
#' data(annotexample)
#'
#' ## Loading other datasets
#' data(ratexample)
#' data(segexample)
#' data(normsegs)
#' 
#' ## How to use annotexample, when segment table does not have columns 
#' ## of integer postions in terms of measuring units(probes), such 
#' ## as "mysegs" below 
#' 
#' mysegs <- segexample[,c(1,5:12)]
#' 
#' ## Analysis limited to chromosomes 1 and 2
#' ## The bsTimes variable should be higher for a real analysis
#' segtable <- CNpreprocessing(segall=mysegs, ratall=ratexample, idCol="ID",
#'     chromCol="chrom", bpStartCol="chrom.pos.start", bpEndCol="chrom.pos.end",
#'     annot=annotexample, annotStartCol="CHROM.POS", annotEndCol="CHROM.POS",
#'     annotChromCol="CHROM", blsize=50, minJoin=0.25, cWeight=0.4, bsTimes=3,
#'     chromRange=1:2, nJobs=1, modelNames="E", 
#'     normalLength=normsegs[,1], normalMedian=normsegs[,2])
#' 
NULL


#' @title Example of copy number log ratio dataset
#' 
#' @description Log ratio data for 5 breast cancer genomes, derived using 
#' Representational Oligonucleotide Microarray Analysis (ROMA), an array-based 
#' hybridization method that uses genomic complexity reduction based 
#' on representations.
#' 
#' @name ratexample
#'
#' @docType data
#'
#' @aliases ratexample
#'
#' @format a log ratio \code{matrix} with 83055 rows, one per 
#' oligonucleotide probe, and 5 columns, one for each breast tumor sample.
#'
#' @return a log ratio \code{matrix} with 83055 rows, one per 
#' oligonucleotide probe, and 5 columns, one for each breast tumor sample.
#' 
#' @details The values are natural log copy number ratios, consistent with 
#' data in \code{segexample} (segmented data for these tumors) and 
#' \code{normsegs}. These copy number ratios are normalized using an 
#' intensity-based lowess curve fitting algorithm.
#' 
#' @seealso
#' \itemize{
#' \item \code{\link{CNpreprocessing}} {for pre-process DNA copy number (CN) 
#' data for detection of CN events.}
#' }
#'
#' @usage data(ratexample)
#'
#' @keywords datasets
#' 
#' @source Hicks, J. et al. Novel patterns of genome rearrangement and their 
#' association with survival in breast cancer. Genome Res. 2006. 16:1465–1479.
#' doi: 10.1101/gr.5460106
#' 
#' @examples
#'
#' ## Loading log ratio dataset
#' data(ratexample)
#'
#' ## Plot the whole genome log ratio data for the first profile "WZ1"
#' ## Note X and Y chromosomes at the far right of the plot
#' plot(ratexample[,"WZ1"], ylab="log ratio", xlab="position")
#' 
NULL


