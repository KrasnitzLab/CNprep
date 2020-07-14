#' @title Parameters validation for the \code{\link{CNpreprocessing}} function
#' 
#' @description Validation of all parameters needed by the public
#' \code{\link{CNpreprocessing}} function.
#' 
#' @param ratall A \code{matrix} whose rows correspond to genomic positions 
#' and columns to copy number profiles. Its matrix elements are functions of 
#' copy number, most often log ratios of copy number to the expected standard 
#' value, such as 2 in diploid genomes.
#' 
#' @param idCol A \code{character} string specifying the name for the 
#' column in \code{segall} tabulating the profile IDs. When not specified,
#' the numerical column of the \code{ratall} object will be used as the
#' profile IDs.
#' 
#' @param startCol A \code{character} string specifying the name of column 
#' in \code{segall} that tabulates the (integer) start postion of each segment 
#' in internal units such as probe numbers for data of CGH microarray origin.
#' 
#' @param endCol A \code{character} string specifying the name of column 
#' in \code{segall} that tabulates the (integer) end postion of each segment 
#' in internal units such as probe numbers for data of CGH microarray origin.
#' 
#' @param medCol A \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' the function (\code{medCol}), a measure of its spread (\code{madCol}) and 
#' its error (\code{errorCol}) for the segment.
#' 
#' @param madCol A \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' a measure of spread (\code{madCol}) related to  
#' the function (\code{medCol}) for the segment.
#' 
#' @param errorCol A \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' the error (\code{errorCol}) related to  
#' the function (\code{medCol}) for the segment. 
#' 
#' @param chromCol A \code{character} string specifying the name for the 
#' column in \code{segall} tabulating the (integer) chromosome number for 
#' each segment.
#' 
#' @param bpStartCol A \code{character} string specifying the name of 
#' column in \code{segall} that tabulates the (integer) genomic start 
#' coordinate of each segment.
#' 
#' @param bpEndCol A \code{character} string specifying the name of 
#' column in \code{segall} that tabulates the (integer) genomic end 
#' coordinate of each segment.
#' 
#' @param annot A matrix or a \code{data.frame} that contains the annotation 
#' for the copy number measurement platform in the study. It is generally 
#' expected to contain columns with names specified by 
#' \code{annotStartCol, annotEndCol, annotChromCol}.
#' 
#' @param annotStartCol A \code{character} string 
#' specifying the name of column in \code{annot} that tabulates the (integer) 
#' genomic start coordinates in case of CGH
#' microarrays.
#' 
#' @param annotEndCol A \code{character} string 
#' specifying the name of column in \code{annot} that tabulates the (integer) 
#' genomic end coordinates in case of CGH
#' microarrays.
#' 
#' @param annotChromCol A \code{character} string 
#' specifying the name of column in \code{annot} that tabulates the chromosome
#' number for each copy number measuring unit, such as a probe in case of CGH
#' microarrays.
#' 
#' @param useEnd A single \code{logical} value specifying whether the segment 
#' end positions as given by the \code{bpEndCol} of \code{segall} are to be 
#' looked up in the \code{annotEndCol} column of \code{annot} 
#' (if \code{useEnd=TRUE}) or in the \code{annotStartCol} column (default). 
#' 
#' @param blsize A single \code{integer} specifying the bootstrap sampling 
#' rate of segment medians to generate input for model-based clustering. The 
#' number of times a segment is sampled is then given by the (integer) 
#' division of the segment length in internal units by \code{blsize}.
#' 
#' @param minJoin A single numeric value between 0 and 1 specifying the 
#' degree of overlap above which two clusters will be joined into one. 
#' 
#' @param nTrial A single positive \code{integer} specifying the number of 
#' times a model-based 
#' clustering is attempted for each profile in order to achieve the 
#' highest Bayesian information criterion (BIC). 
#' 
#' @param bestBIC A single \code{numeric} value for initalizing BIC 
#' maximization. A large negative value is recommended. 
#' 
#' @param modelNames A \code{vector} of \code{character} strings specifying 
#' the names of models to be used in model-based clustering (see package 
#' \code{mclust} for further details).
#' 
#' @param cWeight A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the minimal share of the central cluster in each profile.
#' 
#' @param bsTimes A single positive \code{double} value specifying the number 
#' of time the median of each segment is sampled in order to predict the 
#' cluster assignment for the segment.
#' 
#' @param chromRange A \code{integer} \code{vector} enumerating chromosomes 
#' from which segments are to be used for initial model-based clustering.
#' 
#' @param nJobs a single positive \code{integer} specifying the number of 
#' worker jobs to create in case of distributed computation. 
#' 
#' @param normalLength An integer \code{vector} specifying the genomic lengths 
#' of segments in the normal reference data. 
#' 
#' @param normalMedian A numeric \code{vector}, 
#' of the same length as \code{normalLength}, specifying the segment values
#' of the normal reference segments.
#' 
#' @param normalMad A numeric \code{vector}, 
#' of the same length as \code{normalLength}, specifying the value spreads 
#' of the normal reference segments. 
#' 
#' @param normalError A numeric \code{vector}, 
#' of the same length as \code{normalLength}, specifying the error values
#' of the normal reference segments.
#' 
#' @return \code{0}. 
#' 
#' @examples
#'
#' data(segexample)
#' data(ratexample)
#' data(normsegs)
#' 
#' ## Return zero as all parameters are valid
#' CNprep:::validateCNpreprocessing(segall=segexample,
#'     ratall=ratexample, idCol="ID", startCol="start", endCol="end", 
#'     chromCol="chrom", bpStartCol="chrom.pos.start", 
#'     bpEndCol="chrom.pos.end", blsize=50, nTrial=10,
#'     useEnd=FALSE, minJoin=0.25, cWeight=0.4, bsTimes=50, chromRange=1:3, 
#'     nJobs=1, modelNames="E", normalLength=normsegs[,1],
#'     normalMedian=normsegs[,2])
#' 
#' @author Astrid Deschênes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateCNpreprocessing <- function(segall, ratall, idCol, 
        startCol, endCol, medCol, madCol, errorCol, 
        chromCol, bpStartCol, bpEndCol, annot, 
        annotStartCol, 
        annotEndCol, annotChromCol, useEnd, blsize, 
        minJoin, nTrial, bestBIC, modelNames, cWeight,
        bsTimes, chromRange, nJobs, normalLength, 
        normalMedian, normalMad,
        normalError, weightall) 
{
    ## Validate that nJobs is an positive integer
    if (!(isSingleInteger(nJobs) || isSingleNumber(nJobs)) ||
        as.integer(nJobs) < 1) {
        stop("nJobs must be a positive integer")
    }

    ## Validate that nTrial is an positive integer
    if (!(isSingleInteger(nTrial) || isSingleNumber(nTrial)) ||
        as.integer(nTrial) < 1) {
        stop("nTrial must be a positive integer")
    }
    
    ## Validate that nJobs is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nJobs) != 1) {
       stop("nJobs must be 1 on a Windows system")
    }
    
    ## Validate that useEnd is logical
    if (!is.logical(useEnd)) {
        stop("useEnd must be a logical value")
    }
    
    ## Validate that bsTimes is an positive integer
    if (!(isSingleInteger(bsTimes) || isSingleNumber(bsTimes)) ||
        as.integer(bsTimes) < 1) {
        stop("bsTimes must be a positive integer")
    }
    
    return(0L)
}