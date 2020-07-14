#' @title Pre-process DNA copy number (CN) data for detection of CN events.
#' 
#' @description The package evaluates DNA copy number data, using both their 
#' initial form (copy number as a noisy function of genomic position) and their
#' approximation by a piecewise-constant function (segmentation), for the 
#' purpose of identifying genomic regions where the copy number differs from 
#' the norm.
#' 
#' @param segall a \code{matrix} or a \code{data.frame} for segmented copy 
#' number profiles. It may have a character column, with a name specified 
#' by \code{idCol}, and/or numeric columns with names specified by 
#' \code{startCol, endCol, medCol, madCol,errorCol}  
#' \code{,chromCol, bpStartCol, bpEndCol}. Each row of \code{segall} 
#' corresponds to a segment belonging to one of the profiles 
#' to be pre-processed.
#' 
#' @param ratall a \code{matrix} whose rows correspond to genomic positions 
#' and columns to copy number profiles. The elements of this matrix are 
#' functions of copy number, most often log ratios of copy number to 
#' the expected standard value, such as 2 in diploid genomes. 
#' 
#' @param idCol a \code{character} string specifying the name for the 
#' column in \code{segall} tabulating the profile IDs. When not specified,
#' the numerical column of the \code{ratall} object will be used as the
#' profile IDs. Default: \code{NULL}.
#' 
#' @param startCol a \code{character} string specifying the name of column 
#' in \code{segall} that tabulates the (integer) start position of each segment 
#' in internal units such as probe numbers for data of CGH microarray origin.
#' Default: \code{NULL}.
#' 
#' @param endCol a \code{character} string specifying the name of column 
#' in \code{segall} that tabulates the (integer) end position of each segment 
#' in internal units such as probe numbers for data of CGH microarray origin.
#' Default: \code{NULL}.
#' 
#' @param medCol a \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' the function (\code{medCol}), a measure of its spread (\code{madCol}) and 
#' its error (\code{errorCol}) for the segment. Default: \code{NULL}.
#' 
#' @param madCol a \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' a measure of spread (\code{madCol}) related to  
#' the function (\code{medCol}) for the segment. Default: \code{NULL}.
#' 
#' @param errorCol a \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' the error (\code{errorCol}) related to  
#' the function (\code{medCol}) for the segment. Default: \code{NULL}.
#' 
#' @param chromCol a \code{character} string specifying the name for the 
#' column in \code{segall} tabulating the (integer) chromosome number for 
#' each segment.
#' 
#' @param bpStartCol a \code{character} string specifying the name of 
#' column in \code{segall} that tabulates the (integer) genomic start 
#' coordinate of each segment.
#' 
#' @param bpEndCol a \code{character} string specifying the name of 
#' column in \code{segall} that tabulates the (integer) genomic end 
#' coordinate of each segment.
#' 
#' @param annot a \code{matrix} or a \code{data.frame} that contains the 
#' annotation for the copy number measurement platform in the study. It is 
#' generally expected to contain columns with names specified by 
#' \code{annotStartCol, annotEndCol, annotChromCol}.
#' 
#' @param annotStartCol a \code{character} string 
#' specifying the name of column in \code{annot} that tabulates the (integer) 
#' genomic start coordinates in case of CGH microarrays.
#' 
#' @param annotEndCol a \code{character} string 
#' specifying the name of column in \code{annot} that tabulates the (integer) 
#' genomic end coordinates in case of CGH microarrays.
#' 
#' @param annotChromCol a \code{character} string 
#' specifying the name of column in \code{annot} that tabulates the chromosome
#' number for each copy number measuring unit, such as a probe in case of CGH
#' microarrays.
#' 
#' @param useEnd a single logical value specifying whether the segment end 
#' positions as given by the \code{bpEndCol} of \code{segall} are to be 
#' looked up in the \code{annotEndCol} column of \code{annot} 
#' (if \code{useEnd=TRUE}) or in the \code{annotStartCol} column (default). 
#' Default: \code{FALSE}.
#' 
#' @param blsize a single \code{integer} specifying the bootstrap sampling 
#' rate of segment medians to generate input for model-based clustering. The 
#' number of times a segment is sampled is then given by the (integer) 
#' division of the segment length in internal units by \code{blsize}.
#' 
#' @param minJoin a single \code{numeric} value between 0 and 1 specifying the 
#' degree of overlap above which two clusters will be joined into one. Default:
#' \code{NULL}. TODO= HAVE a default value, not NULL.
#' 
#' @param nTrial a single positive \code{integer} specifying the number of 
#' times a model-based 
#' clustering is attempted for each profile in order to achieve the 
#' highest Bayesian information criterion (BIC). Default: \code{10}.
#' 
#' @param bestBIC a single \code{numeric} value for initalizing the
#' Bayesian information criterion (BIC) 
#' maximization. A large negative value is recommended. Default: \code{-1e7}.
#' 
#' @param modelNames a \code{vector} of \code{character} strings specifying 
#' the names of models to be used in model-based clustering (see package 
#' \code{mclust} for further details). The default is \code{"E"}.
#' 
#' @param cWeight A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the minimal share of the central cluster in each profile.
#' 
#' @param bsTimes a single positive \code{double} value specifying the number 
#' of time the median of each segment is sampled in order to predict the 
#' cluster assignment for the segment. Default: \code{NULL}. TODO: select a 
#' default value that is not null.
#' 
#' @param chromRange a \code{integer} \code{vector} enumerating chromosomes 
#' from which segments are to be used for initial model-based clustering.
#' Default: \code{NULL}.
#' 
#' @param nJobs a single positive \code{integer} specifying the number of 
#' worker jobs to create in case of distributed computation. 
#' Default: \code{1} and always \code{1} for Windows.
#' 
#' @param normalLength an integer \code{vector} specifying the genomic lengths 
#' of segments in the normal reference data. Default: \code{NULL}.
#' 
#' @param normalMedian a numeric \code{vector}, 
#' of the same length as \code{normalLength}, specifying the segment values
#' of the normal reference segments. Default: \code{NULL}.
#' 
#' @param normalMad a numeric \code{vector}, 
#' of the same length as \code{normalLength}, specifying the value spreads 
#' of the normal reference segments. Default: \code{NULL}.
#' 
#' @param normalError a numeric \code{vector}, 
#' of the same length as \code{normalLength}, specifying the error values
#' of the normal reference segments. Default: \code{NULL}.
#' 
#' @param weightall A \code{matrix} whose rows correspond to genomic positions 
#' and columns to copy number profiles (same as ratall). Its matrix elements 
#' are functions of weight of the genomic region.
#' Default: \code{NULL}.
#' 
#' @return If \code{keepCLust} is FALSE: the input \code{segall} \code{data.frame} to which some or all of 
#' the following columns may be bound, depending on the availability of input:
#' \itemize{
#' \item{segmedian}{ a \code{numeric}, the median function of copy number}
#' \item{segmad}{ a \code{numeric}, the MAD for the function of copy number}
#' \item{mediandev}{ a \code{numeric}, the median function of copy number 
#' relative to its central value}
#' \item{segerr}{ a \code{numeric}, the error estimate for the 
#' function of copy number}
#' \item{centerz}{ a \code{numeric} between \code{0} and \code{1}, the 
#' probability that the segment is in the central cluster}
#' \item{marginalprob}{ a \code{numeric}, the marginal probability for 
#' the segment in the central cluster}
#' \item{maxz}{ TODO}
#' \item{maxzmean}{ TODO}
#' \item{maxzsigma}{ TODO}
#' \item{samplesize}{ TODO}
#' \item{negtail}{ the probability of finding the deviation as observed or 
#' larger in a collection of central segments}
#' \item{negtailnormad}{ the probability of finding the deviation/MAD as 
#' observed or larger in a collection of central segments}
#' \item{negtailnormerror}{ a \code{numeric}, the probability of finding 
#' the deviation/error as observed or larger in a collection of 
#' central segments}
#' }
#' else if \code{keepCLust} is TRUE :
#' a \code{list} with the component
#' \itemize{
#' \item{maxz}{ segall define before}
#' \item{resCLust}{TODO}
#' }
#'
#' @details Depending on the availability of input, the function will 
#' perform the following operations for each copy number profile.
#' 
#' If raw data are available in addition to segment start and end positions,
#' median and MAD of each segment will be computed. For each profile, bootstrap 
#' sampling of the segment median values will be performed, and the sample will 
#' be used to estimate the error in the median for each segment. 
#' Model-dependent clustering (fitting to a gaussian mixture) of the sample 
#' will be performed. The central cluster (the one nearest the expected 
#' unaltered value) will be identified and, if necessary, merged with adjacent 
#' clusters in order to comprise the minimal required fraction of the data. 
#' Deviation of each segment from the center, its probability to belong to the 
#' central cluster and its marginal probability in the central cluster will be 
#' computed.
#' 
#' If segment medians or median deviations are available or have been computed, 
#' and, in addition, genomic lengths and average values are given for a 
#' collection of segments with unaltered copy number, additional estimates will 
#' be performed. If median values are available for the unaltered segments, the
#' marginal probability of the observed median or median deviation in the 
#' unaltered set will be computed for each segment. Likewise, marginal 
#' probabilities for median/MAD and/or median/error will be computed if these 
#' statistics are available. 
#' 
#' 
#' @examples
#'
#' ## Load needed datasets
#' data(segexample)
#' data(ratexample)
#' data(normsegs)
#' 
#' ## Small toy example
#' segtable <- CNpreprocessing(segall=segexample[segexample[,"ID"]=="WZ1",],
#'     ratall=ratexample, idCol="ID", startCol="start", endCol="end", 
#'     chromCol="chrom", bpStartCol="chrom.pos.start", 
#'     bpEndCol="chrom.pos.end", blsize=50, 
#'     minJoin=0.25, cWeight=0.4, bsTimes=50, chromRange=1:3, nJobs=1,
#'     modelNames="E", normalLength=normsegs[,1],
#'     normalMedian=normsegs[,2])
#'     
#' \dontrun{
#' ## Example 1: 5 whole genome analysis, choosing the right format of arguments
#' segtable <- CNpreprocessing(segall=segexample,ratall=ratexample, idCol="ID", 
#'    "start","end", chromCol="chrom",bpStartCol="chrom.pos.start",
#'    bpEndCol="chrom.pos.end", blsize=50, minJoin=0.25, cWeight=0.4, 
#'    bsTimes=50, chromRange=1:22, nJobs=4,
#'    modelNames="E", normalLength=normsegs[,1], normalMedian=normsegs[,2])
#'    
#' ## Example 2: how to use annotexample, when segment table does not have 
#' columns of integer positions in terms of  measuring units(probes), such as 
#' "mysegs" below
#' mysegs <- segexample[,c(1,5:12)]
#' 
#' data(annotexample)
#' 
#' segtable <- CNpreprocessing(segall=mysegs,ratall=ratexample, idCol="ID",
#'     chromCol="chrom", bpStartCol="chrom.pos.start",bpEndCol="chrom.pos.end",
#'     annot=annotexample, annotStartCol="CHROM.POS",annotEndCol="CHROM.POS",
#'     annotChromCol="CHROM", blsize=50, minJoin=0.25, cWeight=0.4, bsTimes=50,
#'     chromRange=1:22, modelNames="E", nJobs=4,
#'     normalLength=normsegs[,1], normalMedian=normsegs[,2])
#' }
#' 
#' @author Alexander Krasnitz
#' @importFrom BiocParallel multicoreWorkers SnowParam SerialParam bplapply bptry bpok
#' @export
CNpreprocessing <- function(segall, ratall=NULL, idCol=NULL, startCol=NULL,
    endCol=NULL, medCol=NULL, madCol=NULL, errorCol=NULL, chromCol=NULL,
    bpStartCol=NULL, bpEndCol=NULL, annot=NULL, annotStartCol=NULL, 
    annotEndCol=NULL, annotChromCol=NULL, useEnd=FALSE, blsize=NULL, 
    minJoin=NULL, nTrial=10, bestBIC=-1e7, modelNames="E", cWeight=NULL,
    bsTimes=NULL, chromRange=NULL, nJobs=1, normalLength=NULL, 
    normalMedian=NULL, normalMad=NULL,
    normalError=NULL, weightall=NULL, keepClust=FALSE) {
    
    configCall <- list(nTrial=nTrial, 
                       bestBIC=bestBIC, 
                       modelNames=modelNames, 
                       cWeight=cWeight,
                       bsTimes=bsTimes, 
                       chromRange=chromRange, 
                       nJobs=nJobs, 
                       normalLength=normalLength, 
                       keepClust=keepClust)

    ## Parameters validation
    validateCNpreprocessing(segall=segall, ratall=ratall, idCol=idCol, 
                            startCol=startCol,endCol=endCol, medCol=medCol, 
                            madCol=madCol, errorCol=errorCol, 
                            chromCol=chromCol, bpStartCol=bpStartCol, 
                            bpEndCol=bpEndCol, annot=annot, 
                            annotStartCol=annotStartCol, 
                            annotEndCol=annotEndCol, 
                            annotChromCol=annotChromCol, useEnd=useEnd, 
                            blsize=blsize, minJoin=minJoin, nTrial=nTrial, 
                            bestBIC=bestBIC, modelNames=modelNames, 
                            cWeight=cWeight, bsTimes=bsTimes, 
                            chromRange=chromRange, nJobs=nJobs,  
                            normalLength=normalLength, 
                            normalMedian=normalMedian, normalMad=normalMad,
                            normalError=normalError, weightall=NULL)

    ## Select the type of parallel environment used for parallel processing
    nbrThreads <- as.integer(nJobs)
    if (nbrThreads == 1 || multicoreWorkers() == 1) {
        coreParam <- SerialParam()
    } else {
        seed <- sample(x=seq_len(999999), size=1)
        coreParam <- SnowParam(workers = nbrThreads, RNGseed = seed)
    }
    
    # Initialise the variables for the final result
    res <- NULL
    resClust <- NULL
    
    ## When the column for the profile ID is not specified, see if it can
    ## deducted from the data
    ## If only one column in the ratall table is numerical, it will be 
    ## used as the profile ID column
    if (is.null(idCol)) {
        cat("Found a single segmented profile with no ID","\n")
        if (!is.null(ratall)) {
            if (sum(apply(ratall, 2, data.class) == "numeric") > 1) {
                stop(paste0("Ambiguity: more than 1 numeric column in ",
                                "\"ratall\" matrix. The \"idCol\" ", 
                                "parameter should be specified.\n"))
            } else {
                idrat <- which(apply(ratall, 2, data.class) == "numeric")
                segall <- data.frame(rep(as.character(idrat), nrow(segall)), 
                                        segall)
                idCol <- "ID"
                dimnames(segall)[[2]][1] <- idCol
                if (!is.null(weightall)) {
                    if(sum(apply(weightall, 2, data.class) == "numeric") > 1) {
                        stop("Ambiguity: more than 1 numeric column ", 
                             "in weight data table\n")
                    }
                }
            }
        } 
    }
    
    if (is.null(ratall)) {
        cat("No raw table, proceeding to comparison\n")
        if(!is.null(weightall)) {
            stop("weightall without ratall\n")
        }
        
    } else {
        profnames <- unique(segall[,idCol])
        
        if (!all(profnames %in% dimnames(ratall)[[2]])) {
            stop("Found unmatched segmented profile IDs\n")
        }
        
        if (!(is.null(weightall))) {
            if (!all(profnames%in%dimnames(weightall)[[2]])) {
                stop("Found unmatched segmented profile IDs in weightall\n")
            }
        }
        
        if (is.null(startCol) | is.null(endCol)) { 
            
            # Will need an annotation table
            
            if (is.null(bpStartCol) | is.null(bpEndCol) | is.null(chromCol)) {
                stop("Unable to proceed: incomplete segment annotation\n")
            }
            
            
            
            if (is.null(chromRange)) {
                chromRange <- sort(unique(segall[,chromCol]))
            }
            
            if (is.null(annot)) {
                stop(paste0("No annotation table; unable to determine ", 
                                "boundary probes/bins\n"))
            }
            
            if (is.null(annotStartCol) | is.null(annotChromCol)) {
                stop(paste0("No start and chrom column names provided for ", 
                                "annotation table\n"))
            }
            
            if (useEnd & is.null(annotEndCol)) {
                stop(paste0("End column name required but not provided in ", 
                            "annotation table\n"))
            }
            
            maxbpstart <- max(c(segall[,bpStartCol], annot[,annotStartCol])) + 1
            maxbpend <- ifelse(useEnd, 
                        max(c(segall[,bpEndCol], annot[,annotEndCol])),
                        max(c(segall[,bpEndCol], annot[,annotStartCol]))) + 1
            
            startprobe <- match((segall[,chromCol] - 1) * maxbpstart + 
                                            segall[,bpStartCol],
                            ceiling((annot[,annotChromCol] - 1) * maxbpstart + 
                                            annot[,annotStartCol]))
            
            endprobe <- ifelse(rep(useEnd, length(startprobe)),
                            match((segall[,chromCol] - 1) * maxbpend + 
                                        segall[,bpEndCol],
                                ceiling((annot[,annotChromCol] - 1) * maxbpend +
                                        annot[,annotEndCol])),
                            match((segall[,chromCol] - 1) * maxbpend + 
                                        segall[,bpEndCol],
                                ceiling((annot[,annotChromCol] - 1) * maxbpend +
                                        annot[,annotStartCol])))
            
            if (!all(!is.na(startprobe) & !is.na(endprobe))) {
                stop("Incomplete start and end annotation of segments\n")
            }
            
            segall <- data.frame(segall, startprobe, endprobe)
            dimnames(segall)[[2]][(ncol(segall)-1):ncol(segall)] <- 
                                                    c("StartProbe", "EndProbe")
            startCol <- "StartProbe"
            endCol <- "EndProbe"
        }
        
        profpack <- vector(mode="list", length=length(profnames))
        names(profpack) <- profnames
        
        for (pn in profnames) {
            profpack[[pn]] <- vector(mode="list", length=4)
            names(profpack[[pn]]) <- c("seg", "rat", "stream", "sub")
            profpack[[pn]]$seg <-
                segall[segall[,idCol] == pn, c(startCol, endCol, chromCol), 
                        drop=FALSE]
            dimnames(profpack[[pn]]$seg)[[2]] <- c("StartProbe", 
                                                    "EndProbe", "chrom")
            profpack[[pn]]$rat <- ratall[,pn]
            profpack[[pn]]$stream <- pn
            profpack[[pn]]$sub <- match(pn, profnames)
            
            if (!(is.null(weightall))) {
                profpack[[pn]]$weight <- weightall[, pn]
            } else {
                profpack[[pn]]$weight <- NULL
            }
        }
        
        rm(ratall)
        gc()

        ## Running each profile id on a separate thread
        processed <- bptry(bplapply(X=profpack, FUN=CNclusterNcenter, 
                                blsize=blsize, minJoin=minJoin, nTrial=nTrial, 
                                bestBIC=bestBIC, modelNames=modelNames, 
                                cweight=cWeight, bstimes=bsTimes, 
                                chromRange=chromRange, 
                                keepClust=keepClust,
                                BPPARAM=coreParam))
        ## Check for errors
        if (!all(bpok(processed))) {
           stop("At least one parallel task has thrown an error.")
        }

        #segall <- cbind(segall, do.call(rbind, processed))
        #saveRDS(processed, "processed.rds")
        if(keepClust){
            segall <- cbind(segall, 
                            do.call(rbind, 
                                    vapply(processed, 
                                           FUN=function(x){return(list(x$seg))},
                                           FUN.VALUE = list(1))))
            resClust <- vapply(processed, 
                               FUN=function(x){return(list(x$clustRes))},
                               FUN.VALUE = list(1))
            names(resClust) <- names(processed)
            
        } else{
            segall <- cbind(segall, do.call(rbind, processed))
        }
        

        dimnames(segall)[[2]][(ncol(segall)-8):ncol(segall)] <-
            c("segmedian", "segmad", "mediandev", "segerr", "centerz",
                "marginalprob", "maxz", "maxzmean", "maxzsigma")
        medCol <- "mediandev"
        madCol <- "segmad"
        errorCol <- "segerr"
        
    }
    
    ## Use normal samples to extract extra information when available
    if (!(is.null(normalLength) | is.null(normalMedian) | is.null(medCol))) {
        if (is.null(bpStartCol) | is.null(bpEndCol)) {  
            ## Try to annotate
            if (is.null(startCol) | is.null(endCol) | is.null(annot) | 
                is.null(annotStartCol) | is.null(annotEndCol)) {
                stop("Insufficient annotation for comparison")
            }
            
            tumorlength <- annot[segall[,endCol], annotEndCol] -
                annot[segall[,startCol], annotStartCol] + 1
        } else {
            tumorlength <- segall[, bpEndCol] - segall[, bpStartCol] + 1
        }
        
        tumormedian <- segall[, medCol]
        
        if (!is.null(madCol)) {
            tumormad <- segall[, madCol]
        }
        
        if (!is.null(errorCol)) {
            tumorerror <- segall[, errorCol]
        }
        
        segall <- cbind(segall, normalComparison(normalMedian, normalLength,
                    tumormedian, tumorlength, normalMad, normalError, tumormad, 
                    tumorerror))
    }
    
    if(keepClust){
        res <- list(segall=segall,
                    resClust=resClust)
    }else{
        res=segall
    }
    
    return(res)
}

