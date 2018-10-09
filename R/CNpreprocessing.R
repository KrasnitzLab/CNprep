#' @title Pre-process DNA copy number (CN) data for detection of CN events.
#' 
#' @description The package evaluates DNA copy number data, using both their 
#' initial form (copy number as a noisy function of genomic position) and their
#' approximation by a piecewise-constant function (segmentation), for the 
#' purpose of identifying genomic regions where the copy number differs from 
#' the norm.
#' 
#' @param segall A \code{matrix} or a \code{data.frame} for segmented copy 
#' number profiles. It may have a character column, with a name specified 
#' by \code{idcol}, and/or numeric columns with names specified by 
#' \code{startcol, endcol, medcol, madcol,errorcol}  
#' \code{,chromcol, bpstartcol, bpendcol}. Each row of \code{segall} 
#' corresponds to a segment belonging to one of the profiles 
#' to be pre-processed.
#' 
#' @param ratall A \code{matrix} whose rows correspond to genomic positions 
#' and columns to copy number profiles. Its matrix elements are functions of 
#' copy number, most often log ratios of copy number to the expected standard 
#' value, such as 2 in diploid genomes.
#' 
#' @param idcol A \code{character} string specifying the name for the 
#' column in \code{segall} tabulating the profile IDs.
#' 
#' @param startcol A \code{character} string specifying the name of column 
#' in \code{segall} that tabulates the (integer) start postion of each segment 
#' in internal units such as probe numbers for data of CGH microarray origin.
#' 
#' @param endcol A \code{character} string specifying the name of column 
#' in \code{segall} that tabulates the (integer) end postion of each segment 
#' in internal units such as probe numbers for data of CGH microarray origin.
#' 
#' @param medcol A \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' the function (\code{medcol}), a measure of its spread (\code{madcol}) and 
#' its error (\code{errorcol}) for the segment.
#' 
#' @param madcol A \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' a measure of spread (\code{madcol}) related to  
#' the function (\code{medcol}) for the segment.
#' 
#' @param errorcol A \code{character} string specifying the 
#' name of column in \code{segall} that, for the function of copy number used 
#' in the study (typically log ratios), tabulates the (numeric) values for 
#' the error (\code{errorcol}) related to  
#' the function (\code{medcol}) for the segment.
#' 
#' @param chromcol A character string specifying the name for the column in 
#' \code{segall} tabulating the (integer) chromosome number for each segment.
#' 
#' @param bpstartcol A character string specifying the name of 
#' column in \code{segall} that tabulates the (integer) genomic start 
#' coordinate of each segment.
#' 
#' @param bpendcol A character string specifying the name of 
#' column in \code{segall} that tabulates the (integer) genomic end 
#' coordinate of each segment.
#' 
#' @param annot A matrix or a data frame that contains the annotation for 
#' the copy number measurement platform in the study. It is generally expected 
#' to contain columns with names specified by 
#' \code{annotstartcol, annotendcol, annotchromcol}.
#' 
#' @param annotstartcol A character string 
#' specifying the name of column in \code{annot} that tabulates the (integer) 
#' genomic start coordinates in case of CGH
#' microarrays.
#' 
#' @param annotendcol A character string 
#' specifying the name of column in \code{annot} that tabulates the (integer) 
#' genomic end coordinates in case of CGH
#' microarrays.
#' 
#' @param annotchromcol A character string 
#' specifying the name of column in \code{annot} that tabulates the chromosome
#' number for each copy number measuring unit, such as a probe in case of CGH
#' microarrays.
#' 
#' @param useend A single logical value specifying whether the segment end 
#' positions as given by the \code{bpendcol} of \code{segall} are to be 
#' looked up in the \code{annotendcol} column of \code{annot} 
#' (if \code{useend=TRUE}) or in the \code{annotstartcol} column (default).
#' 
#' @param blsize A single \code{integer} specifying the bootstrap sampling 
#' rate of segment medians to generate input for model-based clustering. The 
#' number of times a segment is sampled is then given by the (integer) 
#' division of the segment length in internal units by \code{blsize}.
#' 
#' @param minjoin A single numeric value between 0 and 1 specifying the 
#' degree of overlap above which two clusters will be joined into one. 
#' 
#' @param ntrial A single integer specifying the number of times a model-based 
#' clustering is attempted for each profile in order to achieve the 
#' highest Bayesian information criterion (BIC).
#' 
#' @param bestbic A single \code{numeric} value for initalizing BIC 
#' maximization. A large negative value is recommended. The default 
#' is \code{-1e7}.
#' 
#' @param modelNames A \code{vector} of \code{character} strings specifying 
#' the names of models to be used in model-based clustering (see package 
#' \code{mclust} for further details). The default is \code{"E"}.
#' 
#' @param cweight A single \code{numeric} value between \code{0} and \code{1} 
#' specifying the minimal share of the central cluster in each profile.
#' 
#' @param bstimes A single \code{double} value specifying the number of 
#' time the median of each segment is sampled in order to predict the cluster 
#' assignment for the segment.
#' 
#' @param chromrange A \code{integer} vector enumerating chromosomes from 
#' which segments are to be used for initial model-based clustering.
#' 
#' @param myseed A single integer value to seed the random number generator.
#' 
#' @param distrib One of \code{"vanilla", "Rparallel"} to specify the 
#' distributed computing option for the cluster assignment step. 
#' For \code{"vanilla"} (default) no distributed computing is performed. For 
#' \code{"Rparallel"} the \code{parallel} package of \code{R} core is used 
#' for multi-core processing.
#' 
#' @param njobs A single integer specifying the number of worker jobs 
#' to create in case of distributed computation. Default: \code{1}.
#' 
#' @param normalength An integer vector specifying the genomic lengths of 
#' segments in the normal reference data.
#' 
#' @param normalmedian A numeric vector, 
#' of the same length as \code{normalength}, specifying the segment values
#' of the normal reference segments.
#' 
#' @param normalmad A numeric vector, 
#' of the same length as \code{normalength}, specifying the value spreads 
#' of the normal reference segments.
#' 
#' @param normalerror A numeric vector, 
#' of the same length as \code{normalength}, specifying the error values
#' of the normal reference segments.
#' 
#' @return The input \code{segall} data frame to which some or all of 
#' the following columns may be bound, depending on the availability of input:
#' \itemize{
#' \item{segmedian}{Median function of copy number}
#' \item{segmad}{MAD for the function of copy number}
#' \item{mediandev}{median function of copy number relative to its 
#' central value}
#' \item{segerr}{error estimate for the function of copy number}
#' \item{segz}{the probability that the segment is in the central cluster}
#' \item{marginalprob}{marginal probability for the segment in the central 
#' cluster}
#' \item{negtail}{the probability of finding the deviation as observed or larger 
#' in a collection of central segments}
#' \item{negtailnormad}{the probability of finding the deviation/MAD as observed 
#' or larger in a collection of central segments}
#' \item{negtailnormerror}{the probability of finding the deviation/error as 
#' observed or larger in a collection of central segments}
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
#' Deviation of each segment from the center, its probabilty to belong to the 
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
#' data(segexample)
#' data(ratexample)
#' data(normsegs)
#' 
#' ## Small toy example
#' segtable<-CNpreprocessing(segall=segexample[segexample[,"ID"]=="WZ1",],
#'     ratall=ratexample, "ID", "start", "end", chromcol="chrom", 
#'     bpstartcol="chrom.pos.start", bpendcol="chrom.pos.end", blsize=50, 
#'     minjoin=0.25, cweight=0.4, bstimes=50, chromrange=1:3, 
#'     distrib="vanilla", njobs=1, modelNames="E", normalength=normsegs[,1],
#'     normalmedian=normsegs[,2])
#'     
#' \dontrun{
#' ## Example 1: 5 whole genome analysis, choosing the right format of arguments
#' segtable<-CNpreprocessing(segall=segexample,ratall=ratexample, "ID", 
#'    "start","end", chromcol="chrom",bpstartcol="chrom.pos.start",
#'    bpendcol="chrom.pos.end",blsize=50, minjoin=0.25,cweight=0.4,bstimes=50,
#'    chromrange=1:22,distrib="Rparallel",njobs=40, modelNames="E",
#'    normalength=normsegs[,1],normalmedian=normsegs[,2])
#'    
#' ## Example 2: how to use annotexample, when segment table does not have 
#' columns of integer postions in terms of  measuring units(probes), such as 
#' "mysegs" below
#' mysegs<-segexample[,c(1,5:12)]
#' data(annotexample)
#' segtable<-CNpreprocessing(segall=mysegs,ratall=ratexample,"ID",
#'     chromcol="chrom", bpstartcol="chrom.pos.start",bpendcol="chrom.pos.end",
#'     annot=annotexample, annotstartcol="CHROM.POS",annotendcol="CHROM.POS",
#'     annotchromcol="CHROM", blsize=50,minjoin=0.25,cweight=0.4,bstimes=50,
#'     chromrange=1:22,distrib="Rparallel", njobs=40,modelNames="E", 
#'     normalength=normsegs[,1],normalmedian=normsegs[,2])
#' }
#' 
#' @author Alexander Krasnitz
#' @importFrom parallel makeCluster clusterEvalQ parLapply stopCluster detectCores
#' @export
CNpreprocessing <- function(segall,ratall=NULL,idcol=NULL,startcol=NULL,
	endcol=NULL,medcol=NULL,madcol=NULL,errorcol=NULL,chromcol=NULL,
	bpstartcol=NULL,bpendcol=NULL,annot=NULL,annotstartcol=NULL,annotendcol=NULL,
	annotchromcol=NULL,useend=F,blsize=NULL,minjoin=NULL,ntrial=10,bestbic=-1e7,
	modelNames="E",cweight=NULL,bstimes=NULL,
	chromrange=NULL,myseed=123,distrib=c("vanilla","Rparallel"),njobs=1,
	normalength=NULL,normalmedian=NULL,normalmad=NULL,normalerror=NULL) {
    
	#try to see what's possible with this input
	if(is.null(idcol)){
		cat("Found a single segmented profile with no ID","\n")
		if(!is.null(ratall)){
			if(sum(apply(ratall,2,data.class)=="numeric")>1)
				stop("Ambiguity: more than 1 numeric column in raw data table\n")
			else{ 
				idrat<-which(apply(ratall,2,data.class)=="numeric")
				segall<-data.frame(rep(as.character(idrat),nrow(segall)),segall)
				idcol<-"ID"
				dimnames(segall)[[2]][1]<-idcol
			}
		}
	}
	if(is.null(ratall))cat("No raw table, proceeding to comparison\n")
	else{
		profnames<-unique(segall[,idcol])	
		if(!all(profnames%in%dimnames(ratall)[[2]]))
			stop("Found unmatched segmented profile IDs\n")
		if(is.null(startcol)|is.null(endcol)){	#will need an annotation table
			if(is.null(bpstartcol)|is.null(bpendcol)|is.null(chromcol))
				stop("Unable to proceed: incomplete segment annotation\n")
			if(is.null(chromrange))chromrange<-sort(unique(segall[,chromcol]))
			if(is.null(annot))
				stop("No annotation table; unable to determine boundary probes/bins\n")
			if(is.null(annotstartcol)|is.null(annotchromcol))
				stop("No start and chrom column names provided for annotation table\n")
			if(useend&is.null(annotendcol))
				stop("End column name required but not provided in annotation table\n")
			maxbpstart<-max(c(segall[,bpstartcol],annot[,annotstartcol]))+1
			maxbpend<-ifelse(useend,max(c(segall[,bpendcol],annot[,annotendcol])),
				max(c(segall[,bpendcol],annot[,annotstartcol])))+1
			startprobe<-match((segall[,chromcol]-1)*maxbpstart+segall[,bpstartcol],
				ceiling((annot[,annotchromcol]-1)*maxbpstart+annot[,annotstartcol]))
			endprobe<-ifelse(rep(useend,length(startprobe)),
				match((segall[,chromcol]-1)*maxbpend+segall[,bpendcol],
				ceiling((annot[,annotchromcol]-1)*maxbpend+annot[,annotendcol])),
				match((segall[,chromcol]-1)*maxbpend+segall[,bpendcol],
				ceiling((annot[,annotchromcol]-1)*maxbpend+annot[,annotstartcol])))
			if(!all(!is.na(startprobe)&!is.na(endprobe)))
				stop("Incomplete start and end annotation of segments\n")
			segall<-data.frame(segall,startprobe,endprobe)
			dimnames(segall)[[2]][(ncol(segall)-1):ncol(segall)]<-c("StartProbe","EndProbe")
			startcol<-"StartProbe"
			endcol<-"EndProbe"
		}
		profpack<-vector(mode="list",length=length(profnames))
		names(profpack)<-profnames
		for(pn in profnames){
			profpack[[pn]]<-vector(mode="list",length=4)
			names(profpack[[pn]])<-c("seg","rat","stream","sub")
			profpack[[pn]]$seg<-
				segall[segall[,idcol]==pn,c(startcol,endcol,chromcol),drop=F]
			dimnames(profpack[[pn]]$seg)[[2]]<-c("StartProbe","EndProbe","chrom")
			profpack[[pn]]$rat<-ratall[,pn]
			profpack[[pn]]$stream<-pn
			profpack[[pn]]$sub<-match(pn,profnames)
		}
		rm(ratall)
		gc()
		distrib<-match.arg(distrib)
		if(distrib=="Rparallel"){
			ncores<-min(njobs,length(profnames),detectCores())
			cl<-parallel::makeCluster(getOption("cl.cores",ncores))
			parallel::clusterEvalQ(cl=cl,expr=requireNamespace("rlecuyer"))
			parallel::clusterEvalQ(cl=cl,expr=requireNamespace("mclust"))
			parallel::clusterEvalQ(cl=cl,expr=requireNamespace("CNprep"))
		}
		processed<-switch(distrib,
			vanilla=lapply(X=profpack,FUN=CNclusterNcenter,blsize=blsize,
				minjoin=minjoin,ntrial=ntrial,bestbic=bestbic,modelNames=modelNames,
				cweight=cweight,bstimes=bstimes,chromrange=chromrange,seedme=myseed),
			Rparallel=parLapply(cl,X=profpack,fun=CNclusterNcenter,blsize=blsize,
				minjoin=minjoin,ntrial=ntrial,bestbic=bestbic,modelNames=modelNames,
				cweight=cweight,bstimes=bstimes,chromrange=chromrange,seedme=myseed))
		if(distrib=="Rparallel")stopCluster(cl)
		segall<-cbind(segall,do.call(rbind,processed))
		dimnames(segall)[[2]][(ncol(segall)-8):ncol(segall)]<-
			c("segmedian","segmad","mediandev","segerr","centerz","marginalprob",
			"maxz","maxzmean","maxzsigma")
		medcol<-"mediandev"
		madcol<-"segmad"
		errorcol<-"segerr"
	}
	if(!(is.null(normalength)|is.null(normalmedian)|is.null(medcol))){
		if(is.null(bpstartcol)|is.null(bpendcol)){	#try to annotate
			if(is.null(startcol)|is.null(endcol)|is.null(annot)|is.null(annotstartcol)
				|is.null(annotendcol))stop("Insufficient annotation for comparison")
			tumorlength<-annot[segall[,endcol],annotendcol]-
				annot[segall[,startcol],annotstartcol]+1
		}
		else tumorlength<-segall[,bpendcol]-segall[,bpstartcol]+1
		tumormedian<-segall[,medcol]
		if(!is.null(madcol))tumormad<-segall[,madcol]
		if(!is.null(errorcol))tumorerror<-segall[,errorcol]
		segall<-cbind(segall,normalComparison(normalmedian,normalength,
  		tumormedian,tumorlength,normalmad,normalerror,tumormad,tumorerror))
	}
	return(segall)
}
