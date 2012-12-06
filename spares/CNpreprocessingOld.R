CNpreprocessing <-
function(segall,ratall=NULL,idcol=NULL,startcol=NULL,
	endcol=NULL,medcol=NULL,madcol=NULL,errorcol=NULL,chromcol=NULL,
	bpstartcol=NULL,bpendcol=NULL,annot=NULL,annotstartcol=NULL,annotendcol=NULL,
	annotchromcol=NULL,useend=F,blsize=NULL,minjoin=NULL,ntrial=10,bestbic=-1e7,
	modelNames="E",cweight=NULL,bstimes=NULL,
	chromrange=NULL,myseed=123,distrib=c("vanilla","Rparallel","sge"),njobs=1,
	normalength=NULL,normalmedian=NULL,normalmad=NULL,normalerror=NULL){
	if(length(grep("-v PATH",getOption("sge.user.options")))==0)
		options(sge.user.options=paste(getOption("sge.user.options"),"-v PATH"))
	#sge.options(sge.save.global=TRUE)
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
				(annot[,annotchromcol]-1)*maxbpstart+annot[,annotstartcol])
			endprobe<-ifelse(useend,
				match((segall[,chromcol]-1)*maxbpend+segall[,bpendcol],
				(annot[,annotchromcol]-1)*maxbpend+annot[,annotendcol]),
				match((segall[,chromcol]-1)*maxbpend+segall[,bpendcol],
				(annot[,annotchromcol]-1)*maxbpend+annot[,annotstartcol]))
			if(!all(!is.na(startprobe)&!is.na(endprobe)))
				stop("Incomplete start and end annotation of segments\n")
			segall<-data.frame(segall,startprobe,endprobe)
			startcol<-"StartProbe"
			endcol<-"EndProbe"
		}
		profpack<-vector(mode="list",length=length(profnames))
		names(profpack)<-profnames
		distrib<-match.arg(distrib)
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
		distrib<-match.arg(distrib)
		if(distrib=="Rparallel"){
			ncores<-min(njobs,length(profnames),detectCores())
			cl<-parallel::makeCluster(getOption("cl.cores",ncores))
			parallel::clusterEvalQ(cl=cl,expr=library(rlecuyer))
			parallel::clusterEvalQ(cl=cl,expr=library(mclust))
			parallel::clusterEvalQ(cl=cl,expr=library(CNprep))
		}
		if(distrib=="sge")ncores<-min(njobs,length(profnames))
		processed<-switch(distrib,
			vanilla=lapply(X=profpack,FUN=CNclusterNcenter,blsize=blsize,
				minjoin=minjoin,ntrial=ntrial,bestbic=bestbic,modelNames=modelNames,
				cweight=cweight,bstimes=bstimes,chromrange=chromrange,seedme=myseed),
			Rparallel=parLapply(cl,X=profpack,fun=CNclusterNcenter,blsize=blsize,
				minjoin=minjoin,ntrial=ntrial,bestbic=bestbic,modelNames=modelNames,
				cweight=cweight,bstimes=bstimes,chromrange=chromrange,seedme=myseed),
			sge=sge.parLapply(X=profpack,FUN=CNclusterNcenter,blsize=blsize,
				minjoin=minjoin,ntrial=ntrial,bestbic=bestbic,modelNames=modelNames,
				cweight=cweight,bstimes=bstimes,chromrange=chromrange,seedme=myseed,
				packages=c("rlecuyer","mclust","CNprep"),njobs=ncores))
		if(distrib=="Rparallel")stopCluster(cl)
		segall<-cbind(segall,t(matrix(nrow=6,data=unlist(lapply(processed,t)))))
		dimnames(segall)[[2]][(ncol(segall)-5):ncol(segall)]<-
			c("segmedian","segmad","mediandev","segerr","segz","marginalprob")
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
