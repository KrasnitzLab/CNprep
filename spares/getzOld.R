getz <-
function(logr,emfit,zgroup,times){ 
	gz<-predict(emfit,newdata=logr)$z
	gz[is.nan(gz)]<-0	#just being honest: we don't know how to assign these
	gz<-gz%*%t(zgroup)
	gz<-matrix(ncol=ncol(gz),
		data=apply(gz,2,cumsum)[seq(from=times,to=nrow(gz),by=times),]/times)
	return(gz-rbind(matrix(nrow=1,data=rep(0,ncol(gz))),
		matrix(ncol=ncol(gz),data=gz[-nrow(gz),])))
}
