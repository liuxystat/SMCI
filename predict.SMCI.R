predict.SMCI <-function(object, method=FALSE,...){
  arg<-list(...)
  P<-length(object$ParEst$Beta)
  
  if(is.null(arg$len)) arg$len<-100
  if(is.null(arg$new.z)) arg$new.z<-rep(0,P)
  new.z<-as.vector(arg$new.z)
  
  mdata<-object$mdata
  
  ti<-unique(c(0,na.omit(mdata$Li),na.omit(mdata$Ri[is.finite(mdata$Ri)])))
  
  if(is.null(arg$tp)){
    tp<-seq(0,max(ti),length.out=arg$len)
  }else{
    tp<-arg$tp
  }
  
  if(is.null(arg$newdata)){
    tu<-seq(-2.5,3,length.out=arg$len)
  }else{
    tu<-as.matrix(arg$newdata)%*%object$ParEst$Gamma
  }
  
  
  if(method=='kernel'){
    index <- as.matrix(object$ParEst$sdata$Xp)%*%object$ParEst$Gamma
    index <- rbind(index,matrix(tu,ncol=1))
    w<-object$ParEst$Eui
    n1<-length(w)+1
    W1 <- as.matrix(data.frame(w,1))
    W2 <- matrix(0,ncol=2,nrow=length(tu))
    W <- rbind(W1,W2)
    n2<-nrow(W)
    bdwth <- object$ParEst$bdwth
    kertype<-object$ParEst$kertype
    K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=bdwth,ckertype=kertype,leave.one.out=F)$ksum
    pxi <- K.sum[1,2,]/NZD(K.sum[2,2,])
    pre_inci<-pxi[n1:n2]
  }
  if(method=='splines'){
    pre_inci<-predict(object$ParEst$incidence,newdata=data.frame(index=tu))
    Floor<-.Machine$double.eps
    pre_inci[which(pre_inci<Floor)] <- Floor
    pre_inci[which(pre_inci>(1-Floor))] <- 1-Floor
  }
  
  ezb<-exp(sum(object$ParEst$Beta*new.z))
  Het.est1<-t(Ispline(tp[tp<=max(ti)],order=object$ParEst$order,knots=object$ParEst$knots))%*%object$ParEst$Eta
  obj<-smooth.spline(tp[tp<=max(ti)],Het.est1)
  Het.est2<-predict(obj,x=tp[tp>max(ti)],deriv=0)$y
  
  if(object$ParEst$r>0)  sut<-(1+object$ParEst$r*c(Het.est1,Het.est2)*ezb)^(-1/object$ParEst$r)
  if(object$ParEst$r==0)  sut<-exp(-c(Het.est1,Het.est2)*ezb)
  
  pred<-list(Survival=data.frame(Time=tp,SurvProb=sut),pre_inci=data.frame(index=tu,pxi=pre_inci))
  class(pred)<-"SIM"
  class(pred)<-"predict.SIM"
  return(pred)
}
