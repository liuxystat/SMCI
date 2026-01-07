SMCI.S <-function(survfun=formula(data),curefun=formula(data),data=parent.frame(),
                  r=0,n.int=5,order=3,best_k=NULL,max.iter=1000,cov.rate=.001,rescale=FALSE){
  call <- match.call()
  temp <- c("", "formula", "data", "na.action")
  avars<-all.vars(survfun)
  
  sdata <- data
  sdata$d1 <- ifelse(sdata[,avars[1]]==0 | is.na(sdata[,avars[1]]), 1, 0)
  sdata$d3 <- ifelse(is.infinite(sdata[,avars[2]]) | is.na(sdata[,avars[2]]), 1, 0)
  sdata$d2 <- 1-sdata$d1-sdata$d3
  sdata$Li <- sdata[,avars[1]]
  sdata$Ri <- sdata[,avars[2]]
  
  survfun1<-as.formula(paste("Surv(",avars[1],",",avars[2],",type='interval2')~",paste(avars[-c(1,2)],collapse="+"),sep=""))
  temp.z<-terms(survfun1, data = sdata)
  temp.x<-terms(curefun, data = sdata)
  attr(temp.x, "intercept") <- 0
  attr(temp.z, "intercept") <- 0
  Xp <- model.matrix(temp.x,sdata)
  colnames(Xp)<-paste("x",1:ncol(Xp),sep="")
  Zp <- model.matrix(temp.z,sdata)
  colnames(Zp)<-paste("z",1:ncol(Zp),sep="")
  
  L<-n.int+order
  P<-ncol(Zp)
  M<-ncol(Xp)
  N<-nrow(sdata)
  sdata$Xp <- Xp
  sdata$Zp <- Zp
  # Estimates
  temp.est<-EM.Ite_sp(sdata=sdata,r=r,n.int=n.int,order=order,best_k=best_k,max.iter=max.iter,cov.rate=cov.rate,rescale=rescale)
  output<-list(survfun=survfun1,curefun=curefun,ParEst=temp.est,call=match.call())
  class(output)<-"SISMC"
  output$mdata<-sdata
  return(output)
}
