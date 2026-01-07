bb.eta0 <-function(bb,sdata,Zp,Eyil,Ewil,Eui,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  ezb<-exp(Zp%*%bb)
  gg<-ee.beta0(bb,sdata,Zp,Eyil,Ewil,Eui,bl.Li,bl.Ri)
  prt1<-(c.mat(log(gg))+r.mat(log(ezb)))*(Eyil+Ewil)
  prt2<-c.mat(gg)*(r.mat(ezb*(1-sdata$d3))*bl.Ri+r.mat(ezb*sdata$d3*Eui)*bl.Li)
  Q<-sum((prt1-prt2))
  return(-Q)
}
