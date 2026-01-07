bb.eta <-function(bb,r,sdata,Zp,Eyil,Ewil,Efi,Euifi,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  ezb<-exp(Zp%*%bb)
  gg<-ee.beta(bb,r,sdata,Zp,Eyil,Ewil,Efi,Euifi,bl.Li,bl.Ri)
  prt1<-(Eyil+Ewil)*(c.mat(log(gg))+r.mat(log(ezb)))
  prt2<-r*c.mat(gg)*(r.mat(ezb*Efi*(1-sdata$d3))*bl.Ri+r.mat(ezb*Euifi*sdata$d3)*bl.Li)
  Q<-sum((prt1-prt2))
  return(-Q)
}