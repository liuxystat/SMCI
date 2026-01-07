ee.beta <-function(bb,r,sdata,Zp,Eyil,Ewil,Efi,Euifi,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  ezb<-exp(Zp%*%bb)
  numer<-apply(Eyil+Ewil,1,sum)
  denom<-apply(r.mat(r*ezb*(1-sdata$d3)*Efi)*bl.Ri+r.mat(r*ezb*sdata$d3*Euifi)*bl.Li,1,sum)
  gg<-numer/denom
  return(gg)
}
