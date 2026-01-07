ee.beta0 <-function(bb,sdata,Zp,Eyil,Ewil,Eui,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  ezb<-exp(Zp%*%bb)
  numer<-apply(Eyil+Ewil,1,sum)
  denom<-apply(r.mat((1-sdata$d3)*ezb)*bl.Ri+r.mat(sdata$d3*ezb*Eui)*bl.Li,1,sum)
  gg<-numer/denom
  return(gg)
}
