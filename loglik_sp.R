loglik_sp <-function(x0,r,sdata,Zp,Xp,bl.Li,bl.Ri,pxi){
  M<-ncol(Xp)
  P<-ncol(Zp)
  g0<-x0[1:M]
  b0<-x0[(M+1):(M+P)]
  e0<-x0[-c(1:(M+P))]
  index <- Xp%*%g0
  ezb<-exp(Zp%*%b0)
  Lambda.Li<-t(bl.Li)%*%matrix(e0,ncol=1)
  Lambda.Ri<-t(bl.Ri)%*%matrix(e0,ncol=1)
  pc0<-(1-sdata$d3)*log(pxi)
  pc0[is.na(pc0)]<-0
  if(r>0){
    pc1<-sdata$d1*(1-(1+r*Lambda.Ri*ezb)^(-1/r))
    pc2<-sdata$d2*((1+r*Lambda.Li*ezb)^(-1/r)-(1+r*Lambda.Ri*ezb)^(-1/r))
    pc3<-sdata$d3*(1-pxi+pxi*(1+r*Lambda.Li*ezb)^(-1/r))
  }else{
    pc1<-sdata$d1*(1-exp(-Lambda.Ri*ezb))
    pc2<-sdata$d2*(exp(-Lambda.Li*ezb)-exp(-Lambda.Ri*ezb))
    pc3<-sdata$d3*(1-pxi+pxi*exp(-Lambda.Li*ezb))
  }
  ll<-sum(pc0+log(pc1+pc2+pc3))
  return(-ll)
}
