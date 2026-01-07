EM.Ite<-function(sdata,r,n.int,order,max.iter,cov.rate,kertype="epanechnikov",rescale = FALSE){
  Xp <- sdata$Xp
  Zp <- sdata$Zp
  P <- ncol(Zp)
  M <- ncol(Xp)
  L <- n.int+order
  N <- nrow(sdata)
  b0 <- rep(0,P)
  e0 <- rep(1,L)
  w <- 1-sdata$d3
  ga.data<-data.frame(Bi=w,Xp)
  tempf<-paste("x",1:M,sep="",collapse="+")
  formu<-paste("Bi~",tempf)
  #gamma.init <- glm(formu,family=binomial(link=logit),data=ga.data)$coefficients
  #g0 <- gamma.init[-1]
  gamma.init <- gam(w~Xp,family=binomial())$coefficients
  g0<-gamma.init[-1]
  CSI <- ifelse(g0[1]>0,1,-1)
  
  r.mat <- function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE) 
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  
  ti <- c(sdata$Li[sdata$d1 == 0], sdata$Ri[sdata$d3 == 0])
  ti.max <- max(ti) + 1e-05
  ti.min <- min(ti) - 1e-05
  #knots <- seq(ti.min, ti.max, length.out = (n.int + 2))
  knots <-quantile(ti,seq(0,1,length.out = (n.int + 2)))
  
  bl.Li<-bl.Ri<-matrix(0,nrow=L,ncol=N)
  bl.Li[,sdata$d1==0]<-Ispline(sdata$Li[sdata$d1==0],order=order,knots=knots)
  bl.Ri[,sdata$d3==0]<-Ispline(sdata$Ri[sdata$d3==0],order=order,knots=knots)
  
  
  dif<-numeric()
  est.par<-matrix(ncol=M+P+L,nrow=max.iter)
  ll<-numeric()
  dd<-1
  ii<-0
  h <- 0.2
  while(dd>cov.rate & ii<max.iter){
    ii<-ii+1
    
    index <- Xp%*%g0
    W <- as.matrix(data.frame(w,1))
    if(ii==1){bdwth <- 0.2}else{bdwth <- h}
    K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=bdwth,ckertype=kertype,leave.one.out=T)$ksum
    pxi <- K.sum[1,2,]/NZD(K.sum[2,2,])
    ezb<-exp(Zp%*%b0)
    Lambda.Li<-t(bl.Li)%*%matrix(e0,ncol=1)
    Lambda.Ri<-t(bl.Ri)%*%matrix(e0,ncol=1)
    
    if(r>0){
      lambdai<-r*ezb*(sdata$d1*Lambda.Ri+sdata$d2*Lambda.Li)
      omegai<-r*ezb*(sdata$d2*(Lambda.Ri-Lambda.Li)+sdata$d3*Lambda.Li)
      lambdail<-c.mat(e0)*r.mat(r*ezb)*(r.mat(sdata$d1)*bl.Ri+r.mat(sdata$d2)*bl.Li)
      omegail<-c.mat(e0)*r.mat(r*ezb)*(r.mat(sdata$d2)*(bl.Ri-bl.Li)+r.mat(sdata$d3)*bl.Li)
      
      ci<-1-(1/(1+sdata$d1*lambdai))^(1/r)
      ci<-ifelse(sdata$d1==1 & ci==0,sqrt(.Machine$double.eps),ci)
      di<-1-((1+sdata$d2*lambdai)/(1+sdata$d2*lambdai+sdata$d2*omegai))^(1/r)
      di<-ifelse(di==0 & sdata$d2==1,sqrt(.Machine$double.eps),di)
      ei<-1-pxi+pxi*(1/(1+sdata$d3*omegai))^(1/r)
      ei<-ifelse(ei==0 & sdata$d3==1,sqrt(.Machine$double.eps),ei)
      
      Eyil<-r.mat(sdata$d1/r/(sdata$d1*ci+1-sdata$d1))*lambdail
      Ewil<-r.mat(sdata$d2/r/(1+sdata$d2*lambdai)/(sdata$d2*di+1-sdata$d2))*omegail
      Eyi<-apply(Eyil,2,sum)
      Ewi<-apply(Ewil,2,sum)
      Eui<-1-sdata$d3+sdata$d3*pxi*(1/(1+sdata$d3*omegai))^(1/r)/ei
      Eui<-ifelse(sdata$d3 == 1&Eui==0, floor(Eui[Eui!=0]), Eui)
      
      pc1<-sdata$d1*(1-(1-ci)^(r+1))/r/(sdata$d1*ci+1-sdata$d1)
      pc2<-sdata$d2*(1-(1-di)^(r+1))/r/(1+sdata$d2*lambdai)/(sdata$d2*di+1-sdata$d2)
      pc3<-sdata$d3*(ei-(ei+pxi-1)*omegai/(1+omegai))/r/ei
      Efi<-pc1+pc2+pc3
      Euifi<-sdata$d3*pxi*(1/(1+omegai))^(1+1/r)/r/ei
    }else{
      lambdai<-ezb*(sdata$d1*Lambda.Ri+sdata$d2*Lambda.Li)
      omegai<-ezb*(sdata$d2*(Lambda.Ri-Lambda.Li)+sdata$d3*Lambda.Li)
      lambdail<-c.mat(e0)*(r.mat(sdata$d1*ezb)*bl.Ri+r.mat(sdata$d2*ezb)*bl.Li)
      omegail<-c.mat(e0)*(r.mat(sdata$d2*ezb)*(bl.Ri-bl.Li)+r.mat(sdata$d3*ezb)*bl.Li)
      
      ci<-1-exp(-sdata$d1*lambdai)
      ci<-ifelse(sdata$d1==1 & ci==0,sqrt(.Machine$double.eps),ci)
      di<-1-exp(-sdata$d2*omegai)
      di<-ifelse(di==0 & sdata$d2==1,sqrt(.Machine$double.eps),di)
      Eyil<-r.mat(sdata$d1/(sdata$d1*ci+1-sdata$d1))*lambdail
      Ewil<-r.mat(sdata$d2/(sdata$d2*di+1-sdata$d2))*omegail
      Eyi<-apply(Eyil,2,sum)      
      Ewi<-apply(Ewil,2,sum)
      
      Eui<-1-sdata$d3+sdata$d3*pxi*exp(-Lambda.Li*ezb)/(1-pxi+pxi*exp(-Lambda.Li*ezb))
      #Eui <- ifelse(sdata$d3 == 1&Eui==0, sqrt(.Machine$double.eps), Eui)
      Eui<-ifelse(sdata$d3 == 1&Eui==0, floor(Eui[Eui!=0]), Eui)
    }
    # M-step
    
    # Incidence
    w <- Eui
    B <- w
    CVh <- function(param){                                   
      h <- param[1]
      h.comp <- function(h){
        index <- Xp%*%(g0)
        W <- as.matrix(data.frame(w,1))
        K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=h,ckertype=kertype,leave.one.out=T)$ksum
        p2 <- K.sum[1,2,]/NZD(K.sum[2,2,])
        
        Floor1 <- sqrt(.Machine$double.eps) 
        #Floor2 <- 0.1
        p2[which(p2<Floor1)] <- Floor1 
        p2[which(p2>1-Floor1)] <- 1 - Floor1
        contrib <- rep(NA,N)
        for(i in 1:dim(Xp)[1]){
          if(B[i]==0) contrib[i] <- log(1-p2[i])
          if(B[i]>0 & B[i]<1) contrib[i] <- B[i]*log(p2[i])+(1-B[i])*log(1-p2[i])
          if(B[i]==1) contrib[i] <- log(p2[i])
        } 
        CV <- -sum(contrib)
        return(CV)  
      }
      if(h>0){return(h.comp(h))}else{return(sqrt(.Machine$double.xmax))}
    }
    
    if(ii==0){h.init <- 0.2}else{h.init <- h}
    h.CV <- optim(par=h.init,CVh,method='Brent',lower=0.1,upper=0.5,control=list(fnscale=1))
    h <- h.CV$par
    bdwth <- h 
    
    
    
    
    l_I <- function(gamma2){    
      Floor1 <- sqrt(.Machine$double.eps)
      Floor2 <- sqrt(.Machine$double.eps)
      theta <- c(CSI,gamma2)
      
      index <- Xp%*%theta
      W <- as.matrix(data.frame(B,1))
      K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=bdwth,ckertype=kertype,leave.one.out=T)$ksum
      p2 <- K.sum[1,2,]/NZD(K.sum[2,2,])
      p2[which(p2<Floor1)] <- Floor1
      p2[which(p2>1-Floor2)] <- 1 - Floor2
      contrib <- rep(NA,N)
      for(i in 1:N){
        if(B[i]==0) contrib[i] <- log(1-p2[i])
        if(B[i]>0 & B[i]<1) contrib[i] <- B[i]*log(p2[i])+(1-B[i])*log(1-p2[i])
        if(B[i]==1) contrib[i] <- log(p2[i])
      }  
      likelihood <- -sum(contrib)  
      
      return(likelihood)
    }
    
    
    starting<-g0[-1]
    
    
    fit1 <- optim(par = starting, fn = l_I, method = 'BFGS', control = list(fnscale = 1))
    
    #Latency
    if(r>0) fit2<-nlm(bb.eta,p=b0,r=r,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,Euifi=Euifi,bl.Li=bl.Li,bl.Ri=bl.Ri)
    if(r==0) fit2<-nlm(bb.eta0,p=b0,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Eui=Eui,bl.Li=bl.Li,bl.Ri=bl.Ri)
    
    if(!is.character(fit1) & !is.character(fit2)){
      g0[1]<-CSI
      g0[2:length(g0)]<-fit1$par
      
      L1<-fit1$value
      b0<-fit2$estimate
      
      if(r==0) e0<-ee.beta0(bb=b0,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Eui=Eui,bl.Li=bl.Li,bl.Ri=bl.Ri)
      if(r>0) e0<-ee.beta(bb=b0,r=r,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,Euifi=Euifi,bl.Li=bl.Li,bl.Ri=bl.Ri)
      
      est.par[ii,]<-c(g0,b0,e0)
      ll<-c(ll,loglik(x0=est.par[ii,],r=r,sdata=sdata,Zp=Zp,Xp=Xp,bl.Li=bl.Li,bl.Ri=bl.Ri,Eui=Eui,bdwth=bdwth,kertype=kertype))
    }
    if(ii>2) dd<-abs(ll[ii]-ll[ii-1])
    
  }
  # After convergence
  loglik_f<-ll[ii]
  AIC<-2*(L+P+M)+2*loglik_f
  if(rescale==TRUE){
    g1<-g0/sqrt(sum(g0^2))
    bdwth<-bdwth/sqrt(sum(g0^2))
  }else{
    g1<-g0
    bdwth<-bdwth
  }
  index <- Xp%*%g1
  W <- as.matrix(data.frame(Eui,1))
  K.sum <- npksum(txdat=index,tydat=W,weights=W,bws=bdwth,ckertype=kertype,leave.one.out=F)$ksum
  pxi <- K.sum[1,2,]/NZD(K.sum[2,2,])
  
  if(dd<cov.rate) convergence <-1 else convergence <- 0
  
  
  return(list(Gamma=g1,Beta=b0,Eta=e0,pre_pxi=pxi,bdwth=bdwth,Eui=Eui,AIC=AIC,loglik_f=-loglik_f,r=r,n.int=n.int,order=order,knots=knots,sdata=sdata,Xp=Xp,Zp=Zp,bl.Li=bl.Li,bl.Ri=bl.Ri,kertype=kertype,convergence=convergence,max.iter=max.iter))
}
