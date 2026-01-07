EM.Ite_sp<-function(sdata,r,n.int,order,max.iter,cov.rate,rescale,best_k=NULL,m_order=3,knotx,fam0=gaussian){
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
  gamma.init <- gam(w ~ Xp, family=fam0)$coefficients
  g0 <- gamma.init[-1]
  CSI <- ifelse(g0[1]>0,1,-1)
  g0 <- g0/abs(g0[1])
  
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
  #h <- 0.5
  while(dd>cov.rate & ii<max.iter){
    ii<-ii+1
    
    index <- Xp%*%g0
    windex <-data.frame(w=w,index=index)
    bic_df <- function(data){
      bic_values <- vector()
      for (k in 3:15) {
        model<-lm(w~bs(index,df=k),data = data)
        bic_values[k-2] <- BIC(model)
      }
      best_k <- which.min(bic_values)+2
      return(best_k)
    }
    
    if(is.null(best_k)){
      best_k <- bic_df(data=windex)
    }else{best_k <- best_k}
    pxi <- lm(w~bs(index,df=best_k),data = windex)$fitted.values
    Floor<-.Machine$double.eps
    pxi[which(pxi<Floor)] <- Floor
    pxi[which(pxi>(1-Floor))] <- 1-Floor
    
    ezb<-exp(Zp%*%b0)
    Lambda.Li<-t(bl.Li)%*%matrix(e0,ncol=1)
    Lambda.Ri<-t(bl.Ri)%*%matrix(e0,ncol=1)
    
    if(r>0){
      lambdai<-r*ezb*(sdata$d1*Lambda.Ri+sdata$d2*Lambda.Li)
      omegai<-r*ezb*(sdata$d2*(Lambda.Ri-Lambda.Li)+sdata$d3*Lambda.Li)
      lambdail<-c.mat(e0)*r.mat(r*ezb)*(r.mat(sdata$d1)*bl.Ri+r.mat(sdata$d2)*bl.Li)
      omegail<-c.mat(e0)*r.mat(r*ezb)*(r.mat(sdata$d2)*(bl.Ri-bl.Li)+r.mat(sdata$d3)*bl.Li)
      
      ci<-1-(1/(1+sdata$d1*lambdai))^(1/r)
      di<-1-((1+sdata$d2*lambdai)/(1+sdata$d2*lambdai+sdata$d2*omegai))^(1/r)
      ei<-1-pxi+pxi*(1/(1+sdata$d3*omegai))^(1/r)
      
      Eyil<-r.mat(sdata$d1/r/(sdata$d1*ci+1-sdata$d1))*lambdail
      Ewil<-r.mat(sdata$d2/r/(1+sdata$d2*lambdai)/(sdata$d2*di+1-sdata$d2))*omegail
      Eyi<-apply(Eyil,2,sum)
      Ewi<-apply(Ewil,2,sum)
      Eui<-1-sdata$d3+sdata$d3*pxi*(1/(1+sdata$d3*omegai))^(1/r)/ei
      
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
      di<-1-exp(-sdata$d2*omegai)
      Eyil<-r.mat(sdata$d1/(sdata$d1*ci+1-sdata$d1))*lambdail
      Ewil<-r.mat(sdata$d2/(sdata$d2*di+1-sdata$d2))*omegail
      Eyi<-apply(Eyil,2,sum)      
      Ewi<-apply(Ewil,2,sum)
      
      Eui<-1-sdata$d3+sdata$d3*pxi*exp(-Lambda.Li*ezb)/(1-pxi+pxi*exp(-Lambda.Li*ezb))
    }
    # M-step
    
    # Incidence
    l_I <- function(zeta0, opt=TRUE) {
      theta <- c(CSI,zeta0) 
      if (rescale==TRUE) {
        theta <- theta/sqrt(sum(theta^2)) 
      }
      a <- Xp%*%theta 
      b <- lm(Eui~bs(a,df=best_k))
      pxi<-b$fitted.values
      pxi[which(pxi<Floor)] <- Floor
      pxi[which(pxi>(1-Floor))] <- 1-Floor
      Q1 <- sum(Eui*log(pxi)+(1-Eui)*log(1-pxi))
      if(opt){
        return(-Q1)
      }else{
        return(list(theta=theta,Q1=Q1,pxi=pxi))
      }
      
    }
    
    starting<-g0[-1]/abs(g0[1])
    
    fit1 <- nlm(l_I,p=starting)
    
    
    #Latency
    if(r>0) fit2<-nlm(bb.eta,p=b0,r=r,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,Euifi=Euifi,bl.Li=bl.Li,bl.Ri=bl.Ri)
    if(r==0) fit2<-nlm(bb.eta0,p=b0,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Eui=Eui,bl.Li=bl.Li,bl.Ri=bl.Ri)
    
    if(!is.character(fit1) & !is.character(fit2)){
      g0 <- l_I(fit1$estimate,opt = FALSE)$theta
      pxi <- l_I(fit1$estimate,opt = FALSE)$pxi
      
      
      L1<--fit1$minimum
      b0<-fit2$estimate
      #b0<-fit2$par
      if(r==0) e0<-ee.beta0(bb=b0,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Eui=Eui,bl.Li=bl.Li,bl.Ri=bl.Ri)
      if(r>0) e0<-ee.beta(bb=b0,r=r,sdata=sdata,Zp=Zp,Eyil=Eyil,Ewil=Ewil,Efi=Efi,Euifi=Euifi,bl.Li=bl.Li,bl.Ri=bl.Ri)
      est.par[ii,]<-c(g0,b0,e0)
      ll<-c(ll,loglik_sp(x0=est.par[ii,],r=r,sdata=sdata,Zp=Zp,Xp=Xp,bl.Li=bl.Li,bl.Ri=bl.Ri,pxi=pxi))
    }
    if(ii>2) dd<-abs(ll[ii]-ll[ii-1])
    
  }
  # After convergence
  if(dd<cov.rate) convergence <-1 else convergence <- 0
  loglik_f<-ll[ii]
  AIC<-2*(L+P+M)+2*loglik_f
  
  index <- Xp%*%g0
  
  incidence <- lm(Eui~bs(index,df=best_k))
  return(list(Gamma=g0,Beta=b0,Eta=e0,pre_pxi=pxi,Eui=Eui,AIC=AIC,loglik_f=-loglik_f,r=r,n.int=n.int,order=order,knots=knots,best_k=best_k,sdata=sdata,Xp=Xp,Xp=Xp,bl.Li=bl.Li,bl.Ri=bl.Ri,best_k=best_k,incidence=incidence,convergence=convergence,max.iter=max.iter))
}
