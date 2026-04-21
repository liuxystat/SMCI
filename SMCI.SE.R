SMCI.SE <- function(object, nboot,method=FALSE) {
  cores <- detectCores()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  dataset <- object$ParEst$sdata
  n<-nrow(dataset)
  tempdata <- dataset
  
  result<-foreach(inb = 1:nboot, .combine = rbind, .packages = c("GORCure","np","pROC","truncnorm","survival","ICsurv","mgcv","splines"),
                  .export = c('NZD','bb.eta','bb.eta0','ee.beta','ee.beta0','loglik_sp','EM.Ite_sp','SISMC2')) %dopar% {
                    
                    
                    Xbs<-list()
                    boots <- function(Xbs)
                    {
                      covg<-0
                      while(covg==0){
                        vec <- sample(1:n, n, replace = TRUE)
                        x<-tempdata[vec,]
                        if(method=='kernel'){
                          bootfit<-try(SMCI.K(survfun=object$survfun, curefun =object$curefun,data=data.frame(x),r=object$ParEst$r,n.int=object$ParEst$n.int,kertype=object$ParEst$kertype,max.iter = object$ParEst$max.iter,rescale = TRUE),silent=TRUE)
                        }else{
                          bootfit<-try(SMCI.S(survfun=object$survfun, curefun =object$curefun,data=data.frame(x),r=object$ParEst$r,n.int=object$ParEst$n.int,best_k=object$ParEst$best_k,degree = object$ParEst$degree,max.iter = object$ParEst$max.iter,rescale = TRUE),silent=TRUE)
                        }
                        
                        if (!is.character(bootfit) ) {
                          if(sum(abs(bootfit$ParEst$Beta))!=0){
                            est=c(bootfit$ParEst$Gamma,bootfit$ParEst$Beta)
                            covg<-1
                          }
                        }
                      }
                      return(est)
                    }
                    out<-boots(Xbs)
                    return(out)
                  }
  stopCluster(cl)
  varem <- sqrt(apply(result, 2, var))
  return(list(est.sd = varem,result=result))
}