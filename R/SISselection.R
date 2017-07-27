SIS.selection<-function(X,Y,pred,scale=F){
  if (scale==T){
    X<-scale(X)
  }
  p=dim(X)[2]
  IndicesSIS <- rep(0,p)
  beta <- rep(0,p)
  for (jj in 1:p){
    beta[jj]<-abs(glm(Y~X[,jj],family = binomial)$coefficients[2])
  }
  IndicesSIS <- sort(beta,index=TRUE,decreasing = TRUE)$ix		
  Xc<-X[,IndicesSIS[1:pred]] 
  Xsis <- Xc
  return(Xsis)
}