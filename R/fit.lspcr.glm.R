fit.lspcr.glm<-function(Y,X,D,ncomp){
  
  k=ncomp
  res <- prcomp(X,center=FALSE,scale=FALSE)
  Xred <- res$x[,1:k]
  V <- res$rotation[,1:k]
  if (ncomp==1)
  {
    Xred <- matrix(Xred,ncol=1)
    V <- matrix(V,ncol=1) ###  XLc%*%V-XLred &&&& t(V)%*%V=Id	
  }
  Xpcr = cbind(D,Xred)
  res=glm(Y~Xpcr,family = binomial)
  reponse<-list(coefficients=res$coefficients,cvg=res$converged,projection=V)
  class(reponse)<-"lspcr.glm"
  reponse
}
