predict.lspcr.glm<-function(object,newX,newD,...){
  etaT=cbind(rep(1,dim(newD)[1]),newD,newX%*%object$projection)%*%object$coefficients
  piT=1/(1+exp(-etaT))
  YT=(etaT>0)+0
  return(list(newY=YT,newPi=piT)) 
}
