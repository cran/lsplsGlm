predict.lspls<-function(fit,newX,newD){
  ZT=(newX-cbind(rep(1,dim(newD)[1]),newD)%*%fit$orthCoef)%*%fit$projection
  etaT=cbind(rep(1,dim(newD)[1]),newD,ZT)%*%fit$coefficients
  YT=(etaT>0)+0
  piT=1/(1+exp(-etaT))
  return (list(newY=YT,newPi=piT))
}

predict.rlspls<-function(fit,newX,newD){
  ZT=(newX-cbind(rep(1,dim(newD)[1]),newD)%*%fit$orthCoef)%*%fit$projection	
  etaT=cbind(newD,ZT)%*%fit$coefficients+fit$intercept
  YT=(etaT>0)+0	
  piT=1/(1+exp(-etaT))
  return (list(newY=YT,newPi=piT))
}

predict.irlspls<-function(fit,newX,newD){
  ZT=(newX-cbind(rep(1,dim(newD)[1]),newD)%*%fit$orthCoef)%*%fit$projection
  etaT=cbind(rep(1,dim(newD)[1]),newD,ZT)%*%fit$coefficients
  YT=(etaT>0)+0
  piT=1/(1+exp(-etaT))
  return(list(newY=YT,newPi=piT))
}

predict.lspls.glm<-function(object,newX,newD,...){
  
  if(class(object)=="LS-PLS-IRLS"){
    pred<-predict.lspls(object,newX=newX,newD=newD)  
    return(pred)
  }
  if (class(object)=="R-LS-PLS"){
    pred<-predict.rlspls(object,newX=newX,newD=newD)  
    return(pred)
  }
  if (class(object)=="IR-LS-PLS"){
    pred<-predict.irlspls(object,newX=newX,newD=newD)  
    return(pred)
  }
}
