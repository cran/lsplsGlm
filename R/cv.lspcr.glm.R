cv.lspcr.glm<-function(Y,X,D,ncompmax,folds=5,proportion=0.9){
  ncompopt<-1
  p.cvg<-rep(0,ncompmax)
  if (ncompmax!=1){
    Erreur<-rep(0,ncompmax)
    for (t in 1:folds){ # 5 tirages de 90 /10
      set.seed(200) 
      ind<-preselected.sample(Y,trunc(proportion*length(Y))) 
      cvXL<-X[ind,]
      cvDL<-D[ind,]
      cvYL<-Y[ind]
      cvDT<-D[-ind,]
      cvXT<-X[-ind,]
      cvYT<-Y[-ind]	
      
      for (k in (1:ncompmax)){
        res<-fit.lspcr.glm(Y=cvYL,X=cvXL,D=cvDL,ncomp=k)
        pred<-predict.lspcr.glm(res,newD=cvDT,newX=cvXT)
        Erreur[k]<-Erreur[k]+sum(abs(pred$newY-cvYT))
        p.cvg[k]<-p.cvg[k]+res$cvg
      }
    }	
    ncompopt<-which.min(Erreur)
  }
  return(list(ncompopt=ncompopt,p.cvg=p.cvg/folds))
}
