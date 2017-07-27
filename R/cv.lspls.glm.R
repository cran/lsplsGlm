cv.lspls<-function(Y,X,D,ncompmax,folds=5,proportion=0.9){
  
  ncompopt<-1
  p.cvg<-c()
  if (ncompmax!=1){
    Erreur<-rep(0,ncompmax)
    for (t in 1:folds){
      set.seed(200)
      ind<-preselected.sample(Y,trunc(proportion*length(Y))) 
      cvXL<-X[ind,]
      cvDL<-D[ind,]
      cvYL<-Y[ind]
      cvXT<-X[-ind,]
      cvDT<-D[-ind,]
      cvYT<-Y[-ind]	
      res <- lspls(Y=cvYL, D=cvDL, X=cvXL, W=diag(rep(1,nrow(cvDL))), ncomp=ncompmax)
      XM2 <- res$predictors
      projection <- res$projection
      orthCoef <- res$orthCoef
      for (k in (1:ncompmax))
      {	
        res=glm(cvYL~XM2[,1:(ncol(cvDL)+k)],family = binomial)	
        ZT=(cvXT-cbind(rep(1,dim(cvDT)[1]),cvDT)%*%orthCoef)%*%projection[,1:k]
        etaT=cbind(rep(1,dim(cvDT)[1]),cvDT,ZT)%*%res$coefficients
        cvYTprev<-(etaT>0)+0
        Erreur[k]<-Erreur[k]+sum(abs(cvYTprev-cvYT))
        p.cvg[k]<-p.cvg[k]+res$cvg
      }
    }	
    ncompopt<-which.min(Erreur)
  }
  return(list(ncompopt=ncompopt,p.cvg=p.cvg/folds))
}

cv.rlspls<-function(Y,X,D,ncompmax,folds=5,proportion=0.9,lambda.grid,penalized,nbrIterMax,threshold){
  lambda.grid <- sort(lambda.grid)
  Erreur<-matrix(0,nrow=ncompmax,ncol=length(lambda.grid))
  Convergence<-rep(1,length(lambda.grid)) #modif le 16/7/15
  p.cvg<-matrix(0,nrow=ncompmax,ncol=length(lambda.grid))
  for (t in 1:folds){
    set.seed(200)
    ind<-preselected.sample(Y,trunc(proportion*length(Y))) 
    cvXL<-X[ind,]
    cvDL<-D[ind,]
    cvYL<-Y[ind]
    cvXT<-X[-ind,]
    cvDT<-D[-ind,]
    cvYT<-Y[-ind]	
    for (l in 1:length(lambda.grid))
    {lam <- lambda.grid[l]
    res <- rirls(Y=cvYL,D=cvDL,X=cvXL,penalized=penalized,lambda=lam,nbrIterMax=nbrIterMax,threshold=threshold) 
    if (res$cvg == 0){Convergence[l]<-0} #modif le 16/7/15
    W <- res$W
    Ylspls <- res$pseudoresponse
    Dlspls <- cvDL
    Xlspls <- cvXL
    for (k in (1:ncompmax))
    {reslspls <- lspls(Y=Ylspls, D=Dlspls, X=Xlspls, W=W, ncomp=k)
    projection <- reslspls$projection
    orthCoef <- reslspls$orthCoef
    ZT=(cvXT-cbind(rep(1,dim(cvDT)[1]),cvDT)%*%orthCoef)%*%projection
    etaT=cbind(cvDT,ZT)%*%reslspls$coefficients+reslspls$intercept
    cvYTprev=(etaT>0)+0	
    Erreur[k,l]<-Erreur[k,l]+sum(abs(cvYTprev-cvYT))	
    p.cvg[k,l]<-res$cvg+p.cvg[k,l]
    }
    }
  }
  
  L=t(matrix(rep(lambda.grid,ncompmax),nrow=length(lambda.grid))) #modif le 16/7/15
  C=matrix(rep(1:ncompmax,length(lambda.grid)),nrow=ncompmax) #modif le 16/7/15
  Erreur <- Erreur[,Convergence==1] #modif le 16/7/15
  L <- L[,Convergence==1] #modif le 16/7/15
  C <- C[,Convergence==1] #modif le 16/7/15
  Erreur <- as.vector(Erreur)		#modif le 16/7/15
  L <- as.vector(L)		#modif le 16/7/15
  C <- as.vector(C)		#modif le 16/7/15
  if(length(Erreur)==0) {
    print("pb de convergence dans rirls")
    lambdaopt <- lambda.grid[length(lambda.grid)]
    ncompopt <- 1
  }
  if(length(Erreur)>0) {
    ii<- which.min(Erreur)
    lambdaopt <- L[ii]
    ncompopt <- C[ii]			
  }
  return(list(lambdaopt=lambdaopt,ncompopt=ncompopt,p.cvg=p.cvg/folds))
}

cv.irlspls<-function(Y,D,X,ncompmax,nbrIterMax=15,threshold=10^(-12),folds,proportion){
  ncompopt<-1
  if (ncompmax!=1){
    Erreur<-rep(0,ncompmax)
    p.cvg<-rep(0,ncompmax)
    for (t in 1:folds){
      set.seed(200)
      ind<-preselected.sample(Y,trunc(proportion*length(Y))) 
      cvXL<-X[ind,]
      cvDL<-D[ind,]
      cvYL<-Y[ind]
      cvXT<-X[-ind,]
      cvDT<-D[-ind,]
      cvYT<-Y[-ind]	
      for (k in (1:ncompmax))
      {
        fit<-fit.irlspls(Y=cvYL,D=cvDL,X=cvXL,ncomp=k,nbrIterMax=nbrIterMax,threshold=threshold)  
        #res <- LSPLSirls(cvYL,D=cvXL,X=cvXLc,ncomp=k,NbrIterMax=15,Threshold=10^(-12))
        pred<-predict.irlspls(fit,newX=cvXT,newD=cvDT)
        Erreur[k]<-Erreur[k]+sum(abs(pred$newY-cvYT))
        p.cvg[k]<-p.cvg[k]+fit$cvg
      }
    }	
    ncompopt<-which.min(Erreur)
  }
  return(list(ncompopt=ncompopt,p.cvg=p.cvg))
}
cv.lspls.glm<-function(Y,X,D,ncompmax,folds=5,proportion=0.9,method=c("LS-PLS-IRLS", "R-LS-PLS", "IR-LS-PLS"),lambda.grid=NULL,penalized=NULL,nbrIterMax=NULL,threshold=NULL){
  
  if (method=="LS-PLS-IRLS"){
    cv<-cv.lspls(Y=Y,X=X,D=D,ncompmax=ncompmax,folds=folds,proportion=proportion)
  }
  if (method=="R-LS-PLS"){
    if(is.null(lambda.grid)==TRUE){
      lambda.grid<-exp(seq(-3*log(10),2*log(10),length.out = 6))
    }
    if(is.null(nbrIterMax)==TRUE){
      nbrIterMax<-15
    }
    if(is.null(threshold)==TRUE){
      threshold<-10^(-12)
    }
    cv<-cv.rlspls(Y=Y,X=X,D=D,ncompmax=ncompmax,folds=folds,proportion=proportion,lambda.grid=lambda.grid,penalized=penalized,nbrIterMax=nbrIterMax,threshold=threshold)
  }
  if (method=="IR-LS-PLS"){
    if(is.null(nbrIterMax)==TRUE){
      nbrIterMax<-15
    }
    if(is.null(threshold)==TRUE){
      threshold<-10^(-12)
    }
    cv<-cv.irlspls(Y=Y,X=X,D=D,ncompmax=ncompmax,folds=folds,proportion=proportion,nbrIterMax=nbrIterMax,threshold=threshold)
  }
  return(cv)
}
