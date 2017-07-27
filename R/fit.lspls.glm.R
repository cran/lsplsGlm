fit.lspls<-function(Y,X,D,W=diag(rep(1,nrow(D))),ncomp){
  k=ncomp
  res <- lspls(Y=Y, D=D, X=X, W=W , ncomp=k)
  XM2 <- res$predictors
  projection <- res$projection
  orthCoef <- res$orthCoef
  res=glm(Y~XM2,family = binomial)
  z<-list(orthCoef=orthCoef,cvg=res$converged,projection=projection,coefficients=res$coefficients)
  class(z)<-"LS-PLS-IRLS"
  z
}

fit.rlspls<-function(Y,X,D,ncomp,lambda,penalized=TRUE,nbrIterMax=15,threshold=10^(-12)){
  k <- ncomp
  res <- rirls(Y=Y,D=D,X=X,penalized=penalized,lambda=lambda,nbrIterMax=nbrIterMax,threshold=threshold) 
  Ylspls <- res$pseudoresponse
  Dlspls <- D
  Xlspls <- X
  reslspls <- lspls(Y=Ylspls, D=Dlspls, X=Xlspls, W=res$W, ncomp=k)
  orthCoef <- reslspls$orthCoef
  projection <- reslspls$projection	
  z<-list(orthCoef=orthCoef, cvg=res$cvg, projection=projection,coefficients=reslspls$coefficients,intercept=reslspls$intercept)
  class(z)<-"R-LS-PLS"
  z
}

fit.irlspls <- function(Y,D,X,ncomp,nbrIterMax=15,threshold=10^(-12)) {
  ##    IN     
  ##########
  ##  D : Design matrix
  ##          matrix n x q 
  ##  X : Spectral design matrix
  ##          matrix n x p  
  ##  Y : vector n 
  ##      response variable {0,1}-valued vector
  ##  ncomp :  number of pls components
  ##          real
  ##  NbrIterMax : Maximal number of iterations
  ##          integer
  ##  Threshold : Used for the stopping rule.
  ##          real
  
  ##    OUT     
  ############
  ##  out :  structure that contains the fields
  ##      Coefficients : vector of the regression coefficients w.r.t the design
  ##      matrix Z=[1 D scoresPLS].
  ##      Cvg : Cvg=1 if the algorithm has converged otherwise 0.
  ##		Loadings matrix pxncompopt that allows to compute the PLS scores : TPLS=Xc%*%loadings
  
  ##  THE NR ALGORITHM
  ##################################
  ##  1. Initialize the parameter
  
  c <- log(3) 
  Eta <- c*Y-c*(1-Y)                
  mu <- 1/(1+exp(-Eta))             
  DiagWNR <- mu*(1-mu)              
  WT <- diag(c(DiagWNR))                  
  pseudoresponse <- Eta+(Y-mu)/mu/(1-mu)
  reslspls <- lspls(Y=pseudoresponse, D=D, X=X, W=WT, ncomp=ncomp)
  predictors <- cbind(rep(1,length(Y)),reslspls$predictors) ## 1 K	
  StopRule <- 0 
  NbrIter <- 0
  Separation <- 0           
  
  ##  2. Newton-Raphson loop
  
  while (StopRule==0)
  {#Increment the iterations number 
    NbrIter <- NbrIter + 1
    #Udapte Eta     
    Gammaold <- c(reslspls$intercept,reslspls$coefficients)
    Eta <- predictors%*%Gammaold  
    Eta[Eta>10]=10
    Eta[Eta< -10]=-10
    #Udapte mu
    mu<-1/(1+exp(-Eta))   
    
    #Udapte Weight   
    DiagWNR <- mu*(1-mu)  
    WT <- diag(c(DiagWNR))
    pseudoresponse <- Eta+(Y-mu)/mu/(1-mu)
    #Compute the new coefficients		
    
    reslspls <- lspls(Y=pseudoresponse, D=D, X=X, W=WT, ncomp=ncomp)
    predictors <- cbind(rep(1,length(Y)),reslspls$predictors) ## 1 K
    Gammanew <- c(reslspls$intercept,reslspls$coefficients)
    
    #Compute the StopRule
    #on the (quasi)-separation detection
    Separation <- as.numeric(sum((Eta>0)+0==Y)==length(Y))
    #on the convergence
    aux <- sum((Gammanew-Gammaold)^2)/sum((Gammaold)^2)
    Cvg <- as.numeric(sqrt(aux)<=threshold)
    StopRule <- max(Cvg,as.numeric(NbrIter>=nbrIterMax),Separation) 
  }
  
  ##  CONCLUSION
  ###############        
  
  z<-list(coefficients=Gammanew,cvg=max(Cvg,Separation),projection = reslspls$projection, orthCoef = reslspls$orthCoef)
  class(z)<-"IR-LS-PLS"
  z
}

fit.lspls.glm<-function(Y,X,D,W=diag(rep(1,nrow(D))),ncomp,method=c("LS-PLS-IRLS", "R-LS-PLS", "IR-LS-PLS"),lambda=NULL,penalized=NULL,nbrIterMax=NULL,threshold=NULL){
  
  if (method=="LS-PLS-IRLS"){
    fit<-fit.lspls(Y=Y,X=X,D=D,ncomp=ncomp,W=W)
  }
  if (method=="R-LS-PLS"){
    
    if(is.null(nbrIterMax)==TRUE){
      nbrIterMax<-15
    }
    if(is.null(threshold)==TRUE){
      threshold<-10^(-12)
    }
    fit<-fit.rlspls(Y=Y,X=X,D=D,ncomp=ncomp,lambda=lambda,penalized = penalized,nbrIterMax=nbrIterMax,threshold = threshold)
  }
  if (method=="IR-LS-PLS"){
    if(is.null(nbrIterMax)==TRUE){
      nbrIterMax<-15
    }
    if(is.null(threshold)==TRUE){
      threshold<-10^(-12)
    }
    fit<-fit.irlspls(Y=Y,X=X,D=D,ncomp=ncomp,nbrIterMax=nbrIterMax,threshold = threshold)
  }
  return(fit)
}
