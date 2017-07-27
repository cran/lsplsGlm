pls <- function (Y,X,W=diag(rep(1,length(Y))),ncomp){
  
  ##    INPUT VARIABLES
  #########################
  ##  X   : matrix n x p
  ##      train data matrix
  ##  y : vector 
  ##      response variable real vector of length n
  ##  W : matrix n x n
  ##    weigth matrix
  ## ncomp : maximal pls components
  
  
  ##    OUTPUT VARIABLES
  ##########################
  ## Coefficients : matrix p x ncomp
  ##                regression coefficients w.r.t. the columns of X
  ##                for number of composants = 1:ncomp
  ## CoefPLS : matrice qui permet de calculer les nouvelles composantes PLS
  
  ## RUN PLS
  ###########
  #Some initializations
  n <- dim(X)[1]
  r <- dim(X)[2]
  
  ## Run PLS
  #Initialize some variables
  PsiAux <- diag(c(rep(1,r)))
  ##scale y and X
  E  <- scale(X)
  f <- scale(Y)
  Omega <- matrix(0,r,ncomp)
  Scores <- matrix(0,n,ncomp)
  TildePsi <- matrix(0,r,ncomp)
  Loadings <- matrix(0,r,ncomp)
  qcoeff <- vector(ncomp,mode="numeric")
  GAMMA <- matrix(0,nrow=r,ncol=ncomp)
  
  #WPLS loop
  for (count in 1:ncomp)
  {Omega[,count]<-t(E)%*%W%*%f/sqrt(sum((t(E)%*%W%*%f)^2))
  #au lieu de Omega[,count]<-t(E)%*%W%*%f c'est pour avoir les mêmes scores que pour NIPALS
  #Score vector
  t<-E%*%Omega[,count]
  c<-t(t)%*%W%*%t
  Scores[,count]<-t #modif : faut-il normer? Scores[,count]<-t/sqrt(c[1,1])
  TildePsi[,count] <- PsiAux%*%Omega[,count]
  #Deflation of X
  Loadings[,count]<-t(t(t)%*%W%*%E)/c[1,1]
  E<-E-t%*%t(Loadings[,count])
  #Deflation of f1
  qcoeff[count]<-t(W%*%f)%*%t/c[1,1]
  f <- f-qcoeff[count]*t
  #Recursve definition of RMatrix
  PsiAux<-PsiAux%*%(diag(c(rep(1,r)))-Omega[,count]%*%t(Loadings[,count]))
  #Express regression coefficients w.r.t. the columns of [X] for ncomp=count
  if (count==1)
    GAMMA[,count]<-TildePsi[,1:count]%*%t(c(qcoeff[1:count]))
  if (count!=1)
    GAMMA[,count]<-TildePsi[,1:count]%*%c(qcoeff[1:count])
  }
  
  siginv <- diag(1/attr(scale(X),"scaled:scale"))
  stdY <- attr(scale(Y),"scaled:scale")
  GAMMAX <- stdY*siginv%*%GAMMA
  ALPHA <- attr(scale(Y),"scaled:center")-t(attr(scale(X),"scaled:center"))%*%GAMMAX
  TildePsiX <- stdY*siginv%*%TildePsi
  ScoresX <- X%*%TildePsiX	
  List <- list(coefficients=GAMMAX,projection=TildePsiX,scores=ScoresX,intercept=ALPHA)
  #alpha=constante
  #GAMMAX => coefficient de T
  #scores => T
  return(List)
  
}