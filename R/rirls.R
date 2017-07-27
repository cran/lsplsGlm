rirls <- function(Y,D,X,lambda=0,penalized=TRUE,nbrIterMax=15,threshold=10^(-12)) {
  ##    IN     
  ##########
  ##  D : Design matrix
  ##          matrix n x q 
  ##  X : Spectral design matrix
  ##          matrix n x p  
  ##  Y : vector n 
  ##      response variable {0,1}-valued vector
  ##  lambda :  coefficient of the ridge penalty
  ##          real
  ##  nbrIterMax : Maximal number of iterations
  ##          integer
  ##  threshold : Used for the stopping rule.
  ##          real
  ##  penalized : if TRUE the parameter associated with D are ridge penalized
  
  ##    OUT     
  ############
  ##  out :  structure that contains the fields
  ##      Gamma : vector of the regression coefficients w.r.t the design
  ##      matrix Z=[1 X].
  ##      Cvg : Cvg=1 if the algorithm has converged otherwise 0.
  ##		W : matrix nxn weight at convergence
  ##		pseudoresponse : matrix nx1 pseudoresponse at convergence
  
  p <- dim(X)[2]
  q <- dim(D)[2]
  n <- dim(X)[1]
  Z <- cbind(rep(1,n),D,X)
  dZ2 <- dim(Z)[2]
  
  R<-matrix(0,dZ2,dZ2) 	
  if(is.null(D)==TRUE)
  {q<-0} 
  if (penalized==TRUE)
  {R[2:dZ2,2:dZ2]<-diag(1,nrow=(p+q))} 
  
  R[(q+2):dZ2,(q+2):dZ2]<-diag(1,nrow=p) 
  
  ##  THE NR ALGORITHM
  ##################################
  ##  1. Initialize the parameter
  
  c <- log(3) 
  Eta <- c*Y-c*(1-Y)                
  mu <- 1/(1+exp(-Eta))             
  DiagWNR <- mu*(1-mu)              
  WT <- diag(c(DiagWNR))         
  H <- t(Z)%*%WT%*%Z+lambda*R                
  trysolve<-try(solve(H,t(Z)%*%(diag(c(DiagWNR))%*%Eta+(Y-mu))), silent = TRUE)
  if (sum(is.nan(trysolve))>0) {
    trysolve <- NULL}
  
  if (is.matrix(trysolve)==FALSE) {
    Gamma <- rep(1,dZ2) 
    StopRule <- 1   
    NbrIter <- nbrIterMax+1
    Illcond <- 1
    Separation <- 0
    Cvg <- 0}
  if (is.matrix(trysolve)==TRUE)  {
    StopRule <- 0 
    NbrIter <- 0
    Illcond <- 0
    Separation <- 0}              
  
  ##  2. Newton-Raphson loop
  
  while (StopRule==0)
  {#Increment the iterations number 
    NbrIter <- NbrIter + 1
    #Udapte Gamma
    if (NbrIter==1)
    {Gamma <- trysolve}
    if (NbrIter>1)
    {Gamma <- Gamma+trysolve}           
    #Udapte Eta            
    Eta <- Z%*%Gamma   
    Eta[Eta>10]=10
    Eta[Eta< -10]=-10
    #Udapte mu
    mu<-1/(1+exp(-Eta))   
    #Udapte Weight   
    DiagWNR <- mu*(1-mu)  
    WT <- diag(c(DiagWNR))
    #Udapte Gradient
    Gradient <- t(Z)%*%(Y-mu)-lambda*R%*%Gamma
    #Udapte H           
    H <- t(Z)%*%WT%*%Z+lambda*R   
    trysolve <- try(solve(H,Gradient), silent = TRUE)  
    if (sum(is.nan(trysolve))>0) {
      trysolve <- NULL}       
    
    #Compute the StopRule
    #on the (quasi)-separation detection (for lambda=0)
    if (lambda==0)
    {Separation <- as.numeric(sum((Eta>0)+0==Y)==n)}
    #on the conditioning of matrix
    Illcond <- as.numeric(is.matrix(trysolve)==FALSE) 
    #on the convergence
    Cvg <- as.numeric(sqrt(sum(abs(Gradient)^2))<=threshold)
    StopRule <- max(Cvg,as.numeric(NbrIter>=nbrIterMax),Separation,Illcond) 
  }
  
  ##  CONCLUSION
  ###############        
  
  list <- list(coefficients=Gamma,cvg=max(Cvg,Separation),W=WT,pseudoresponse=(Eta+(Y-mu)/mu/(1-mu)))
  return(list)
}

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
}