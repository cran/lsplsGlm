lspls<-function(Y, D, X, W=diag(rep(1,nrow(D))), ncomp){
  #Function LS-PLS
  
  ##    INPUT VARIABLES
  #########################
  ##  D  : design matrix n x p (donn?es cliniques)
  ##
  ##  y : vector
  ##      response variable real vector of length n
  ##
  ##  X : spectral matrix  (donnees genetiques) n x c
  ##
  ##  W : weigth matrix n x n
  ##  	
  ## ncomp : maximal pls components
  ##
  
  ########################################
  ##    OUTPUT VARIABLES
  #########################
  ##
  ##
  ## coefficients.PLS : matrix with final coefficients predictions
  ##
  ## predictors : matrix with variables (design matrix) and scores from PLS regression
  ##                on the residual (OLS Step between Y on the design matrix D) and the centred X matrix (spectral matrix)
  ##
  ## scores :  matrix of pls scores used in the final regression
  ##
  ## loadings : matrix of pls loadings
  ##
  ## prediction : matrix DT 
  ##
  ## orthCoef : coefficients matrix c x p returned also in the function lspls (Bjorn-Helge Mevik) to be used to predict new predictors (solve(((t(D)%*%D))%*%D%*%X_cent)
  ## For new matrix test 
  ## Ttest=(Xtest-cbind(rep(1,dim(Dtest)[1]),Dtest)%*%modpls$orthCoef)%*%modpls$loadin
  ## Ypred=cbind(Dtest,Ttest)%*%modpls$coefficients+modpls$alpha in the gaussian case
  ##
  ############################################
  #estimation of the regression coefficients
  Dc <- cbind(rep(1,dim(D)[1]),D)
  reg_coef<- solve(t(Dc)%*%W%*%Dc)%*%t(Dc)%*%W%*%Y
  
  #computation of the residuals
  residual<-Y-D%*%reg_coef[-1]-reg_coef[1]
  
  #orthogonalisation de X_cent par rapport a D
  orthCoef<-solve(t(Dc)%*%W%*%Dc)%*%t(Dc)%*%W%*%X  
  X_cent_orth<-X-Dc%*%orthCoef
  
  #pls fit		
  fit<-pls(Y=residual,X=X_cent_orth,W=W,ncomp=ncomp) 
  
  #Merge the matrix design D ans the score matrices into a new matrices  K
  K<-cbind(D, fit$scores)
  
  #OLS regression of Y on the combined K matrices
  Kc <- cbind(rep(1,dim(K)[1]),K)
  reg_coef<- solve(t(Kc)%*%W%*%Kc)%*%t(Kc)%*%W%*%Y
  
  return(list(coefficients=reg_coef[-1],predictors=K, scores=fit$scores,projection=fit$projection, orthCoef=orthCoef,intercept=reg_coef[1]))
}
