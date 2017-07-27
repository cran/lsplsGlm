\name{fit.lspls.glm}
\alias{fit.lspls.glm}

\title{
Fitting LS-PLS for generalized model for logistic regression  
}
\description{
Fits a Least Square-Partial Least Square for logistic regression model. There are 3 extensions.
}
\usage{
fit.lspls.glm(Y, X, D, W=diag(rep(1,nrow(D))), ncomp, method=c("LS-PLS-IRLS", 
"R-LS-PLS", "IR-LS-PLS"),lambda=NULL, penalized = NULL, nbrIterMax = NULL, 
threshold= NULL)
}

\arguments{
  \item{Y}{
 a vector of length \code{n} giving the classes of the \code{n} observations. The classes
    must be coded as 1 or 0.
}
  \item{X}{
    a data matrix (\code{nxp}) of genes. NAs and Inf are not allowed. Each row
    corresponds to an observation and each column to a gene.
}
  \item{D}{
    a data matrix (\code{nxq}) of clinical data. NAs and Inf are not allowed. Each row
    corresponds to an observation and each column to a clinical variable.
}
  \item{W}{
    weight matrix. If \code{W} is the identity matrix then the function will fit a standard LS-PLS model.
}
  \item{ncomp}{
 a positive integer. \code{ncomp} is the number of selected components.
}
  \item{method}{
 one of the 3 extensions of LS-PLS for logistic regression models. 
}
  \item{penalized}{
 if TRUE the parameter associated with \code{D} is ridge penalized. To use only if \code{method} equals to "R-LS-PLS". 
}
  \item{lambda}{
 coefficient of ridge penalty. If \code{penalized = TRUE}, \code{lambda} is the penalty for \code{D}. To use only if \code{method} equals to "R-LS-PLS".
}

 \item{nbrIterMax}{
 maximal number of iterations. To use only if \code{method} equals to "R-LS-PLS" or "IR-LS-PLS".
}
  \item{threshold}{
 used for the stopping rule. To use only if \code{method} equals to "R-LS-PLS" or "IR-LS-PLS".
}

}
\details{
This function fits LS-PLS models. With the argument "method" the user can choose one of the three extensions of LS-PLS for logistic regression (LS-PLS-IRLS, R-LS-PLS, IR-LS-PLS). 
For more details see references.

}
\value{

  \item{coefficients }{vector of length \code{q+p} associate to \code{cbind(D,X)}.}
  \item{cvg }{the 0-1 value indicating convergence of the algorithm (1 for convergence, 0 otherwise).}
  \item{orthCoef }{coefficients matrix (\code{pxq}) returned also in the function \code{\link{lspls}} to be used to compute new predictors.}
  \item{projection }{the projection matrix used to convert \code{X} to scores.}
  \item{intercept }{the constant of the model.}
}
\references{
Caroline Bazzoli, Sophie Lambert-Lacroix. Classification using LS-PLS with logistic regression based on both clinical 
and gene expression variables. 2017. <hal-01405101>
}
\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}

\seealso{\code{\link{cv.lspls.glm}},\code{\link{predict.lspls.glm}}.}

\examples{
#Data
data(BreastCancer)
#vector of responses
Y<-BreastCancer$Y
#Genetic data
X<-BreastCancer$X
#Clinical data
D<-BreastCancer$D
#SIS selection
X<-scale(X)
\donttest{
X<-SIS.selection(X=X,Y=Y,pred=50)
#Cross validation, 90% of our datasets is used to compose learning samples

#method LS-PLS-IRLS
ncompopt.lsplsirls<-cv.lspls.glm(Y=Y,X=X,D=D,folds=5,ncompmax=5,proportion=0.9,
		                                 method="LS-PLS-IRLS")$ncompopt
#method R-LS-PLS
cv<-cv.lspls.glm(Y=Y,X=X,D=D,ncompmax=5,proportion=0.9,method="R-LS-PLS",
	           lambda.grid=exp(log(10^seq(-3,2,0.7))),penalized=TRUE,
	           nbrIterMax=15,threshold=10^(-12))
ncompopt.rlspls<-cv$ncompopt
lambdaopt.rlspls<-cv$lambdaopt
#method IR-LS-PLS
ncompopt.irlspls<-cv.lspls.glm(Y=Y,X=X,D=D,ncompmax=5,proportion=0.9,
	method="IR-LS-PLS",nbrIterMax=15,threshold=10^(-12))$ncompopt         

#fitting model
fit.lsplsirls<-fit.lspls.glm(Y=Y,X=X,D=D,ncomp=ncompopt.lsplsirls,method="LS-PLS-IRLS")
fit.rlspls<-fit.lspls.glm(Y=Y,X=X,D=D,ncomp=ncompopt.rlspls,method="R-LS-PLS",
		         lambda=lambdaopt.rlspls,penalized=TRUE,nbrIterMax=15,
                         threshold=10^(-12))
fit.irlspls<-fit.lspls.glm(Y=Y,X=X,D=D,ncomp=ncompopt.irlspls,method="IR-LS-PLS",
                                                nbrIterMax=15,threshold=10^(-12))

}
}
