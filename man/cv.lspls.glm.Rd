\name{cv.lspls.glm}
\alias{cv.lspls.glm}

\title{
Cross-validation for LS-PLS model for logistic regression  
}
\description{
Finds the optimal number of component for one of the three extesions of LS-PLS. Moreover it finds the \code{lambda} optimal for the R-LS-PLS method.
}
\usage{
cv.lspls.glm(Y, X, D, ncompmax, folds = 5, proportion = 0.9, method=c("LS-PLS-IRLS", 
"R-LS-PLS", "IR-LS-PLS"),lambda.grid=NULL,penalized=NULL, 
nbrIterMax=NULL,threshold=NULL)
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
  \item{ncompmax}{
 a positive integer. \code{ncompmax} is the maximal number of selected components.
}
  \item{folds}{
a positive integer indicating the number of folds in K-folds cross-validation procedure.
}
  \item{proportion}{
proportion of the dataset in the learning sample. \code{proportion} has to be between 0 and 1.
}
  \item{method}{
one of the three extensions of LS-PLS for logistic regression models (LS-PLS-IRLS, R-LS-PLS, IR-LS-PLS).
}
  \item{lambda.grid}{vector of positif real (grid for ridge parameter). To use only if \code{method} equals to "R-LS-PLS". By default \code{lambda.grid = exp(log(10^seq(-3,2,0.7)))} }

  \item{penalized}{if TRUE the parameter associated with D is ridge penalized. To use only if \code{method} equals to "R-LS-PLS".}

   \item{nbrIterMax}{
 maximal number of iterations. To use only if \code{method} equals to "R-LS-PLS" or "IR-LS-PLS".
}
  \item{threshold}{
 used for the stopping rule. To use only if \code{method} equals to "R-LS-PLS" or "IR-LS-PLS".
}

}

\details{
This function finds the optimal number of component and the optimal lambda for a LS-PLS regression. At each cross validation run, \code{X}, \code{D} and \code{Y} are split into one training set and one test set (of proportion \code{proportion} and \code{1-proportion}). Then for each component between 1 and \code{ncompmax} (and for each value of \code{lambda.grid} if \code{method} equals to R-LS-PLS) classification error rate is determined. At the end we choose the \code{lambda} and the \code{ncomp} for which the classification error rate is minimal. This function returns also \code{p.cvg}. It's a vector of size \code{ncompmax} which contains convergence proportion for each number of component between 1 and \code{ncompmax}. For the method R-LS-PLS, p.cvg is a matrix of size \code{ncompmax x length(lambda.grid}).
}
\value{
 \item{ncompopt }{the optimal number of component.}
 \item{lambdaopt }{lambda optimal.}		
 \item{p.cvg }{convergence proportion.}

}

\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}

\seealso{\code{\link{fit.lspls.glm}}.}

\examples{

#data
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
	                          lambda.grid=exp(log(10^seq(-3,2,0.7))), 
                                  penalized=TRUE,nbrIterMax=15,
	                          threshold=10^(-12))
ncompopt.rlspls<-cv$ncompopt
lambdaopt.rlspls<-cv$lambdaopt
#method IR-LS-PLS
ncompopt.irlspls<-cv.lspls.glm(Y=Y,X=X,D=D,ncompmax=5,proportion=0.9,method="IR-LS-PLS",
                               nbrIterMax=15,threshold=10^(-12))$ncompopt    

}
}

