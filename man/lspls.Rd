\name{lspls}
\alias{lspls}
\title{
Weighted LS-PLS gaussian regression
}
\description{
Performs a weighted Least Square-Partial Least Square gaussian regression for both clinical and genetic data.
}
\usage{
  lspls(Y, D, X, W=diag(rep(1,nrow(D))), ncomp)
}
\arguments{
  \item{Y}{
    a vector of length \code{n} giving the classes of the \code{n} observations. \code{Y} contains continuous values.
}
  \item{X}{
   a data matrix (\code{nxp}) of genes. NAs and Inf are not allowed. Each row
    corresponds to an observation and each column to a gene.
  }
  \item{D}{
   a data matrix (\code{nxq}) of clinical data. NAs and Inf are not allowed. Each row
    corresponds to an observation and each colone to a clinical variable.
}

  \item{W}{
    weight matrix, if \code{W} is the identity matrix then the function will perform a standard LS-PLS regression.
  }
  \item{ncomp}{
    a positive integer. \code{ncomp} is the number of selected components.
  }
}
\details{
This function is a combination of Least Squares (LS) and Partial Least Square (PLS)[1]. This is an iterative procedure: the first
step is to use OLS on \code{D} to predict \code{Y}. New estimates for the residuals of \code{Y} on \code{D} are calculated from this regression and the algorithm is repeated until convergence. 
Here we use the orthogonalised variant. To do that we create a new matrix which is the projection of the matrix \code{X} into a space orthogonal to the space spanned by the design variables of \code{D}.
The standard PLS regression is then used on this new matrix instead of \code{X} [2].
}
\value{
\item{predictors }{matrix which combines \code{D} and scores from PLS regression}
\item{projection }{the projection matrix used to convert \code{X} to scores.}
\item{orthCoef }{the coefficients matrix of size \code{pxq} to be used to compute new predictors.}
\item{coefficients }{an array of PLS regression coefficients (\code{(p+1)xncomp})}
\item{intercept }{the constant of the model.}
}
\references{
[1] Jørgensen, K., Segtnan, V., Thyholt, K., and Næs, T. (2004). A comparison of
methods for analysing regression models with both spectral and designed variables.
Journal of Chemometrics, 18(10), 451-464.

[2] Caroline Bazzoli, Sophie Lambert-Lacroix. Classification using LS-PLS with logistic regression based on both clinical 
and gene expression variables. 2017. <hal-01405101>
}
\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}

\examples{
#X simulation
meanX<-sample(1:300,50)
sdeX<-sample(50:150,50)
X<-matrix(nrow=60,ncol=50)
for (i in 1:50){
  X[,i]<-rnorm(60,meanX[i],sdeX[i])
}

#D simulation
meanD<-sample(1:30,5)
sdeD<-sample(1:15,5)
D<-matrix(nrow=60,ncol=5)
for (i in 1:5){
  D[,i]<-rnorm(60,meanD[i],sdeD[i])
}

#Y simulation
Y<-rnorm(60,30,10)

# Learning sample
index<-sample(1:length(Y),round(2*length(Y)/3))
XL<-X[index,]
DL<-D[index,]
YL<-Y[index]

#fit the model
fit<-lspls(YL,X=XL,D=DL,ncomp=3,W=diag(rep(1,length(YL))))

#Testing sample
newX=X[-index,]
newD<-D[-index,]

#predictions with the constant of the model
a.coefficients<-c(fit$intercept,fit$coefficients)

#predictions
newZ=(newX-cbind(rep(1,dim(newD)[1]),newD)\%*\%fit$orthCoef)\%*\%fit$projection
newY=cbind(rep(1,dim(newD)[1]),newD,newZ)\%*\%a.coefficients

}


