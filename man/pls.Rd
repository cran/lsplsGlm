\name{pls}
\alias{pls}
\title{
  Weighted PLS gaussian regression
}
\description{
  Performs a weighted Partial Least Square gaussian regression.
}
\usage{
pls(Y,X,W = diag(rep(1, length(Y))),ncomp)
}
\arguments{
  \item{Y}{
  a vector of length \code{n} giving the classes of the \code{n} observations. \code{Y} contains continuous values.
  }
  \item{X}{
  a data matrix (\code{nxp}) of genes. NAs and Inf are not allowed. Each row
  corresponds to an observation and each column to a gene.
  }
  \item{W}{
    weight matrix, if \code{W} is the identity matrix then the function will perform a standard PLS regression.
  }

 \item{ncomp}{
    a positive integer. \code{ncomp} is the number of PLS components.
  }
}

\details{
This function performs a weighted PLS gaussian regression. It takes as input a vector of response, a data matrix about genes, a number of component and a weight matrix. 
If weight matrix is the identity matrix then it performs a standard PLS regression.
}
\value{

\item{coefficients }{an array of regression coefficients (\code{(p+1)xncomp}).}
\item{projection }{the projection matrix, used to convert \code{X} to scores.}
\item{scores }{the scores matrix \code{(nxncomp)} of PLS regression.}
\item{intercept }{the constant of the model.}
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

#Y simulation
Y<-rnorm(60,30,10)

# Learning sample
index<-sample(1:length(Y),round(2*length(Y)/3))
XL<-X[index,]
YL<-Y[index]

#fit the model
fit<-pls(Y=YL,X=XL,ncomp=3)

#Testing sample
newX=X[-index,]

#predictions with the constant of the model
a.coefficients<-rbind(fit$intercept,fit$coefficients)

#predictions
newY=cbind(rep(1,dim(newX)[1]),newX)\%*\%a.coefficients

}

