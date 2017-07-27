\name{fit.lspcr.glm}
\alias{fit.lspcr.glm}

\title{
  Fitting a LS-PCR model for logistic regression
}
\description{
Fits a model the combination of two methods: Ordinary Least Square (OLS) and Principal Component Regression (PCR) to fit both clinical and gene expression data.
}
\usage{
  fit.lspcr.glm(Y,X,D,ncomp)
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
  \item{ncomp}{
    a positive integer. \code{ncomp} is the number of selected components.
  }
}
\details{
 This function combines two methods, the first one is the Principal Components Regression on genetic data to reduce the dimension using \code{prcomp} from \{\code{stats}\} package. The second one is the logistic regression on the concatenation of the \code{ncomp} first selected axes and clinical data (\code{D}) to explain \code{Y}. To do that we use \code{glm} from \{\code{stats}\} package.
}
\value{
  \item{coefficients }{coefficients of logistic regression.}
  \item{cvg }{the 0-1 value indicating convergence of the algorithm (1 for convergence, 0 otherwise).}
  \item{projection }{projection matrix used to convert \code{X} to scores.}
}

\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}

\seealso{\code{\link{cv.lspcr.glm}},\code{\link{predict.lspcr.glm}}.}

\examples{ 
#Data
data(BreastCancer)
#Vector of response
Y<-BreastCancer$Y
#Genetic data
X<-BreastCancer$X
#Clinical data
D<-BreastCancer$D
#Apply fit.lspcr.glm with ncomp=5 using the 76 first patients
fit<-fit.lspcr.glm(Y=Y[1:76],X=X[1:76,],D=D[1:76,],ncomp=5)

#using projection to predict class of 2 last patients
newX<-X[77:78,]
newD<-D[77:78,]

#New Score matrix
newScores<-newX\%*\%fit$projection

#prediction
newEta=cbind(rep(1,dim(newD)[1]),newD,newScores)\%*\%fit$coefficients
newPi=1/(1+exp(-newEta))
newY=as.numeric(newEta>0)
}

