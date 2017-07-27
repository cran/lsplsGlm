\name{cv.lspcr.glm}
\alias{cv.lspcr.glm}

\title{
Cross-validation for LS-PCR model for logistic regression  
}
\description{
Finds the optimal number of component for LS-PCR model for logistic regression.
}
\usage{
cv.lspcr.glm(Y, X, D, ncompmax, folds = 5, proportion = 0.9)
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
}
\details{
This function finds the optimal number of component for a LS-PCR model. At each cross validation run, \code{X}, \code{D} and \code{Y} are split into one training set 
and one test set (of proportion \code{proportion} and \code{1-proportion}). Then the classification error rate is computed for each value of \code{ncomp} between 1 and \code{ncompmax}. At the end we choose the number of component for which the classification error rate is minimal. This function returns also \code{p.cvg}. It's a vector of size \code{ncompmax} which contains convergence proportion of the logistic regression for each number of component between 1 and \code{ncompmax}.
}
\value{
 \item{ncompopt }{the optimal number of component.}		
 \item{p.cvg }{convergence proportion of the logistic regression.}
}

\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}
\seealso{\code{\link{fit.lspcr.glm}}.}

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
#cross validation to find the optimal number of component
cv<-cv.lspcr.glm(Y=Y,X=X,D=D,folds=5,ncompmax=5,proportion=0.9)
ncompopt<-cv$ncompopt
}
}

