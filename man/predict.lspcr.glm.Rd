\name{predict.lspcr.glm}
\alias{predict.lspcr.glm}

\title{
Predict method for LS-PCR fits.
}
\description{
Obtains predictions and prediction probabilities from a fitted LS-PCR object.
}
\usage{
\method{predict}{lspcr.glm}(object,newX,newD,...)
}

\arguments{
  \item{object}{
results from \code{\link{fit.lspcr.glm}} function.
}
  \item{newX}{
new matrix of genetic data.
}

  \item{newD}{
new matrix of clinical data. 
}
\item{...}{further arguments.  Currently not used.}

}
\details{
This function is used to obtained predicted values using a model fitting with \code{\link{fit.lspcr.glm}}. It returns predictions and prediction probabilites in case that the user wants to choose his own decision rule. By default (in the prediction vector) if probabilities are under 0.5 observations are in the 0 class and if probabilities are upper than 0.5 observations are in the 1 class.
}
\value{

\item{newY }{the \code{newY} containing the \code{n} predicted values of the response variables for the observations from \code{newX} and \code{newD}.}
\item{newPi }{the \code{newPi} containing the \code{n} probabilities of the response variables for the observations from \code{newX} and \code{newD}.}
\item{cvg }{the 0-1 value indicating convergence of the algorithm (1 for convergence, 0 otherwise).}
}

\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix 
}

\seealso{\code{\link{fit.lspcr.glm}}.}

\examples{
#Data
data(BreastCancer)
#Vector of response
Y<-BreastCancer$Y
#Genetic data
X<-BreastCancer$X
#Clinical data
D<-BreastCancer$D

#Learn dataset and test dataset (80/20)
index<-preselected.sample(Y,round(length(Y)*0.8))
XL<-X[index,]
XT<-X[-index,]
DL<-D[index,]
DT<-D[-index,]
YL<-Y[index]
#Apply fit.lspcr.glm with 5 components
fit<-fit.lspcr.glm(Y=YL,X=XL,D=DL,ncomp=5)
#predictions
pred<-predict.lspcr.glm(fit,newD=DT,newX=XT)
}

