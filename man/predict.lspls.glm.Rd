\name{predict.lspls.glm}
\alias{predict.lspls.glm}

\title{
Predict method for LS-PLS model fits.
}
\description{
Obtains predictions and prediction probabilities from a fitted LS-PLS object.
}
\usage{
\method{predict}{lspls.glm}(object,newX,newD,...)
}

\arguments{
  \item{object}{
results from \code{\link{fit.lspls.glm}} function.
}

  \item{newX}{
new matrix of clinical data. 
}

  \item{newD}{
new matrix of genetic data.
}
\item{...}{further arguments.  Currently not used.}
}
\details{
This function is used to obtained predicted values using a model fitting with \code{\link{fit.lspls.glm}}. It returns predictions and prediction probabilites in case that the user wants to choose his own decision rule. By default (in the prediction vector) if probabilities are under 0.5 observations are in the 0 class and if probabilities are upper than 0.5 observations are in the 1 class.
}
\value{
\item{newY }{the \code{newY} containing the \code{n} predicted values of the response variables for the observations from \code{newX} and \code{newD}.}
\item{newPi }{the \code{newPi} containing the \code{n} probabilities of the response variables for the observations from \code{newX} and \code{newD}.}
\item{cvg }{the 0-1 value indicating convergence of the algorithm (1 for convergence, 0 otherwise).}
}

\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}

\seealso{\code{\link{fit.lspls.glm}}.}

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

#split data, 2/3 of our datasets are used to compose learning samples
index<-preselected.sample(Y,trunc(2*length(Y)/3)) 
XL<-X[index,]
DL<-D[index,]
YL<-Y[index]

#1/3 is four our test samples
XT<-X[-index,]
DT<-D[-index,]
YT<-Y[-index]

#cross validation to find the optimal number of component
#method LS-PLS-IRLS
ncompopt.lsplsirls<-cv.lspls.glm(Y=YL,X=XL,D=DL,folds=5,ncompmax=5,
                      proportion=0.7,method="LS-PLS-IRLS")$ncompopt
#method R-LS-PLS
cv<-cv.lspls.glm(Y=YL,X=XL,D=DL,ncompmax=5,proportion=0.7,method="R-LS-PLS",
                      lambda.grid=exp(log(10^seq(-3,2,0.7))),penalized=TRUE,
		      nbrIterMax=15,threshold=10^(-12))
ncompopt.rlspls<-cv$ncompopt
lambdaopt.rlspls<-cv$lambdaopt
#method IR-LS-PLS
ncompopt.irlspls<-cv.lspls.glm(Y=YL,X=XL,D=DL,ncompmax=5,proportion=0.7,method="IR-LS-PLS",
                                                 nbrIterMax=15,threshold=10^(-12))$ncompopt  

#fits
fit.lsplsirls<-fit.lspls.glm(Y=YL,X=XL,D=DL,ncomp=ncompopt.lsplsirls,
                                                method="LS-PLS-IRLS")
fit.rlspls<-fit.lspls.glm(Y=YL,X=XL,D=DL,ncomp=ncompopt.rlspls,method="R-LS-PLS",
                                          lambda=lambdaopt.rlspls,penalized=TRUE,
                                          nbrIterMax=15,threshold=10^(-12))
fit.irlspls<-fit.lspls.glm(Y=YL,X=XL,D=DL,ncomp=ncompopt.irlspls,method="IR-LS-PLS",
                                                   nbrIterMax=15,threshold=10^(-12))

#predictions
pred.lsplsrils<-predict.lspls.glm(fit.lsplsirls,newX = XT,newD = DT)
pred.rlspls<-predict.lspls.glm(fit.rlspls,newX = XT,newD = DT)
pred.irlspls<-predict.lspls.glm(fit.irlspls,newX = XT,newD = DT)

#Confusion mmatrix
table(YT,pred.lsplsrils$newY)
table(YT,pred.rlspls$newY) 
table(YT,pred.irlspls$newY) 
}
}

