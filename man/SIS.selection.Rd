\name{SIS.selection}
\alias{SIS.selection}
\title{
  Sure Independence Screening
}
\description{
  SIS has been performed to select relevant gene expression variables. SIS ranks the importance of features according to their magnitude of marginal regression coefficients.
}
\usage{
SIS.selection(X,Y, pred, scale = F)
}
\arguments{
  \item{X}{
    a data matrix (\code{nxp}) of genes. NAs and Inf are not allowed. Each row
    corresponds to an observation and each column to a gene.
  }
  \item{Y}{
    a vector of length n giving the classes of the n observations. The classes
    must be coded as 1 or 0.
  }
  \item{pred}{
    number of relevant variable to select, \code{pred} has to be lower than p.
  }
  \item{scale}{
    If scale=TRUE, \code{X} will be scaled.
  }
}
\details{
  Sure Independence Screening (SIS) has been performed to select relevant gene expression
  variables \code{pred} such as \code{pred} < \code{p}. SIS refers to ranking features according to marginal
  utility, namely, each feature is used independently as a predictor to decide its usefulness
  for predicting the response. Precisely SIS ranks the importance of features according to
  their magnitude of marginal regression coefficients.
}

\value{
  Return a matrix (\code{nxpred}) with only the \code{pred} most relevant gene and all the observations
}
\references{
 Fan, J. and Lv, J. (2008). Sure independence screening for ultrahigh dimensional
feature space. Journal of the Royal Statistical Society, 70, 849-911.
}
\author{
  Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}

\examples{
data("BreastCancer")
X<-scale(BreastCancer$X)
Y<-BreastCancer$Y
\donttest{
Xsis<-SIS.selection(X,Y,50)
}
}
