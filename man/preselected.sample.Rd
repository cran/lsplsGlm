\name{preselected.sample}
\alias{preselected.sample}
\title{
 Selected randomized controlled random sample
}
\description{
  Creates a random sample keeping the proportions of 1 and 0 that are in \code{label}.
}
\usage{
  preselected.sample(label, ns)
}

\arguments{
  \item{label}{
    vector of 0 and 1.
    }
  \item{ns}{
    sample size, contain \code{ns} random index of \code{label}.\code{ns} has to be lower than  length of \code{label}.
    }
}
\details{
  The aim of this method is to select randomly \code{ns} index of \code{label}. This function returns a vector of size \code{ns} composed of random rank of \code{label}.
}
\value{
  \item{index.learn }{a vector of size \code{ns} containing random index of \code{label}.}
}

\author{
Caroline Bazzoli, Thomas Bouleau, Sophie Lambert-Lacroix
}

\examples{
#load data
data(BreastCancer)
index<-preselected.sample(BreastCancer$Y,10)
index
BreastCancer$Y[index]      
}
