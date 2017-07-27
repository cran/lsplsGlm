\name{CentralCNS}
\alias{CentralCNS}

\title{Gene expression and clinical data used to predict tumors of Central Nervous System from children}
\usage{
data(CentralCNS)
}

\description{
Gene expression (\code{60x7129}) and clinical data (\code{60x4}) are used to predict the response of childhood malignant embryonal tumors of Central Nervous System (CNS) to therapy.
}

\details{
The dataset is composed of 60 patients samples, 21 patients died and 39 survived within 24 months. There are two matrices: genetic data X of size \code{60x7129} and clinical data D of size \code{60x4}. There is also  a response variable \code{Y} of size \code{60x1}.
}

\value{
  \item{X}{a data gene matrix (\code{60x7129}) giving the expression levels of 7129 
  genes for the 60 patients. Each row corresponds to a patient, each column to a
  gene.}
  \item{D}{a data matrix (\code{60x4}) of clinical data. Each row
    corresponds to a patient and each column to a clinical variable.}
  \item{Y}{a numeric vector of length 60 giving the condition of the patient (1 dead, 0 alive).} 
  

}

\source{S. L. Pomeroy, P. Tamayo, and M. Gaasenbeek. Prediction of central nervous system
embryonal tumour outcome based on gene expression. Nature, 415:436-442, 2002.}


\examples{
# load dataset
data(CentralCNS)

# how many patients and how many genes ?
dim(CentralCNS$X)

# how many patients of class 0 and 1 respectively ?
sum(CentralCNS$Y==0)
sum(CentralCNS$Y==1)
}
