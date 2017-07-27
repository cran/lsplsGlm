\name{BreastCancer}
\alias{BreastCancer}

\title{Gene expression and clinical data used to predict the presence of subclinical metastases for breast cancer patients}
\usage{
data(BreastCancer)
}

\description{
Genetic and clinical data about 78 primary breast cancers. BreastCancer is a list of two matrice (\code{X} and \code{D}) and one vector (\code{Y}).

  
}

\details{
The object is composed of a list of two matrices: genetic data X of size \code{78x4348} and clinical data D of size \code{78x7}. There is also  a response variable Y of size \code{78x1}. The clinical  data contain information on 78 primary breast cancers (34 from patients who developed metastases within 5 years and 44 from patients who continue to be disease-free after a period of at least 5 years) which have
been selected from patients who were lymph node negative and under 55 years of age at diagnosis.
}

\value{
  \item{X}{a data gene matrix (\code{78x4348}) giving the expression levels of 4348 
  genes for the 78 patients. Each row corresponds to a patient, each column to a
  gene.}
  \item{D}{a data matrix (\code{78x7}) of clinical data. Each row
    corresponds to a patient and each column to a clinical variable.}
  \item{Y}{a numeric vector of length 78 giving the presence of subclinical metastases (1 for presence, 0 otherwise).} 
  

}

\examples{
# load dataset
data(BreastCancer)

# how many patients and how many genes ?
dim(BreastCancer$X)

# how many patients of class 0 and 1 respectively ?
sum(BreastCancer$Y==0)
sum(BreastCancer$Y==1)
}
