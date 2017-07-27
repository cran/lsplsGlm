preselected.sample <- function (label,ns){
  #sample of size ns with the same proportion of 1
  #label must be 0 or 1
  nl <- length(label)
  p2 <- sum(label==1)/nl
  n2 <- trunc(ns*p2)
  n1 <- ns-n2
  ind <-sort(c(sample((1:nl)[label==0],n1,replace=FALSE),sample((1:nl)[label==1],n2,replace=FALSE)))
  return(ind)	
}