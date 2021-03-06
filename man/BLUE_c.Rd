% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLUE.R
\name{BLUE_c}
\alias{BLUE_c}
\title{BLUEs of global location and scale parameters}
\usage{
BLUE_c(alpha_c, B_c, X_c)
}
\arguments{
\item{alpha_c}{the expectation of a combined standardized vector of ordered summary statistics, i.e. equation (3.21) in Yang et al. (2018).}

\item{B_c}{the variance-covariance matrix of a combined standardized vector of ordered summary statistics, i.e. equation (3.22) in Yang et al. (2018).}

\item{X_c}{a combined vector of ordered summary statistics.}
}
\description{
To obtain the global or overall best linear unbiased estimator (BLUE) of location and scale parameters (Yang et al., 2018).
}
\examples{
n1<-30 # sample sizes of three included studies
n2<-45
n3<-67
X1<-c(3,1.2) # the mean and standard deviation
X2<-c(1,4,10) # the sample mean, minimum and maximum values
X3<-c(1.5,3,5.5,8,12) # the sample mean, first and third quartiles, and minimum and maximum values
X_c<-c(X1[1],X2,X3)

alpha1<-0  #Approximate by the CLT.
B1<-1/sqrt(n1)
alpha2<-BLUE_s(X2,n2,"S1")$alpha
B2<-BLUE_s(X2,n2,"S1")$B
alpha3<-BLUE_s(X3,n3,"S3")$alpha
B3<-BLUE_s(X3,n3,"S3")$B

alpha_c<-c(alpha1,alpha2,alpha3)
B_c<-Matrix::bdiag(B1,B2,B3)

BLUE_c(alpha_c,B_c,X_c)

}
\references{
Yang X, Hutson AD, and Wang D. (2018). A generalized BLUE approach for combining location and scale information in a meta-analysis (Submitted).
}
