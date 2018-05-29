#' Luo et al.'s mean estimators
#'
#' Use Luo et al.'s methods (Luo et al., 2016) to estimate the sample mean based on sample quantiles derived summaries for a single study.
#' @param X a vector of ordered summary statistics
#' @param n the sample size
#' @param type a character string indicating which type of summary statistics is reported. The options for the \strong{type} argument are: 
#' \itemize{
#'    \item "S1" for the sample mean, minimum and maximum values
#'    \item "S2" for the sample mean, first and third quartiles
#'    \item "S3" for the sample mean, first and third quartiles, and minimum and maximum values
#' }
#
#' @examples
#' X<-c(1,4,10)
#' n<-30
#' type<-"S1"
#' Luo.mean(X,n,type)
#' @references
#' Luo D, Wan X, Liu J, and Tong T. (2016). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical Methods in Medical Research}, arXiv:1505.05687.
#' @export
Luo.mean<-function(X,n,type){
  
  if(type=="S1"){ # X is {a,m,b}
    sum_stat_index=c(1,floor(n/2)+1,n)
    nsum_stat=length(sum_stat_index)
    
    muhat=(4/(4+n^0.75))*(X[1]+X[3])/2+(n^0.75/(4+n^0.75))*X[2]
  }
  
  if(type=="S2") { # X is {q1,m,q3}
    sum_stat_index=c(floor(n/4)+1,floor(n/2)+1,floor(3*n/4)+1)
    nsum_stat=length(sum_stat_index)
    
    muhat=(0.7+(0.39/n))*(X[1]+X[3])/2+(0.3-(0.39/n))*X[2]
  }
  
  if(type=="S3") {  # X is 5 # sum stat
    sum_stat_index=c(1,floor(n/4)+1,floor(n/2)+1,floor(3*n/4)+1,n)
    nsum_stat=length(sum_stat_index)
    
    muhat=(2.2/(2.2+n^0.75))*(X[1]+X[5])/2+(0.7-(0.72/(n^0.55)))*(X[2]+X[4])/2+(0.3+(0.72/n^0.55)-(2.2/(2.2+n^0.75)))*X[3]
  }
  
  out=list(muhat=muhat,
           sum_stat_index=sum_stat_index,sum_stat=X)
  return(out)
}


#' Wan et al.'s standard deviation estimators
#'
#' Use Wan et al.'s methods (Wan et al., 2014) to estimate the sample standard deviation based on sample quantiles derived summaries for a single study.
#' @param X a vector of ordered summary statistics
#' @param n the sample size
#' @param type a character string indicating which type of summary statistics is reported. The options for the \strong{type} argument are: 
#' \itemize{
#'    \item "S1" for the sample mean, minimum and maximum values
#'    \item "S2" for the sample mean, first and third quartiles
#'    \item "S3" for the sample mean, first and third quartiles, and minimum and maximum values
#' }
#' 
#' @examples
#' X<-c(1,4,10)
#' n<-30
#' type<-"S1"
#' Wan.std(X,n,type)
#' @references
#' Wan X,Wang W, Liu J, and Tong T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC Medical Research Methodology}, \strong{14}:135.
#' @export
Wan.std<-function(X,n,type){
  
  if(type=="S1"){ # X is {a,m,b}
    sum_stat_index=c(1,floor(n/2)+1,n)
    nsum_stat=length(sum_stat_index)
    
    tao=2*stats::qnorm((n-0.375)/(n+0.25))
    sigmahat=(X[3]-X[1])/tao
  }
  
  if(type=="S2"){ # X is {q1,m,q3}
    sum_stat_index=c(floor(n/4)+1,floor(n/2)+1,floor(3*n/4)+1)
    nsum_stat=length(sum_stat_index)
    
    eta=2*stats::qnorm((0.75*n-0.125)/(n+0.25))
    sigmahat=(X[3]-X[1])/eta
  }
  
  if(type=="S3"){ # X is 5 # sum stat
    sum_stat_index=c(1,floor(n/4)+1,floor(n/2)+1,floor(3*n/4)+1,n)
    nsum_stat=length(sum_stat_index)
    
    tao=2*stats::qnorm((n-0.375)/(n+0.25))
    eta=2*stats::qnorm((0.75*n-0.125)/(n+0.25))
    sigmahat=0.5*(X[5]-X[1])/tao+0.5*(X[4]-X[2])/eta
  }
  
  out=list(sigmahat=sigmahat,
           sum_stat_index=sum_stat_index,sum_stat=X)
  return(out)
}


#' BLUEs of individual location and scale parameters
#'
#' To obtain the best linear unbiased estimator (BLUE) of location and scale parameters based on any set of order statistics (Yang et al., 2018), where the underlying distribution is assumed to be normal.
#' @param X a vector of ordered summary statistics
#' @param n the sample size
#' @param type a character string indicating which type of summary statistics is reported. The options for the \strong{type} argument are: 
#' \itemize{
#'    \item "S1" for the sample mean, minimum and maximum values
#'    \item "S2" for the sample mean, first and third quartiles
#'    \item "S3" for the sample mean, first and third quartiles, and minimum and maximum values
#'    \item "tertiles" for tertiles, "quintiles" for quintiles, and "deciles" for deciles
#' }
#' 
#' @examples
#' X<-c(1,4,10)
#' n<-30
#' type<-"S1"
#' BLUE_s(X,n,type)
#'
#' X<-c(5,8)
#' n<-45
#' type<-"tertiles"
#' BLUE_s(X,n,type)
#' @references
#' Yang X, Hutson AD, and Wang D. (2018). A generalized BLUE approach for combining location and scale information in a meta-analysis (Submitted).
#' @export
BLUE_s=function(X,n,type){   #weights: by Normal assumption
  
  if (type=="S1") sum_stat_index=c(1,floor(n/2)+1,n) else # X is {a,m,b}
    if (type=="S2") sum_stat_index=c(floor(n/4)+1,floor(n/2)+1,floor(3*n/4)+1) else # X is {q1,m,q3}
      if (type=="S3") sum_stat_index=c(1,floor(n/4)+1,floor(n/2)+1,floor(3*n/4)+1,n) else # X is 5 # sum stat
        if (type=="tertiles") sum_stat_index=c(floor(n/3)+1,floor(2*n/3)+1)
        if (type=="quintiles") sum_stat_index=c(floor(n/5)+1,floor(2*n/5)+1,floor(3*n/5)+1,floor(4*n/5)+1)
        if (type=="deciles") sum_stat_index=c(floor(n/10)+1,floor(2*n/10)+1,floor(3*n/10)+1,floor(4*n/10)+1,floor(5*n/10)+1,floor(6*n/10)+1,floor(7*n/10)+1,floor(8*n/10)+1,floor(9*n/10)+1)
        nsum_stat=length(sum_stat_index)
        
        InverseErfc<-function(u) (-1/sqrt(2))*stats::qnorm(u/2)
        Qr_d1=function(r) exp(InverseErfc(2*r/(n+1))^2)*sqrt(2*pi)
        Qr_d2=function(r) -2*sqrt(2)*exp(2*InverseErfc(2*r/(n+1))^2)*pi*InverseErfc(2*r/(n+1))
        Qr_d3=function(r) 2*sqrt(2)*exp(3*InverseErfc(2*r/(1+n))^2)*pi^(3/2)*(1+4*InverseErfc(2*r/(1+n))^2)
        Qr_d4=function(r) -4*sqrt(2)*exp(4*InverseErfc(2*r/(1+n))^2)*pi^2*InverseErfc(2*r/(1+n))*(7+12*InverseErfc(2*r/(1+n))^2)
        Ex_r<-function(rr) sapply(rr,function(rr) stats::qnorm(rr/(n+1))+Qr_d2(rr)*(rr/(n+1))*(1-(rr/(n+1)))/(2*(n+2))+
                                    ((rr/(n+1))*(1-(rr/(n+1)))/(n+2)^2)*(((1-(rr/(n+1)))-(rr/(n+1)))*Qr_d3(rr)/3+(rr/(n+1))*(1-(rr/(n+1)))*Qr_d4(rr)/8) )
        Var_r<-function(rr) sapply(rr,function(rr) (rr/(n+1))*(1-(rr/(n+1)))*Qr_d1(rr)^2/(n+2)+
                                     ((rr/(n+1))*(1-(rr/(n+1)))/(n+2)^2)*(2*((1-(rr/(n+1)))-(rr/(n+1)))*Qr_d1(rr)*Qr_d2(rr)+(rr/(n+1))*(1-(rr/(n+1)))*(Qr_d1(rr)*Qr_d3(rr)+Qr_d2(rr)^2/2)) )
        covar_rs<-function(rr,ss) ((rr/(n+1))*(1-(ss/(n+1)))*Qr_d1(rr)*Qr_d1(ss)/(n+2)+
                                     ((rr/(n+1))*(1-(ss/(n+1)))/(n+2)^2)*(((1-(rr/(n+1)))-(rr/(n+1)))*Qr_d2(rr)*Qr_d1(ss)+((1-(ss/(n+1)))-(ss/(n+1)))*Qr_d1(rr)*Qr_d2(ss)+
                                                                            (rr/(n+1))*(1-(rr/(n+1)))*Qr_d3(rr)*Qr_d1(ss)/2+(ss/(n+1))*(1-(ss/(n+1)))*Qr_d3(ss)*Qr_d1(rr)/2+(rr/(n+1))*(1-(ss/(n+1)))*Qr_d2(rr)*Qr_d2(ss)/2 ))
        
        alpha<-sapply(sum_stat_index,function(i) Ex_r(i))
        scovar<-matrix(nrow=nsum_stat,ncol=nsum_stat)
        ii=0
        for (rr in sum_stat_index){
          jj=0
          ii=ii+1
          for (ss in sum_stat_index) {
            jj=jj+1
            if (rr<ss) {scovar[ii,jj]=covar_rs(rr,ss)}
            else if (rr>ss) {scovar[ii,jj]=scovar[jj,ii]}
            else {scovar[ii,jj]=Var_r(rr)}
          }
        }
        
        ones<-rep(1,nsum_stat)
        A<-cbind(ones,alpha)
        omega<-solve(scovar)
        delta<-det(t(A)%*%omega%*%A)
        gamma<-omega%*%(ones%*%t(alpha)-alpha%*%t(ones))%*%omega/delta
        
        c_mu=-t(alpha)%*%gamma
        c_sigma=t(ones)%*%gamma
        
        muhat=c_mu%*%X
        sigmahat=c_sigma%*%X
        Var_mu=t(alpha)%*%omega%*%alpha*sigmahat^2/delta
        Var_sigma=t(ones)%*%omega%*%ones*sigmahat^2/delta
        
        out=list(alpha=alpha,B=scovar,
                 muhat=muhat,sigmahat=sigmahat,Var_mu=Var_mu,Var_sigma=Var_sigma)
        return(out)
}


#' BLUEs of global location and scale parameters
#'
#' To obtain the global or overall best linear unbiased estimator (BLUE) of location and scale parameters (Yang et al., 2018).
#' @param alpha_c the expectation of a combined standardized vector of ordered summary statistics, i.e. equation (3.21) in Yang et al. (2018).
#' @param B_c the variance-covariance matrix of a combined standardized vector of ordered summary statistics, i.e. equation (3.22) in Yang et al. (2018).
#' @param X_c a combined vector of ordered summary statistics.
#' @import Matrix
#' @examples
#' n1<-30 # sample sizes of three included studies
#' n2<-45
#' n3<-67
#' X1<-c(3,1.2) # the mean and standard deviation
#' X2<-c(1,4,10) # the sample mean, minimum and maximum values
#' X3<-c(1.5,3,5.5,8,12) # the sample mean, first and third quartiles, and minimum and maximum values
#' X_c<-c(X1[1],X2,X3)
#'
#' alpha1<-0  #Approximate by the CLT.
#' B1<-1/sqrt(n1)
#' alpha2<-BLUE_s(X2,n2,"S1")$alpha
#' B2<-BLUE_s(X2,n2,"S1")$B
#' alpha3<-BLUE_s(X3,n3,"S3")$alpha
#' B3<-BLUE_s(X3,n3,"S3")$B
#' 
#' alpha_c<-c(alpha1,alpha2,alpha3)
#' B_c<-Matrix::bdiag(B1,B2,B3)
#'
#' BLUE_c(alpha_c,B_c,X_c)
#'
#' @references
#' Yang X, Hutson AD, and Wang D. (2018). A generalized BLUE approach for combining location and scale information in a meta-analysis (Submitted).
#' @export
BLUE_c=function(alpha_c,B_c,X_c){
  N=length(X_c)
  ones_c<-rep(1,N)
  A_c<-cbind(ones_c,alpha_c)
  omega_c<-solve(B_c)
  delta_c<-det(t(A_c)%*%omega_c%*%A_c)
  gamma_c<-omega_c%*%(ones_c%*%t(alpha_c)-alpha_c%*%t(ones_c))%*%omega_c/delta_c
  
  c_mu=-t(alpha_c)%*%gamma_c
  c_sigma=t(ones_c)%*%gamma_c
  muhat=as.numeric(c_mu%*%X_c)
  sigmahat=as.numeric(c_sigma%*%X_c)
  Var_mu=as.numeric(t(alpha_c)%*%omega_c%*%alpha_c*sigmahat^2/delta_c)
  Var_sigma=as.numeric(t(ones_c)%*%omega_c%*%ones_c*sigmahat^2/delta_c)
  
  out=list(muhat=muhat,sigmahat=sigmahat,Var_mu=Var_mu,Var_sigma=Var_sigma)
  return(out)
}