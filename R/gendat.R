#' @title Generate Simulated Datasets for Cox Proportional Hazards Model
#'
#' @description
#' This function generates multiple datasets for survival analysis based on a Cox proportional hazards model.
#' The baseline hazard function follows either a Weibull or an exponential distribution, depending on the values of \code{lambda}.
#' The function ensures that the maximum observed time in both the control and treatment groups is checked for censoring.
#' If the maximum time is not censored, it is forced to be censored to maintain the desired censoring rate.
#'
#' @param dimension Integer. The number of datasets to be generated.
#' @param n Integer. The number of samples within each cluster.
#' @param K Integer. The number of clusters (groups) within each dataset.
#' @param lambda Numeric vector. A two-element vector specifying the parameters for the baseline distribution:
#' \itemize{
#'   \item If \code{lambda = c(a, b)}, where \code{a > 1}, the baseline follows a Weibull distribution.
#'   \item If \code{lambda = c(1, b)}, the baseline follows an exponential distribution.
#' }
#' @param b1 Vector. The regression coefficient for the covariates, affecting the hazard function. We suggest that the maximum of \code{b1} should be lower than 2.
#' @param theta Numeric. A parameter controlling the dependency structure between survival times within clusters.
#' Higher values indicate stronger within-cluster correlation.
#' @param censrate Numeric. The target censoring rate for the dataset.
#' @param type Character. If \code{type = 'bin'}, the covariates are generated as  binary variables; if \code{type = 'cont'}  continuous covariates are generated.
#'
#' @importFrom stats pnorm rbinom rnorm runif uniroot
#' @return A list containing:
#' \itemize{
#'   \item \code{data} - A list of data frames, each containing a generated dataset.
#'   \item \code{censoringrates} - A numeric vector representing the censoring rate for each dataset.
#'   \item \code{mean(censoringrates)} - The mean censoring rate across all datasets.
#' }
#'
#' @export
#' @examples
#' # Generate binary covariate datasets with 1 datasets, 10 clusters, and 6 samples per cluster
#' print(gendat(type = 'bin', dimension = 1, K = 6, n = 10, lambda = c(1, 2),
#'       b1 = c(log(2),-log(2)), theta = 8, censrate = 0.5))


gendat<-function(type='bin',dimension=10,K=30,n=2,lambda=c(1,2),b1=c(log(2),-0.1),theta=8,censrate=0.3){
datasets<-list()
censoringrates=rep(0,dimension)
for(lll in 1:dimension){
  dim=length(b1)
  tt=matrix(0,K,n)
  xxx=matrix(0,n*K,dim)
  for(i in 1:K){
    u=runif(n)
    xmat=matrix(0,dim,n)
    if(type=='bin'){
      for(j in 1:dim){
        xtemp=rbinom(n,1,0.5)
        xmat[j,]=xtemp
      }
    }
    else if(type=='cont'){
      for(j in 1:dim){
        xtemp=rnorm(n)
        xmat[j,]=xtemp
      }
    }
    else{stop('Invalid type of covariates')}

    #x2[i,]=x
    #x2[(1+(i-1)*dim):i*dim,]=xmat

    xxx[(1+(i-1)*n):(i*n),]=t(xmat)
    a=rep(0,n-1)
    tt[i,1]=(1/lambda[2])*(-exp(-b1%*%xmat[,1])*log(1-u[1]))^(1/lambda[1])
    for(m in 2:n)
    {
      a[m-1]=exp(theta*(lambda[2]*tt[i,m-1])^(lambda[1])*exp(b1%*%xmat[,m-1]))
      tt[i,m]=(1/lambda[2])*((theta^(-1))*log((m-1)-sum(a)+(sum(a)-(m-2))*
        (1-u[m])^(-(theta^(-1)+m-1)^(-1)))*exp(-b1%*%xmat[,m]))^((lambda[1])^(-1))
    }
  }

  for(i in 1:K){
    for(j in 1:dim){

    }
  }
  b1xxx=xxx %*% b1
  f2<-function(x1,y1)
  {
    sum((x1*lambda[2]*exp(b1xxx))^(-1)*(1-exp(-lambda[2]*exp(b1xxx)*x1)))-y1*K*n
  }
  cr<-uniroot(f2, c(0.01,10),y1=censrate)$root
  censor=runif(K*n,0,cr)
  c1=rep(0,K*n)
  t2=rep(0,K*n)
  t1=as.vector(t(tt))
  for(i in 1:(K*n))
  {
    if (t1[i]>censor[i])
    { c1[i]=0
    t2[i]=censor[i]
    }
    else
    { c1[i]=1
    t2[i]=t1[i]
    }
  }
  id=rep((1:K),each=n)
  cluster1 <- data.frame(cens = c1, xxx, id, time=t2)
  datasets[[lll]]<-cluster1
  censoringrates[lll]=sum(c1==0)/length(c1)
}
datasets_re=datasets
result_dde<-list(
  j3=datasets_re,
  j1=censoringrates,
  j2=mean(censoringrates)
)
names(result_dde)<-c('data','censoringrates','mean(censoringrates)')
ret=match.call()
gendat_result <- list(Call=ret,
  datasets = result_dde
)
return (gendat_result)
}
