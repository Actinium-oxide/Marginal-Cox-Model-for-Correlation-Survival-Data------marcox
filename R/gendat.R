#' @title Generate Simulated Datasets for Cox Proportional Hazards Model
#'
#' @description
#' This function generates multiple datasets for survival analysis based on a Cox proportional hazards model.
#' The baseline hazard function follows either a Weibull or an exponential distribution, depending on the values of \code{lambda}.
#' The function ensures that the maximum observed time in both the control and treatment groups is checked for censoring.
#' If the maximum time is not censored, it is forced to be censored to maintain the desired censoring rate.
#' Users can also specify \code{result_data_length} to control the number of rows in each dataset.
#' @param dimension Integer. The number of datasets to be generated.
#' @param K Integer. The number of samples within each cluster.
#' @param n Integer. The number of clusters (groups) within each dataset.
#' @param lambda Numeric vector. A two-element vector specifying the parameters for the baseline distribution:
#' \itemize{
#'   \item If \code{lambda = c(a, b)}, where \code{a > 1}, the baseline follows a Weibull distribution.
#'   \item If \code{lambda = c(1, b)}, the baseline follows an exponential distribution.
#' }
#' @param b1 Numeric. The regression coefficient for the covariate, affecting the hazard function.
#' @param theta Numeric. A parameter controlling the dependency structure between survival times within clusters.
#' Higher values indicate stronger within-cluster correlation.
#' @param censrate Numeric. The target censoring rate for the dataset.
#' @param binary Logical. If \code{TRUE}, the covariate is generated as a binary variable; otherwise, a continuous covariate is generated.
#' @param result_data_length Integer or \code{NULL}. If specified, truncates each dataset to the given number of rows.
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
#' # Generate datasets with 5 datasets, 2 clusters, and 100 samples per cluster
#' gendat(binary = FALSE, dimension = 5, K = 100, n = 2, lambda = c(1, 2),
#'       b1 = log(2), theta = 8, censrate = 0.3, result_data_length = 20)


gendat<-function(binary=T,dimension=10,K=30,n=2,lambda=c(1,2),b1=log(2),theta=5,censrate=0.5,result_data_length=NULL){
datasets<-list()
censoringrates=rep(0,dimension)
for(lll in 1:dimension){
  tt=matrix(0,K,n)
  x2=matrix(0,K,n)
  for(i in 1:K){
    u=runif(n)
    if(binary){x=rbinom(n,1,0.5)}else{x=rnorm(n)}
    x2[i,]=x
    a=rep(0,n-1)
    tt[i,1]=(1/lambda[2])*(-exp(-b1*x[1])*log(1-u[1]))^(1/lambda[1])
    for(m in 2:n)
    {
      a[m-1]=exp(theta*(lambda[2]*tt[i,m-1])^(lambda[1])*exp(b1*x[m-1]))
      tt[i,m]=(1/lambda[2])*((theta^(-1))*log((m-1)-sum(a)+(sum(a)-(m-2))*
        (1-u[m])^(-(theta^(-1)+m-1)^(-1)))*exp(-b1*x[m]))^((lambda[1])^(-1))
    }
  }
  xxx=as.vector(t(x2))
  f2<-function(x1,y1)
  {
    sum((x1*lambda[2]*exp(b1*xxx))^(-1)*(1-exp(-lambda[2]*exp(b1*xxx)*x1)))-y1*K*n
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
if (is.null(result_data_length)==FALSE){
for(i in 1:length(datasets)){
  datasets_re[[i]]=datasets[[i]][1:result_data_length,]
}}
result_dde<-list(
  j3=datasets_re,
  j1=censoringrates,
  j2=mean(censoringrates)
)
names(result_dde)<-c('data','censoringrates','mean(censoringrates)')
cat('Call:\n')
print(match.call())
print('lambda')
print(lambda)
print(paste('b1',b1,sep='='))
print(paste('theta',theta,sep='='))
print(result_dde)
}
