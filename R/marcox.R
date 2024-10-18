#' @title Macrox Analysis for Cox Proportional Hazards Models
#'
#' @description
#' This function performs Macrox analysis for Cox proportional hazards models, incorporating clustered data and handling time-dependent covariates.
#' It estimates coefficients, standard errors, and p-values based on the specified formula and dataset.
#' The \code{formula} parameter can include dummy variables created using the \code{factor()} function, and the \code{dat} parameter must be initialized using the \code{init()} function.
#'
#' @param formula A formula object specifying the model to be fitted. Must include a \code{Surv()} function to define the survival object.
#' Dummy variables can be created using the \code{factor()} function.
#' If any \code{type} variables (categorical variables) are included, they must be passed as strings when using \code{factor()}.
#' @param dat A list containing the dataset initialized by the \code{init()} function. This list should include the processed data frame and any original column names.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{coef} - The estimated coefficients for the covariates.
#'   \item \code{exp(coef)} - The exponentiated coefficients, representing the hazard ratios.
#'   \item \code{se(coef)} - The standard errors of the estimated coefficients.
#'   \item \code{z} - The z-values for the coefficients.
#'   \item \code{p} - The p-values associated with the coefficients.
#' }
#'
#' @export
#'
#' @details
#' The \code{marcox()} function is used to analyze survival data, specifically for clustered or time-dependent scenarios.
#' The \code{formula} must be constructed using the \code{Surv()} function to properly define the time-to-event data, and can include covariates specified directly or created via the \code{factor()} function.
#' The \code{dat} parameter should be prepared using the \code{init()} function, which sets up the dataset with necessary variables and identifiers.
#'
#' This function outputs the analysis results, including estimates, standard errors, z-values, and p-values, and provides a summary of the model fit.

#' @examples
#' # Assuming `init` has been used to prepare the dataset
#' # dat <- init(kidney_data, div = 2)
#' # formula <- Surv(time,cens)~sex+factor('type')
#' # marcox(formula, dat)
marcox<-function(formula,dat){
  cluster2<<-dat[[1]]
  t2<<-cluster2[,as.character(formula[[2]][[2]])]
  c1<<-cluster2[,as.character(formula[[2]][[3]])]
  Y1<<-matrix(cluster2[,as.character(formula[[2]][[3]])],ncol=1)
  new_id<<-cluster2[,'id']
  new_uid<<-sort(unique(new_id))
  Kn<<- dim(cluster2)[1]
  n<<- rep(0,length(new_uid))
  for(i in 1:length(new_uid)){
    n[i]<<- sum(id==i)
  }
  cens<<-c1
  t11<<-sort(t2)
  c11<<-c1[order(t2)]
  tt1<<-unique(t11[c11==1])
  kk<<-length(table(t11[c11==1]))
  dd<<-as.matrix(table(t11[c11==1]))
  gg1<<-rep(1,length(c1))
  g11<<-gg1[order(t2)]
  Y1<<-matrix(cluster2[,as.character(formula[[2]][[3]])],ncol=1)
  typelist<<-c()
  typedumlist<<-c()
  col_origin<<-dat[[2]]
  col_origin_1<<-''
  col_origin_1<<-formula[[3]]
  index<<-strsplit(deparse(col_origin_1),'\\+')[[1]]
  cov_temp<<-c()
  xxx_1<<-c()
  xxx_2<<-c()
  xxx_21=c()
  for (i in 1:length(index)){index[i]<<-gsub(' ','',index[i])}
  for (j in index){
    if (grepl('factor',j)==FALSE){cov_temp<<-c(cov_temp,j)}
    else{
      tm=eval(parse(text=j))
      tm_1=tm[[1]]
      typelist<<-c(typelist,tm_1[1])
      typedumlist<<-c(typedumlist,tm_1[2:length(tm_1)])
      }}
  for (i in index){
    if(grepl('factor',i)==TRUE){
    tm_2=eval(parse(text=i))[[2]]
    xxx_21<-c(xxx_21,tm_2)
  }
  }
  for (i in 1:length(xxx_21)){
    xxx_2<<-cbind(xxx_2,xxx_21[[i]])
  }
    if(is.null(cov_temp)==FALSE){xxx_1 <<- as.matrix(cluster2[,cov_temp[1]])

    if (length(cov_temp)>=2){
      for (i in 2:length(cov_temp)){
        xxx_1<<-cbind(xxx_1,cluster2[,cov_temp[i]])
      }
    }
    }
  if (is.null(typelist)==FALSE){xxx<-cbind(xxx_1,xxx_2)}
  else {if(is.null(xxx_1)){xxx<<-xxx_2}
    else{xxx<<-xxx_1}}
  betainit<<-matrix(rep(0,dim(xxx)[2]),ncol=1)
  beta2<<-matrix(rep(0,dim(xxx)[2]),ncol=1)
  gSSS1<<-rep(0,kk)
  KK1<<-1
  x111<<-as.matrix(xxx[order(t2),])
  rownames(x111)<<-rep(1:dim(x111)[1])
  xxx<<-as.matrix(xxx)
  X1<<-xxx
  Kn_ls<<-1:Kn
#######################################
  if(TRUE){
  repeat{
    gSS <<- rep(0,kk)
    gSS1 <<- rep(1,kk)
    gSS[1]<<-dd[1]/(sum(g11[min((1:(Kn))[t11==tt1[1]]):(Kn)]*exp(x111[(min((1:(Kn))[t11==tt1[1]]):(Kn)),]%*%betainit)))
    for (i in 1:(kk-1))
    {
      gSS[i+1]<<-gSS[i]+dd[i+1]/(sum(g11[min((1:(Kn))[t11==tt1[i+1]]):(Kn)]*exp(x111[min((1:(Kn))[t11==tt1[i+1]]):(Kn),]%*%betainit)))
    }
    gSS1<<-exp(-gSS)
    gSS2<-rep(0,Kn)
    gSS3<-rep(0,Kn)
    for(i in 1:(Kn))
    {  kk1=1
    if(t2[i]<tt1[1])
    {
      gSS2[i]=1
      gSS3[i]=0.00000001
    }
    else {
      if(t2[i]>=tt1[kk])
      {
        gSS2[i]=0
        gSS3[i]=gSS[kk]
      }
      else {
        repeat{
          if(t2[i]>=tt1[kk1]) {kk1=kk1+1}
          else break
        }
        { gSS2[i]=(gSS1[kk1-1])^(exp(xxx[i,]%*%betainit))
          gSS3[i]=gSS[kk1-1]
        }
      }
    }
    }
    Lambda<<-gSS3
    W1<<-diag(Lambda)
    SK1=1
    beta1=matrix(rep(0,dim(xxx)[2]),ncol=1)
    repeat{
      mu<<-exp(X1%*%betainit)
      newY1<<-c1/Lambda
      res<<-as.vector((newY1-mu)/sqrt(mu))
      rres=0
      pphi<<-(sum(res^2)/(sum(n)-dim(X1)[2]))
      resm=matrix(0,ncol=K,nrow=max(n))
      for(i in 1:K)
      {
        resm[1:n[i],i]=res[id==i]
      }
      res<<-t(resm)
      for(i in 1:K)
      {
        if(n[i]==1)
        { rres=rres+res[i,1]
        }
        else
        {
          for(j in 1:(n[i]-1)){
            rres=rres+res[i,j]*sum(res[i,(j+1):n[i]])
          }
        }
      }
      rho<<-(pphi^(-1))*rres/(sum(n*(n-1))/2-dim(X1)[2])
      SK=1
      repeat{
        D1<<-diag(mu[id==1])%*%diag(rep(1,n[1]))%*%(X1[id==1,])
        for(i in 2:K){
          D1<<-rbind(D1,diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(X1[id==i,])) }
        S1<<-newY1-mu
        R1<<-matrix(rho,n[1],n[1])
        diag(R1)<<-1
        V1<<-sqrt(diag(mu[id==1]))%*%R1%*%sqrt(diag(mu[id==1]))*pphi
        for(i in 2:K)
        {
          R1<<-matrix(rho,n[i],n[i])
          diag(R1)<<-1
          V1<<-bdiag(V1,sqrt(diag(mu[id==i]))%*%R1%*%sqrt(diag(mu[id==i]))*pphi)
        }
        V1<<-as.matrix(V1)
        Z1<<-D1%*%betainit+S1
        ##############################
        geebeta<<-solve(t(D1)%*%solve(V1)%*%W1%*%D1)%*%t(D1)%*%solve(V1)%*%W1%*%Z1
        #################################
        if(any(abs(geebeta-betainit)>1e-6) && (SK<=500))
        {
          betainit<<-geebeta
          mu<<-exp(X1%*%betainit)
          SK<<-SK+1
        }
        else break
      }
      if(any(abs(betainit-beta1)>0.000001) && (SK1<30))
      {
        beta1<-betainit
        SK1<-SK1+1
      }
      else break
    }
    if (any(abs(betainit-beta2)>0.000001) || any(abs(gSS1-gSSS1)>0.000001) )
    {
      beta2<<-betainit
      gSSS1<<-gSS1
      KK1<<-KK1+1
    }
    else  break
  }
  }
  #############################
  if (FALSE){
  interface=kernal(dd_c=dd, g11_c=g11, t11_c=t11, tt1_c=tt1,
                     t2_c=t2, betainit_h=betainit, x111_c=x111, xxx_c=xxx,
                     c1_c=c1, n_c=n, id_c=id,
                     Kn_ls_c=Kn_ls, kk_c=kk, Kn_c=Kn, K_c=K)
  betainit<<-interface[1:length(betainit)]
  rho<<-interface[length(betainit)+1]
  }
  #############################
  betacorr<<-rho
  betascale<<-1
  be<<-betainit
  gS<<-c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])
  gg1<<-rep(1,Kn)
  xxxx<<-xxx
  z2<<-cbind(xxx)
  c2<<-c1
  mu2<<-exp(z2%*%betainit)
  ABC1<<-rep(0,K)
  VA1<<-matrix(0,dim(z2)[2],dim(z2)[2])
  for(v in 1:(dim(z2)[2]))
  {
    for(w in 1:(dim(z2)[2]))
    {
      for(i in 1:K)
      {
        Q1<<-matrix(betacorr,n[i],n[i])
        diag(Q1)<<-1
        IQ1<<-solve(Q1)
        B2<<-matrix(0,n[i],n[i])
        z22<<-matrix(z2[id==i,],nrow=n[i],)
        A2<<-t(z22[,v])
        c22<<-c2[id==i]
        Lam22<<-Lambda[id==i]
        mu22<<-mu2[id==i]
        BB1<<-(mu22^(1/2))%*%((t(mu22))^(-1/2))*IQ1
        for(s in 1:n[i])
        {
          for(l in 1:n[i])
          {
            B2[s,l]<-(1/2)*(z22[s,w]-z22[l,w])*BB1[s,l]
          }
        }
        C2<<-(c22/Lam22)-mu22
        D2<<-BB1
        E2<<-z22[,w]*mu22
        G2<<-diag(Lam22)
        ABC1[i]<<-A2%*%(B2%*%G2%*%C2-D2%*%G2%*%E2)
      }
      VA1[v,w]<<-sum(ABC1)*(betascale^(-1))
      ABC1<<-rep(0,K)
    }
  }
  sdm<<-VA1
  BBC<<-matrix(0,kk,dim(xxxx)[2])
  for(j in 1:dim(z2)[2])
  {
    for(s in 1:(kk))
    {
      BCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)
      BBC[s,j]=sum(exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)*(exp(-BCm)+BCm*exp(-BCm)-1)/((1-exp(-BCm))^2)*xxxx[(c1==1)&(t2==tt1[s]),j])+
        sum(gg1[t2>=tt1[s]]*exp(xxxx[t2>=tt1[s],]%*%be)*xxxx[t2>=tt1[s],j])
    }
  }
  CCC<<-rep(0,(kk))
  for(s in 1:(kk))
  {
    CCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)
    CCC[s]<<-sum(exp(2*(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)-CCm)/(1-exp(-CCm))^2)
  }
  BC<<-matrix(0,dim(xxxx)[2],kk)
  for(r in 1:dim(xxxx)[2])
  {
    for(s in 1:(kk))
    {
      elem=0
      for(i in 1:K)
      {
        mu22=mu2[id==i]
        xxx1=xxx[id==i,r]
        t21=t2[id==i]

        Q1=matrix(betacorr,n[i],n[i])
        diag(Q1)=1
        IQ1=solve(Q1)

        for(j in 1:n[i])
        {
          if(t21[j]>=tt1[s])
          elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*IQ1[,j])*mu22[j]*(betascale^(-1))
        }
      }
      BC[r,s]<<-elem
    }
  }
  M22=-sdm
  M23=BC
  M32=BBC
  M33=diag(CCC)
  M<<-rbind(cbind(M22,M23),cbind(M32,M33))
  fdm<-0
  for(i in 1:K)
  {
    xxx1=xxxx[id==i,]
    gg11=gg1[id==i]
    c111=c1[id==i]
    t21=t2[id==i]
    g11=Lambda[id==i]
    z22=z2[id==i,]
    mu22=mu2[id==i,]
    mu22m=diag(mu22,n[i],n[i])
    G2=diag(g11)
    c22=c2[id==i]
    C2=(c22/g11)-mu22
    Q1=matrix(betacorr,n[i],n[i])
    diag(Q1)=1
    fdv<<-t(mu22m%*%z22)%*%solve(sqrt(mu22m)%*%Q1%*%sqrt(mu22m)*(betascale))%*%G2%*%C2
    t22=t2[id==i]
    eqalpha=rep(0,kk)
    for(j in 1:kk)
    {
      Aalpha=(1:n[i])[(t22==tt1[j])&(c111==1)]
      Balpha=(1:n[i])[t22>=tt1[j]]
      if(length(Balpha)==0)
      {  eqalpha[j]<-0  }
      if(length(Balpha)!=0 & length(Aalpha)==0)
      {  eqalpha[j]<- -sum(gg11[Balpha]*mu22[Balpha]) }
      else
        eqalpha[j]<-sum(mu22[Aalpha]/(1-exp(-gS[j]*mu22[Aalpha])))-sum(gg11[Balpha]*mu22[Balpha])
    }
    fdm<-fdm+t(t(c(fdv,eqalpha)))%*%t(c(fdv,eqalpha))
  }
  vcmR<-solve(M)%*%fdm%*%t(solve(M))
  V1<-diag(vcmR)[1:dim(xxx)[2]]
  V2<-(diag(solve(M)))[1:dim(xxx)[2]]
  sandv<<-V1
  naivv<<-V2
  z=betainit/sandv
  p_value<<-2*(1-pnorm(abs(z)))
  result<<-data.frame(
    x1=c(betainit),
    x2=c(exp(betainit)),
    x3=c(sqrt(sandv)),
    x4=c(betainit/sqrt(sandv)),
    x5=c(p_value))
  colnames(result)<<-c('coef','exp(coef)','se(coef)','z','p')
  k=0
  if(is.null(cov_temp)==FALSE){
    rownames(result)[1:length(cov_temp)]<<-cov_temp
    if(is.null(typelist)==FALSE){
        for (i in 1:length(typedumlist)){
          rownames(result)[i+length(cov_temp)]<<-typedumlist[i]
        }
      }
    }
  cat('Call:\n')
  print(match.call())
  print(result)
}
