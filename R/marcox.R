#' @title Macrox Analysis for Cox Proportional Hazards Models
#' @description
#' This function performs Macrox analysis for Cox proportional hazards models, incorporating clustered data
#' and handling time-dependent covariates. It estimates coefficients, standard errors, and p-values based on
#' the specified formula and dataset. In addition, an optional penalization functionality is provided via the
#' \code{penalize} parameter. When \code{penalize} is set to TRUE, the function applies a SCAD-type penalty with
#' tuning parameter selection via a GCV criterion to perform variable selection.
#'
#' @param formula A model formula that uses the \code{Surv()} function to define the survival outcome. It should include
#' both continuous and categorical covariates, where categorical variables must be specified using the \code{factor()} function.
#' @param dat A list containing the dataset as prepared by the \code{init()} function. This list must include the processed data frame
#' and original column names to ensure proper matching of variables.
#' @param penalize A logical flag indicating whether penalization should be applied for variable selection. If set to \code{TRUE},
#' the function will incorporate a SCAD-type penalty using a Generalized Cross Validation (GCV) criterion for tuning.
#' @param pen_fla A numeric value representing the penalty factor used in the SCAD penalty. The default value is typically around 3.7.
#' @param pen_tp An integer specifying the number of tuning parameter values to be evaluated in the penalization process.
#' @param weight_v A numeric value or vector used to adjust the effect of the penalty during the optimization.
#'
#' @return A data frame (if \code{penalize=FALSE}) or a list (if \code{penalize=TRUE}) containing the following components:
#' \itemize{
#'   \item \strong{For the standard Macrox analysis (\code{penalize=FALSE})}:
#'     \itemize{
#'       \item \code{coef} - The estimated regression coefficients.
#'       \item \code{exp(coef)} - The exponentiated coefficients (hazard ratios).
#'       \item \code{se(coef)} - The standard errors of the estimated coefficients.
#'       \item \code{z} - The z-statistics for testing the significance of the coefficients.
#'       \item \code{p} - The p-values associated with the coefficients.
#'     }
#'
#'   \item \strong{For the penalized analysis (\code{penalize=TRUE})}:
#'     \itemize{
#'       \item \code{result_pen} - A data frame with columns:
#'         \itemize{
#'           \item \code{pp} - The penalized regression coefficients.
#'           \item \code{var} - The variance estimates of the penalized coefficients.
#'           \item \code{CIlow} - The lower bounds of the 95% confidence intervals.
#'           \item \code{CIupp} - The upper bounds of the 95% confidence intervals.
#'         }
#'       \item \code{npp} - A numeric vector containing additional diagnostics including the correlation parameter (\code{rho}),
#'             the overdispersion factor (\code{pphi}), and the number of iterations (\code{KK1}).
#'       \item \code{GCV} - A matrix of Generalized Cross Validation (GCV) values corresponding to the evaluated tuning parameters.
#'     }
#' }
#'
#' @details
#' The \code{marcox()} function is specifically designed for survival data analysis using Cox proportional hazards models. It handles both clustered and time-dependent covariates effectively.
#' The survival outcome must be defined using the \code{Surv()} function in the model formula, and covariates can be included directly or by converting categorical variables with the \code{factor()} function.
#'
#' When the \code{penalize} parameter is enabled, the function applies a SCAD-type penalty to perform variable selection. A range of tuning parameters,
#' determined by \code{pen_tp}, is evaluated via a GCV criterion to select the optimal penalty strength. This process helps in reducing model complexity by
#' shrinking insignificant coefficients toward zero.
#'
#' @examples
#'   dat <- init(kidney_data, div = 2)
#'   formula <- Surv(time, cens) ~ sex + factor('type')
#'   # Standard Macrox analysis without penalization (default)
#'   result1 <- marcox(formula, dat)
#'
#'   # Penalized analysis for variable selection
#'   result2 <- marcox(formula, dat, penalize=TRUE, pen_fla=3.7, pen_tp=30, weight_v=1)
#'
#' @export
marcox<-function(formula,dat,penalize=FALSE,pen_fla=3.7,pen_tp=30,weight_v=NULL){
#preprocessing



############################
#  penalizing       part   #
############################

if(penalize) {
RR=1
  cluster2<<-dat[[1]]
  t2<<-cluster2[,as.character(formula[[2]][[2]])]
  c1<<-cluster2[,as.character(formula[[2]][[3]])]
  Y1<<-matrix(cluster2[,as.character(formula[[2]][[3]])],ncol=1)
  cens<<-c1
  t11<<-sort(t2)
  c11<<-c1[order(t2)]
  tt1<<-unique(t11[c11==1])
  kk<<-length(table(t11[c11==1]))
  dd<<-as.matrix(table(t11[c11==1]))
  gg1<<-rep(1,length(c1))
  g11<<-gg1[order(t2)]
  col_origin<-dat[[2]]
  col_origin_1<-''
  col_origin_1<-deparse(formula[[3]])
  col_origin_1=paste(col_origin_1,collapse = '')
  col_origin_1=gsub(' ','',col_origin_1)
  col_origin_1<<-gsub('[\r\n]+','',col_origin_1)
  index<<-strsplit(col_origin_1,'\\+')[[1]]
  cov_temp<<-c()
  xxx_1<<-c()
  xxx_2<<-c()
  xxx_21=c()
  typelist<<-c()
  typedumlist<<-c()
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
  xxx<<-as.matrix(xxx)
  dimcov<-dim(xxx)[2]
  n<<-length(id)/K
  x111=xxx[order(t2),]
  X1=xxx
  FLa=pen_fla
  xx1=x111[c11==1,]
  form_t=update(formula, . ~ . + cluster(id))
  coxre=coxph(form_t,cluster2)$coef
  betainit<<-coxre
  n_e<<-ncol(xxx)
  tplambda<-seq(0,1,length=pen_tp)
  GCV<-NULL
  for(ll in 1:length(tplambda)){
    cat('Evaluating Lambda:',ll*100/pen_tp,'%','\n')
    tlambda=tplambda[ll]
    # tlambda<-c(rep(0,5),rep(tplambda[ll],13),rep(0,5),rep(tplambda[ll],13))
    beta2=matrix(rep(0,n_e),ncol=1)

    gSSS1=rep(0,kk)
    KK1=1

    repeat{

      gSS=rep(0,kk)
      gSS1=rep(1,kk)

      gSS[1]=dd[1]/(sum(g11[min((1:(K*n))[t11==tt1[1]]):(K*n)]*exp(x111[(min((1:(K*n))[t11==tt1[1]]):(K*n)),]%*%betainit)))

      for (i in 1:(kk-1))
      {
        gSS[i+1]=gSS[i]+dd[i+1]/
          (sum(g11[min((1:(K*n))[t11==tt1[i+1]]):(K*n)]*exp(x111[min((1:(K*n))[t11==tt1[i+1]]):(K*n),]%*%betainit)))
      }

      gSS1=exp(-gSS)

      gSS2=rep(0,K*n)
      gSS3=rep(0,K*n)
      gSS4=rep(0,K*n)  #baseline hazard function for all survival time

      for(i in 1:(K*n))
      {  kk1=1

      if(t2[i]<tt1[1])
      {
        gSS2[i]=1
        gSS3[i]=0.00000001
        gSS4[i]=0.00000001
      }
      else {
        if(t2[i]>=tt1[kk])
        {
          # gSS2[i]=0
          gSS2[i]=0.000001
          gSS3[i]=gSS[kk]
          gSS4[i]=gSS[kk]-gSS[kk-1]
        }
        else {
          repeat{
            if(t2[i]>=tt1[kk1]) kk1=kk1+1
            else break
          }
          {      gSS2[i]=(gSS1[kk1-1])^(exp(xxx[i,]%*%betainit))
            gSS3[i]=gSS[kk1-1]
            if(kk1==2) gSS4[i]=gSS[kk1-1]
            else gSS4[i]=gSS[kk1-1]-gSS[kk1-2]
          }
        }
      }
      }

      Lambda<-gSS3
      W1<<-diag(Lambda)

      SK1=1
      beta1=matrix(rep(0,n_e),ncol=1)

      repeat{
        mu=exp(X1%*%betainit)
        newY1=c1/Lambda

        res=as.vector((newY1-mu)/sqrt(mu))
        rres=0

        pphi=(sum(res^2)/(K*n-dim(X1)[2]))
        #  pphi=1

        res=matrix(res,ncol=K)
        res=t(res)
        for(i in 1:K)
        {
          for(j in 1:(n-1))
            rres=rres+res[i,j]*sum(res[i,(j+1):n])
        }

        rho=(pphi^(-1))*rres/(K*n*(n-1)/2-dim(X1)[2])
        #  rho=0

        SK=1

        repeat{
          D1=diag(mu[id==1,])%*%diag(rep(1,n))%*%(X1[id==1,])
          for(i in 2:K)
          { D1=rbind(D1,diag(mu[id==i,])%*%diag(rep(1,n))%*%(X1[id==i,])) }

          S1=newY1-mu

          R1=matrix(rho,n,n)
          diag(R1)=1

          V1=sqrt(diag(mu[id==1,]))%*%R1%*%sqrt(diag(mu[id==1,]))*pphi
          for(i in 2:K)
          { V1=bdiag(V1,sqrt(diag(mu[id==i,]))%*%R1%*%sqrt(diag(mu[id==i,]))*pphi) }

          V1=as.matrix(V1)

          #  Z1=D1%*%betainit+S1
          #  geebeta=solve(t(D1)%*%solve(V1)%*%W1%*%D1)%*%t(D1)%*%solve(V1)%*%W1%*%Z1

          absbeta=abs(as.vector(betainit))
          qofbeta=tlambda*(absbeta<=tlambda)+(FLa*tlambda-absbeta)*(FLa*tlambda>=absbeta)*(absbeta>tlambda)/(FLa-1)

          Ebeta=diag(qofbeta/(0.000001+absbeta))

          geebeta=betainit+ginv(t(D1)%*%ginv(V1)%*%W1%*%D1+K*Ebeta)%*%(t(D1)%*%ginv(V1)%*%W1%*%S1-K*Ebeta%*%betainit)

          if(any(abs(geebeta-betainit)>1e-6) && (SK<=500))
          {
            betainit<<-geebeta
            mu=exp(X1%*%betainit)
            SK=SK+1
          }
          else break

        }



        if(any(abs(betainit-beta1)>0.000001) && (SK1<30))

        {  beta1=betainit
        # mu=exp(X1%*%betainit)
        SK1=SK1+1
        }

        else break
      }


      if (any(abs(betainit-beta2)>0.000001) || any(abs(gSS1-gSSS1)>0.000001) )
      {
        beta2<-betainit
        gSSS1<-gSS1
        # Lambda<-gSS3
        KK1<-KK1+1
      }


      else  break
    }


    sumofbeta=diag(qofbeta/absbeta)
    eoflambda=sum(diag(solve(-t(D1)%*%solve(V1)%*%W1%*%D1+sumofbeta)%*%(-t(D1)%*%solve(V1)%*%W1%*%D1)))

    gS=c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])
    # lofbeta=sum(log(gS*exp(xx1%*%betainit)))-sum(gSS3*(xxx%*%betainit))
    # lofbeta=sum(log(rep(gS,dd[,1])*exp(xx1%*%betainit)))-sum(gSS3*(xxx%*%betainit))
    # lofbeta=(sum(log(rep(gS,dd[,1])*exp(xx1%*%betainit)))-sum(gSS3*exp(xxx%*%betainit)))
    lofbeta=(sum(xx1%*%betainit)-sum(gSS3*exp(xxx%*%betainit)))

    GCV[ll]=-lofbeta/(K*(1-eoflambda/K)^2)

  }
  GCV1<-GCV[1:10]
  #print((1:10)[GCV1==min(GCV1)])
  goodlambda<-(tplambda[1:10])[GCV1==min(GCV1)]
  # tlambda<-goodlambda*sqrt(sandvnopenalty)
  # tlambda<-goodlambda
  tlambda<-goodlambda*weight_v
  betainit<<-coxre

  KK1=1
  gSSS1=rep(0,kk)
  # Lambda=t2
  beta2=rep(0,dimcov)
  countpp=1

  repeat{
    countpp<-countpp+1
    gSS=rep(0,kk)
    gSS1=rep(1,kk)

    gSS[1]=dd[1]/(sum(g11[min((1:(K*n))[t11==tt1[1]]):(K*n)]*exp(x111[(min((1:(K*n))[t11==tt1[1]]):(K*n)),]%*%betainit)))

    for (i in 1:(kk-1))
    {
      gSS[i+1]=gSS[i]+dd[i+1]/
        (sum(g11[min((1:(K*n))[t11==tt1[i+1]]):(K*n)]*exp(x111[min((1:(K*n))[t11==tt1[i+1]]):(K*n),]%*%betainit)))
    }

    gSS1=exp(-gSS)

    gSS2=rep(0,K*n)
    gSS3=rep(0,K*n)
    gSS4=rep(0,K*n)  #baseline hazard function for all survival time

    for(i in 1:(K*n))
    {  kk1=1

    if(t2[i]<tt1[1])
    {
      gSS2[i]=1
      gSS3[i]=0.00000001
      gSS4[i]=0.00000001
    }
    else {
      if(t2[i]>=tt1[kk])
      {
        # gSS2[i]=0
        gSS2[i]=0.000001
        gSS3[i]=gSS[kk]
        gSS4[i]=gSS[kk]-gSS[kk-1]
      }
      else {
        repeat{
          if(t2[i]>=tt1[kk1]) kk1=kk1+1
          else break
        }
        {      gSS2[i]=(gSS1[kk1-1])^(exp(xxx[i,]%*%betainit))
          gSS3[i]=gSS[kk1-1]
          if(kk1==2) gSS4[i]=gSS[kk1-1]
          else gSS4[i]=gSS[kk1-1]-gSS[kk1-2]
        }
      }
    }
    }



    Lambda<-gSS3
    W1<<-diag(Lambda)

    SK1=1
    beta1=matrix(rep(0,n_e),ncol=1)

    repeat{
      mu=exp(X1%*%betainit)
      newY1=c1/Lambda

      res=as.vector((newY1-mu)/sqrt(mu))
      rres=0

      pphi=(sum(res^2)/(K*n-dim(X1)[2]))
      #  pphi=1

      res=matrix(res,ncol=K)
      res=t(res)
      for(i in 1:K)
      {
        for(j in 1:(n-1))
          rres=rres+res[i,j]*sum(res[i,(j+1):n])
      }

      rho=(pphi^(-1))*rres/(K*n*(n-1)/2-dim(X1)[2])
      #  rho=0

      SK=1

      repeat{

        D1=diag(mu[id==1,])%*%diag(rep(1,n))%*%(X1[id==1,])
        for(i in 2:K)
        { D1=rbind(D1,diag(mu[id==i,])%*%diag(rep(1,n))%*%(X1[id==i,])) }

        S1=newY1-mu

        R1=matrix(rho,n,n)
        diag(R1)=1

        V1=sqrt(diag(mu[id==1,]))%*%R1%*%sqrt(diag(mu[id==1,]))*pphi
        for(i in 2:K)
        { V1=bdiag(V1,sqrt(diag(mu[id==i,]))%*%R1%*%sqrt(diag(mu[id==i,]))*pphi) }

        V1=as.matrix(V1)

        #  Z1=D1%*%betainit+S1
        #  geebeta=solve(t(D1)%*%solve(V1)%*%W1%*%D1)%*%t(D1)%*%solve(V1)%*%W1%*%Z1

        absbeta=abs(as.vector(betainit))
        qofbeta=tlambda*(absbeta<=tlambda)+(FLa*tlambda-absbeta)*(FLa*tlambda>=absbeta)*(absbeta>tlambda)/(FLa-1)

        Ebeta=diag(qofbeta/(0.000001+absbeta))

        geebeta=betainit+solve(t(D1)%*%solve(V1)%*%W1%*%D1+K*Ebeta)%*%(t(D1)%*%solve(V1)%*%W1%*%S1-K*Ebeta%*%betainit)

        if(any(abs(geebeta-betainit)>1e-6) && (SK<=500))
        {
          betainit<<-geebeta
          mu=exp(X1%*%betainit)
          SK=SK+1
        }
        else break

      }



      if(any(abs(betainit-beta1)>0.000001) && (SK1<30))

      {  beta1=betainit
      # mu=exp(X1%*%betainit)
      SK1=SK1+1
      }

      else break
    }


    if (any(abs(betainit-beta2)>0.000001) || any(abs(gSS1-gSSS1)>0.000001) )
    {
      beta2<-betainit
      gSSS1<-gSS1
      # Lambda<-gSS3
      KK1<-KK1+1
    }


    else  break
    if (countpp>9999){
      break
    }
  }


  ############################
  ##    variance estimate   ##
  ############################

  ###
  #  for beta variance
  ###

  adjustedbeta<-betainit
  adjustedbeta[abs(betainit)<=0.001]=0
  adjustedmu=exp(X1%*%adjustedbeta)
  nonzerobeta<-(1:length(adjustedbeta))[adjustedbeta!=0]

  betacorr=rho
  betascale=1
  # betacorr=0

  be=adjustedbeta
  gS=c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])

  gg1=rep(1,K*n)
  xxxx=xxx

  Q1=matrix(betacorr,n,n)
  diag(Q1)=1
  IQ1=solve(Q1)

  z2=xxx[,nonzerobeta]

  c2=c1

  mu2=exp(xxx%*%adjustedbeta)

  B2=matrix(0,n,n)
  ABC1=rep(0,K)
  VA1=matrix(0,dim(z2)[2],dim(z2)[2])

  for(v in 1:length(nonzerobeta))
  {
    for(w in 1:length(nonzerobeta))
    {
      for(i in 1:K)
      {
        z22=z2[id==i,]

        A2=t(z22[,v])

        c22=c2[id==i]
        Lam22=Lambda[id==i]
        # t22=t2[id==i]

        mu22=mu2[id==i,]

        BB1=(mu22^(1/2))%*%((t(mu22))^(-1/2))*IQ1

        for(s in 1:n)
        {
          for(l in 1:n)
          {
            B2[s,l]=(1/2)*(z22[s,w]-z22[l,w])*BB1[s,l]
          }
        }

        C2=(c22/Lam22)-mu22
        D2=BB1
        E2=z22[,w]*mu22
        G2=diag(Lam22)
        ABC1[i]=A2%*%(B2%*%G2%*%C2-D2%*%G2%*%E2)

      }
      VA1[v,w]=sum(ABC1)*(betascale^(-1))

      ABC1=rep(0,K)
    }

  }

  sdm=VA1
  newEbeta=diag((qofbeta/(0.000001+absbeta))[nonzerobeta])
  newsdm=sdm+K*newEbeta

  BBC=matrix(0,kk,length(nonzerobeta))

  for(s in 1:(kk))
  {
    BCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)

    for(j in 1:length(nonzerobeta))
    {
      BBC[s,j]=sum(exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)*(exp(-BCm)+BCm*exp(-BCm)-1)/((1-exp(-BCm))^2)*xxxx[(c1==1)&(t2==tt1[s]),j])+
        sum(gg1[t2>=tt1[s]]*exp(xxxx[t2>=tt1[s],]%*%be)*xxxx[t2>=tt1[s],j])
    }
  }

  CCC=rep(0,(kk))

  for(s in 1:(kk))
  {
    CCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)
    CCC[s]=sum(exp(2*(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)-CCm)/(1-exp(-CCm))^2)
  }

  BC=matrix(0,length(nonzerobeta),kk)

  for(r in 1:length(nonzerobeta))
  {

    for(s in 1:(kk))
    {
      elem=0
      for(i in 1:K)
      {
        mu22=mu2[id==i,]
        xxx1=xxx[id==i,r]
        t21=t2[id==i]

        for(j in 1:n)
        {
          if(t21[j]>=tt1[s])
            elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*IQ1[,j])*mu22[j]*(betascale^(-1))
        }
      }
      BC[r,s]=elem
    }
  }

  M22=-newsdm
  M23=BC
  M32=BBC
  M33=diag(CCC)

  M=rbind(cbind(M22,M23),cbind(M32,M33))

  ###### first derivative vector (fdv) #####

  fdm=0

  for(i in 1:K)
  {
    xxx1=xxxx[id==i,]

    gg11=gg1[id==i]
    c111=c1[id==i]
    t21=t2[id==i]

    g11i=Lambda[id==i]

    z22=z2[id==i,]
    mu22=mu2[id==i,]
    mu22m=diag(mu22)

    G2=diag(g11i)
    c22=c2[id==i]
    C2=(c22/g11i)-mu22

    fdv=t(mu22m%*%z22)%*%solve(sqrt(mu22m)%*%Q1%*%sqrt(mu22m)*(betascale))%*%G2%*%C2

    d11=rep(0,(kk))

    for(s in 1:(kk))
    {
      d11[s]=sum(exp(xxx1[(c111==1)&(t21==tt1[s]),]%*%be)/(1-exp(-gS[s]*exp(xxx1[(c111==1)&(t21==tt1[s]),]%*%be))))-
        sum(gg11[t21>=tt1[s]]*exp(xxx1[t21>=tt1[s],]%*%be))
    }


    fdm=fdm+t(t(c(fdv,d11)))%*%t(c(fdv,d11))

  }


  vcmR=solve(M)%*%fdm%*%t(solve(M))
  V1=diag(vcmR)[1:length(nonzerobeta)]
  sandv=V1

    result_pen<<-NULL
    result_pen<<-cbind(adjustedbeta,sandv,adjustedbeta[nonzerobeta]-1.96*sqrt(sandv),adjustedbeta[nonzerobeta]+1.96*sqrt(sandv))
    colnames(result_pen)<-c('pp','var','CIlow','CIupp')
    npp<<-c(rho,pphi,KK1)
    GCV_re<<-matrix(0,1,pen_tp)
    GCV_re<<-GCV
    cat('Call:\n')
    print(match.call())
    print(result_pen)
    print('npp')
    print(npp)
    print('GCV')
    print(GCV_re)

}



  else{
    if(T){
      n<<- rep(0,length(new_uid))
      for(i in 1:length(new_uid)){
        n[i]<<- sum(id==i)
      }
      cluster2<<-dat[[1]]
      t2<<-cluster2[,as.character(formula[[2]][[2]])]
      c1<<-cluster2[,as.character(formula[[2]][[3]])]
      Y1<<-matrix(cluster2[,as.character(formula[[2]][[3]])],ncol=1)
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
      col_origin<-dat[[2]]
      col_origin_1<-''
      col_origin_1<-deparse(formula[[3]])
      col_origin_1=paste(col_origin_1,collapse = '')
      col_origin_1=gsub(' ','',col_origin_1)
      col_origin_1<<-gsub('[\r\n]+','',col_origin_1)
      index<<-strsplit(col_origin_1,'\\+')[[1]]
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
      xxx<<-as.matrix(xxx)
      beta2<<-matrix(rep(0,dim(xxx)[2]),ncol=1)
      gSSS1<<-rep(0,kk)
      KK1<<-1
      x111<<-as.matrix(xxx[order(t2),])
      rownames(x111)<<-rep(1:dim(x111)[1])
      xx1<<-x111[c11==1,]

      X1<<-xxx
      Kn_ls<<-1:Kn
    }
  n<<- rep(0,length(new_uid))
  for(i in 1:length(new_uid)){
    n[i]<<- sum(id==i)
  }
  betainit<<-matrix(rep(0,dim(xxx)[2]),ncol=1)



  ############################C++ Reserve############################
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
          geebeta<<-solve(t(D1)%*%solve(V1)%*%W1%*%D1)%*%t(D1)%*%solve(V1)%*%W1%*%Z1
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
          SK1<<-SK1+1
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
  #############################C++ Reserve#########################
  if (FALSE){
    interface=kernal(dd_c=dd, g11_c=g11, t11_c=t11, tt1_c=tt1,
                     t2_c=t2, betainit_h=betainit, x111_c=x111, xxx_c=xxx,
                     c1_c=c1, n_c=n, id_c=id,
                     Kn_ls_c=Kn_ls, kk_c=kk, Kn_c=Kn, K_c=K)
    betainit<<-interface[1:length(betainit)]
    rho<<-interface[length(betainit)+1]
  }
  #############################C++ Reserve#########################
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
  else{
    for (i in 1:length(typedumlist)){
      rownames(result)[i]<<-typedumlist[i]
    }
  }
  cat('Call:\n')
  print(match.call())
  print(result)

  }
  }

