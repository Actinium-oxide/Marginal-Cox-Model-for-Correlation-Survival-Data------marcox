#' @title Analysis for Cox Proportional Hazards Models
#' @description
#' This function performs marcox analysis for Cox proportional hazards models, incorporating clustered data
#' and handling time-dependent covariates. It estimates coefficients, standard errors, and p-values based on
#' the specified formula and dataset. In addition, an optional penalization functionality is provided via the
#' \code{penalize} parameter. When \code{penalize} is set to TRUE, the function applies a SCAD-type penalty with
#' tuning parameter selection via a GCV criterion to perform variable selection.
#'
#' @param formula A model formula that uses the \code{Surv()} function to define the survival outcome. It should include
#' both continuous and categorical covariates, where categorical variables must be specified using the \code{factormar()} function.
#' @param dat A list containing the dataset as prepared by the \code{init()} function. This list must include the processed data frame
#' and original column names to ensure proper matching of variables.
#' @useDynLib marcox, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' @import Matrix
#' @import survival
#' @return A data frame (if \code{penalize=FALSE}) or a list (if \code{penalize=TRUE}) containing the following components:
#' \itemize{
#'   \item \strong{For the standard marcox analysis}:
#'     \itemize{
#'       \item \code{coef} - The estimated regression coefficients.
#'       \item \code{exp(coef)} - The exponentiated coefficients (hazard ratios).
#'       \item \code{se(coef)} - The standard errors of the estimated coefficients.
#'       \item \code{z} - The z-statistics for testing the significance of the coefficients.
#'       \item \code{p} - The p-values associated with the coefficients.
#'       \item \code{Correlation]} - The correlation coefficient of the data.
#'     }
#'
#'
#' }
#'
#' @details
#' The \code{marcox()} function is specifically designed for survival data analysis using Cox proportional hazards models. It handles both clustered and time-dependent covariates effectively.
#' The survival outcome must be defined using the \code{Surv()} function in the model formula, and covariates can be included directly or by converting categorical variables with the \code{factormar()} function.
#'
#' When the \code{penalize} parameter is enabled, the function applies a SCAD-type penalty to perform variable selection. A range of tuning parameters,
#' determined by \code{pen_tp}, is evaluated via a GCV criterion to select the optimal penalty strength. This process helps in reducing model complexity by
#' shrinking insignificant coefficients toward zero.
#'
#' @examples
#'   dat <- init(kidney_data, div = 2)
#'   formula <- Surv(time, cens) ~ sex + factormar('type', d_v=c(1,2,3))
#'   result1 <- marcox(formula, dat)
#' @export
marcox<-function(formula,dat){
  cluster2<-dat[[1]]
  col_num<-dim(cluster2)[2]
  new_id<-cluster2$id
  id<-new_id
  new_uid<-sort(unique(new_id))
  K<-length(unique(id))
  Kn<- dim(cluster2)[1]
      n<- rep(0,K)
      for(i in 1:K){
        n[i]<- sum(id==new_uid[i])
      }
      t2<-cluster2[,as.character(formula[[2]][[2]])]
      c1<-cluster2[,as.character(formula[[2]][[3]])]
      Y1<-matrix(cluster2[,as.character(formula[[2]][[3]])],ncol=1)
      cens<-c1
      t11<-sort(t2)
      c11<-c1[order(t2)]
      tt1<-unique(t11[c11==1])
      kk<-length(table(t11[c11==1]))
      dd<-as.matrix(table(t11[c11==1]))
      gg1<-rep(1,length(c1))
      g11<-gg1[order(t2)]
      typelist<-c()
      typedumlist<-c()
      col_origin<-dat[[2]]
      col_origin_1<-''
      col_origin_1<-deparse(formula[[3]])
      col_origin_1=paste(col_origin_1,collapse = '')
      col_origin_1=gsub(' ','',col_origin_1)
      col_origin_1<-gsub('[\r\n]+','',col_origin_1)
      index<-strsplit(col_origin_1,'\\+')[[1]]
      cov_temp<-c()
      xxx_1<-c()
      xxx_21=matrix(0,Kn,1)
      for (i in 1:length(index)){index[i]<-gsub(' ','',index[i])}
      for (j in index){
        if (grepl('factormar\\(',j)==FALSE){cov_temp<-c(cov_temp,j)}
        else{
          para=gsub(' ','',j)
          para=gsub('factormar\\(','',para)
          para=gsub('typename=','',para)
          para=gsub('d_v=','',para)
          if(grepl('c\\(',para)){
            d_vec <- substring(para,first = regexpr('c\\(',para)[[1]],last = nchar(para)-1)
            typenm <- substring(para,first=2,last =  regexpr('c\\(',para)[[1]]-3)
            factorls <- list(typename=typenm,d_v=d_vec,cluster22=cluster2)
          }else{
            para=substring(para,first = 2,last = nchar(para)-2)
            factorls=list(typename=para,cluster22=cluster2)
          }
          tm=do.call(factormar,factorls)
          #tm=eval(parse(text=j))
          tm_1=tm[[1]]
          typelist<-c(typelist,tm_1[1])
          typedumlist<-c(typedumlist,tm_1[2:length(tm_1)])

          tm_2 <- do.call(factormar,factorls)[[2]]
          xxx_21<-cbind(xxx_21,tm_2)
        }}
      xxx_21=xxx_21[,-1]
      # for (i in index){
      #   if(grepl('factormar',i)==TRUE){
      #     tm_2=eval(parse(text=i))[[2]]
      #     xxx_21<-c(xxx_21,tm_2)
      #   }
      # }

        xxx_2<-xxx_21

      if(is.null(cov_temp)==FALSE){
        xxx_1 <- as.matrix(cluster2[,cov_temp[1]])

        if (length(cov_temp)>=2){
          for (i in 2:length(cov_temp)){
          xxx_1<-cbind(xxx_1,cluster2[,cov_temp[i]])
        }
      }
      }
      if (is.null(typelist)==FALSE){xxx<-cbind(xxx_1,xxx_2)}
      else {if(is.null(xxx_1)){xxx<-xxx_2}
        else{xxx<-xxx_1}}
      xxx<-as.matrix(xxx)
      covnum<-dim(xxx)[2]
      beta2<-matrix(rep(0,dim(xxx)[2]),ncol=1)
      gSSS1<-rep(0,kk)
      KK1<-1
      x111<-as.matrix(xxx[order(t2),])
      rownames(x111)<-rep(1:dim(x111)[1])
      xx1<-x111[c11==1,]
      SK1<-1
      X1<-xxx
      Kn_ls<-1:Kn

  betainit<-matrix(rep(0,dim(xxx)[2]),ncol=1)



    repeat{
      gSS1 <- rep(1,kk)
      temp1<-sum(g11[min((1:(Kn))[t11==tt1[1]]):(Kn)]*exp(x111[(min((1:(Kn))[t11==tt1[1]]):(Kn)),]%*%betainit))
      gSS <- rep(1,kk)
      gSS[1]<-dd[1]/temp1
      for (i in 1:(kk-1))
      {
        gSS[i+1]<-gSS[i]+dd[i+1]/(sum(g11[min((1:(Kn))[t11==tt1[i+1]]):(Kn)]*exp(x111[min((1:(Kn))[t11==tt1[i+1]]):(Kn),]%*%betainit)))
      }

      gSS1<-exp(-gSS)
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
      Lambda<-gSS3
      W1<-diag(Lambda)

      beta1=matrix(rep(0,dim(xxx)[2]),ncol=1)
      # repeat{
      #   mu<<-exp(X1%*%betainit)
      #   newY1<<-c1/Lambda
      #   res<<-as.vector((newY1-mu)/sqrt(mu))
      #   rres=0
      #   pphi<<-(sum(res^2)/(sum(n)-dim(X1)[2]))
      #   for(i in 1:K) {
      #     group_res <- res[id == new_uid[i]]
      #     ni <- length(group_res)
      #     if(ni == 1) {
      #       rres <- rres + group_res[1]
      #     } else {
      #       for(j in 1:(ni-1)) {
      #         rres <- rres + group_res[j] * sum(group_res[(j+1):ni])
      #       }
      #     }
      #   }
      #   rho<<-(pphi^(-1))*rres/(sum(n*(n-1))/2-dim(X1)[2])
      #   SK=1
      #
      #   repeat{
      #     D1<<-matrix(0,Kn,ncol(xxx))
      #     temp2=1
      #     for(i in 1:K){
      #       id_eq=which(id==new_uid[i])
      #       mu_gp=as.matrix(mu[id_eq])
      #       D1[temp2:(temp2+n[i]-1),]<<-diag(as.vector(mu_gp)+1e-6)%*%diag(rep(1,n[i]))%*%(X1[id==new_uid[i],])
      #       temp2=temp2+n[i]
      #     }
      #
      #     S1<<-newY1-mu
      #     V1<<-matrix(0,Kn,Kn)
      #     temp3=1
      #     for(i in 1:K){
      #       id_eq=which(id==new_uid[i])
      #       mu_gp=as.vector(mu[id_eq])
      #       R1<-matrix(rho,n[i],n[i])
      #       diag(R1)<-1
      #       V1[temp3:(temp3+n[i]-1),temp3:(temp3+n[i]-1)]<-sqrt(diag(pmax(mu_gp, 1e-8)))%*%R1%*%sqrt(diag(pmax(mu_gp, 1e-8)))*pphi
      #       temp3=temp3+n[i]
      #     }
      #     V1<<-as.matrix(V1)
      #     Z1<<-D1%*%betainit+S1
      #
      #     sol_V1<<-as.matrix(matSolCpp(V1, block=TRUE, K, as.integer(n), tol_cond=1e+3))
      #
      #     #sol_V1<<-mat_sol(V1)
      #
      #     t_D1=t(D1)
      #     Z1<<-D1%*%betainit+S1
      #     geebeta<<-as.matrix(matSolCpp(t_D1%*%sol_V1%*%W1%*%D1,block=F,K, as.integer(n), tol_cond=1e+3))%*%t_D1%*%sol_V1%*%W1%*%Z1
      #
      #     if(any(abs(geebeta-betainit)>1e-6) && (SK<=500))
      #     {
      #
      # betainit <<- geebeta
      #       mu<<-exp(X1%*%betainit)
      #       SK<<-SK+1
      #     }
      #     else break
      #   }
      #   if(any(abs(betainit-beta1)>0.000001) && (SK1<30))
      #   {
      #     beta1<<-betainit
      #     SK1<<-SK1+1
      #   }
      #   else break
      # }
      #######################





      res_iter <- marcox_iterCpp(
        X1, betainit, Lambda, c1, W1, id, new_uid, n,
        tol = 1e-6, maxIter = 30, maxInner = 500,
        pphi = 1.0, rho = 0.0
      )
      betainit <- as.matrix(res_iter$betainit)
      mu       <- as.matrix(res_iter$mu)
      rho      <- res_iter$rho
      pphi     <- res_iter$pphi
      SK1      <- res_iter$SK1


      if (any(abs(betainit-beta2)>0.000001) || any(abs(gSS1-gSSS1)>0.000001) )
      {
        beta2<-betainit
        gSSS1<-gSS1
        KK1<-KK1+1
      }
      else  break
    }

  cat('Estimation of Beta Has Completed.\n')




  # betacorr<-rho
  # betascale<-1
  # be<-betainit
  # gS<-c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])
  # xxxx<-xxx
  # z2<-xxx
  # c2<-c1
  # mu2<-exp(z2%*%betainit)
  # ABC1<-rep(0,K)
  # VA1<-matrix(0,covnum,covnum)
  # for(v in 1:(covnum))
  # {
  #   for(w in 1:(covnum))
  #   {
  #     for(i in 1:K)
  #     {
  #       Q1<-matrix(betacorr,n[i],n[i])
  #       diag(Q1)<-1
  #       IQ1<-solve(Q1)
  #       B2<-matrix(0,n[i],n[i])
  #       z22<-matrix(z2[id==new_uid[i],],nrow=n[i],)
  #       A2<-t(z22[,v])
  #       c22<-c2[id==new_uid[i]]
  #       Lam22<-Lambda[id==new_uid[i]]
  #       mu22<-mu2[id==new_uid[i]]
  #       BB1<-(mu22^(1/2))%*%((t(mu22))^(-1/2))*IQ1
  #       for(s in 1:n[i])
  #       {
  #         for(l in 1:n[i])
  #         {
  #           B2[s,l]<-(1/2)*(z22[s,w]-z22[l,w])*BB1[s,l]
  #         }
  #       }
  #       C2<-(c22/Lam22)-mu22
  #       D2<-BB1
  #       E2<-z22[,w]*mu22
  #       G2<-diag(Lam22)
  #       ABC1[i]<-A2%*%(B2%*%G2%*%C2-D2%*%G2%*%E2)
  #     }
  #     VA1[v,w]<-sum(ABC1)*(betascale^(-1))
  #     ABC1<-rep(0,K)
  #   }
  # }
  # sdm<-VA1
  # BBC<-matrix(0,kk,covnum)
  # for(j in 1:covnum)
  # {
  #   for(s in 1:(kk))
  #   {
  #     BCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)
  #     BBC[s,j]=sum(exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)*(exp(-BCm)+BCm*exp(-BCm)-1)/((1-exp(-BCm))^2)*xxxx[(c1==1)&(t2==tt1[s]),j])+
  #       sum(gg1[t2>=tt1[s]]*exp(xxxx[t2>=tt1[s],]%*%be)*xxxx[t2>=tt1[s],j])
  #   }
  # }
  # CCC<-rep(0,(kk))
  # for(s in 1:(kk))
  # {
  #   CCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)
  #   CCC[s]<-sum(exp(2*(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)-CCm)/(1-exp(-CCm))^2)
  # }
  # BC<-matrix(0,covnum,kk)
  # for(r in 1:covnum)
  # {
  #   for(s in 1:(kk))
  #   {
  #     elem=0
  #     for(i in 1:K)
  #     {
  #       mu22=mu2[id==new_uid[i]]
  #       xxx1=xxx[id==new_uid[i],r]
  #       t21=t2[id==new_uid[i]]
  #
  #       Q1=matrix(betacorr,n[i],n[i])
  #       diag(Q1)=1
  #       IQ1=solve(Q1)
  #
  #       for(j in 1:n[i])
  #       {
  #         if(t21[j]>=tt1[s])
  #           elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*IQ1[,j])*mu22[j]*(betascale^(-1))
  #       }
  #     }
  #     BC[r,s]<-elem
  #   }
  # }
  # M22=-sdm
  # M23=BC
  # M32=BBC
  # M33=diag(CCC)
  # M<<-rbind(cbind(M22,M23),cbind(M32,M33))
  # fdm<-0
  # for(i in 1:K)
  # {
  #   xxx1=xxxx[id==new_uid[i],]
  #   gg11=gg1[id==new_uid[i]]
  #   c111=c1[id==new_uid[i]]
  #   t21=t2[id==new_uid[i]]
  #   g11=Lambda[id==new_uid[i]]
  #   z22=z2[id==new_uid[i],]
  #   mu22=mu2[id==new_uid[i],]
  #   mu22m=diag(mu22,n[i],n[i])
  #   G2=diag(g11)
  #   c22=c2[id==new_uid[i]]
  #   C2=(c22/g11)-mu22
  #   Q1=matrix(betacorr,n[i],n[i])
  #   diag(Q1)=1
  #   sol_temp=as.matrix(matSolCpp((sqrt(mu22m)%*%Q1%*%sqrt(mu22m)*(betascale)),block=FALSE,K,n,1e+3))
  #   fdv<-t(mu22m%*%z22)%*%sol_temp%*%G2%*%C2
  #   t22=t2[id==new_uid[i]]
  #   eqalpha=rep(0,kk)
  #   for(j in 1:kk)
  #   {
  #     Aalpha=(1:n[i])[(t22==tt1[j])&(c111==1)]
  #     Balpha=(1:n[i])[t22>=tt1[j]]
  #     if(length(Balpha)==0)
  #     {  eqalpha[j]<-0  }
  #     if(length(Balpha)!=0 & length(Aalpha)==0)
  #     {  eqalpha[j]<- -sum(gg11[Balpha]*mu22[Balpha]) }
  #     else
  #       eqalpha[j]<-sum(mu22[Aalpha]/(1-exp(-gS[j]*mu22[Aalpha])))-sum(gg11[Balpha]*mu22[Balpha])
  #   }
  #   fdm<-fdm+t(t(c(fdv,eqalpha)))%*%t(c(fdv,eqalpha))
  # }
  # sol_M=as.matrix(matSolCpp(M,block=FALSE,K,n,1e+3))
  #
  # vcmR<-sol_M%*%fdm%*%t(sol_M)
  # V1<-diag(vcmR)[1:covnum]
  # #V2<-(diag(solve(M)))[1:dim(xxx)[2]]
  # sandv<<-V1
  #naivv<<-V2



  betainit <- as.numeric(betainit)
  gSS      <- as.numeric(gSS)
  kk       <- as.integer(kk)
  covnum   <- as.integer(covnum)
  K        <- as.integer(K)
  n       <- as.integer(n)
  new_uid <- as.integer(new_uid)
  id      <- as.integer(id)
  xxx <- as.matrix(xxx)
  mode(xxx) <- "numeric"
  c1     <- as.numeric(c1)
  t2     <- as.numeric(t2)
  tt1    <- as.numeric(tt1)
  gg1    <- as.numeric(gg1)
  Lambda <- as.numeric(Lambda)

  sandv <- sandwich_rcpp(rho, 1, betainit, gSS, kk, covnum, K, n, new_uid,
                             xxx, c1, t2, tt1, gg1, Lambda, id)
  z=betainit/sqrt(sandv)
  p_value = 2 * (1 - pnorm(abs(z)))
  result<-data.frame(
    x1=c(betainit),
    x2=c(exp(betainit)),
    x3=c(sqrt(sandv)),
    x4=c(z),
    x5=c(p_value))
  colnames(result)<-c('coef','exp(coef)','se(coef)','z','p')
  k=0
  if(is.null(cov_temp)==FALSE){
    rownames(result)[1:length(cov_temp)]<-cov_temp
    if(is.null(typelist)==FALSE){
        for (i in 1:length(typedumlist)){
          rownames(result)[i+length(cov_temp)]<-typedumlist[i]
        }
      }
  }
  else{
    for (i in 1:length(typedumlist)){
      rownames(result)[i]<-typedumlist[i]
    }
  }
  cat('Call:\n')
  print(match.call())
  print(result)
  cat('Correlation:', rho)

  }
