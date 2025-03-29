#' @title Analysis for Cox Proportional Hazards Models
#' @description
#' This function performs marcox analysis for Cox proportional hazards models, incorporating clustered data
#' and handling time-dependent covariates. It estimates coefficients, standard errors, and p-values based on
#' the specified formula and dataset.
#'
#' @param formula A model formula that uses the \code{Surv()} function to define the survival outcome. It should include
#' both continuous and categorical covariates, where categorical variables must be specified using the \code{factormar()} function.
#' @param dat A list containing the dataset as prepared by the \code{init()} function. This list must include the processed data frame
#' and original column names to ensure proper matching of variables.
#' @param method The method employed to solve the correlation coefficient:
#' \itemize{
#'     \item Exchangeable correlation structure: \code{method = 'exc'}<<Default>>
#'     \item Autoregressive(AR-1): \code{method = 'ar1'}
#'     \item k-dependent: \code{method = 'kd'}
#'     \item Toeplitz: \code{method = 'toep'}
#'     \item Independent: \code{method = 'indp'}
#'     \item Unstructured: \code{method = 'uns'}
#'     }
#' @param k_value The k value only for k-dependent structure. The default value is 1.
#'
#' @useDynLib marcox, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' @import Matrix
#' @import survival
#' @return A data frame containing the following components:
#' \itemize{
#'
#'       \item \code{coef} - The estimated regression coefficients.
#'       \item \code{exp(coef)} - The exponentiated coefficients (hazard ratios).
#'       \item \code{se(coef)} - The standard errors of the estimated coefficients.
#'       \item \code{z} - The z-statistics for testing the significance of the coefficients.
#'       \item \code{p} - The p-values associated with the coefficients.
#'       \item \code{correlation]} - correlation coefficients of the data.
#'
#' }
#' @details
#' The \code{marcox()} function is specifically designed for survival data analysis using Cox proportional hazards models. It handles both clustered and time-dependent covariates effectively.
#' The survival outcome must be defined using the \code{Surv()} function in the model formula, and covariates can be included directly or by converting categorical variables with the \code{factormar()} function.
#' @importFrom utils getFromNamespace
#' @examples
#'   dat <- init(kidney_data, div = 2)
#'   formula <- Surv(time, cens) ~ sex + factormar('type', d_v=c(1,2,3))
#'   result1 <- marcox(formula, dat, method = 'exc')
#' @export
#'
marcox<-function(formula,dat,method='exc',k_value=1){

    res <- marcox.fit(formula,dat,method,k_value)
    callmsg <- match.call()
    # list(
    #   coef=betainit,
    #   se=sandv,
    #   zvalue=z,
    #   pvalue=p_value,
    #   corr_exc=rho,
    #   corr_toeporkd=rho_vec_k,
    #   corr_uns=rhomat,
    #   cens1times=length(tt1),
    #   cens0times=Kn-length(tt1)
    # )
    betainit <- res$coef
    sandv <- res$se
    z <- res$zvalue
    p_value <- res$pvalue
    cens1 <- res$cens1times
    cens0 <- res$cens0times
    rho <- res$corr_exc
    rho_vec_k <- res$corr_toeporkd
    rhomat <- res$corr_uns
    cov_temp <- res$cov
    typedumlist <- res$dummylist
    methodd <- res$mtd
    typelist <- res$typelist
    phi=res$phi
    loglik=res$LR
    df=res$df
    ppp=res$ppp

    result <- data.frame(
      x1 = c(betainit),
      x2 = c(exp(betainit)),
      x3 = c(sqrt(sandv)),
      x4 = c(z),
      x5 = c(p_value)
    )
    colnames(result) <- c("coef", "exp(coef)", "se(coef)", "z", "p")

    if(!is.null(cov_temp)) {
      rownames(result)[1:length(cov_temp)] <- cov_temp
      if(!is.null(typelist)) {
        for(i in seq_along(typedumlist)){
          rownames(result)[i + length(cov_temp)] <- typedumlist[i]
        }
      }
    } else {
      for(i in seq_along(typedumlist)){
        rownames(result)[i] <- typedumlist[i]
      }
    }
    if(methodd %in% c(0,1,5)){
      correlation_val <- rho
    } else if(methodd %in% c(2,3)){
      correlation_val <- rho_vec_k
    } else if(methodd == 4){
      correlation_val <- rhomat
    }
    marcox_result <- list(
      call=callmsg,
      Estimation =  result,
      .correlation = correlation_val,
      .phi = phi,
      .cens1=cens1,
      .cens0=cens0,
      .loglik=loglik,
      .df=df,
      .p_value=ppp
    )

    class(marcox_result) <- "marcox"

    return(marcox_result)
}














