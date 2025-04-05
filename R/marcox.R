#' @title Analysis for Cox Proportional Hazards Models
#' @description
#' This function performs marcox analysis for Cox proportional hazards models, incorporating clustered data
#' and handling time-dependent covariates. It estimates coefficients, standard errors, and p-values based on
#' the specified formula and dataset.
#'
#' @param formula A model formula that uses the \code{Surv()} function to define the survival outcome. It should include
#' both continuous and categorical covariates, where categorical variables must be specified using the \code{factormar()} function.
#' @param data The file path or the dataset(matrix) to be analyzed. If a file path is provided, the file will be loaded into a matrix.
#' The file should be in a tabular format (e.g., .csv, .txt).
#' @param sep Character. The \code{sep} parameter specifies the character that separates
#' the fields in each line of the file. For instance, for a comma-separated file, set \code{sep = ","},
#' and for a tab-separated file, set \code{sep = "\t"}.
#' @param col_id Character. The name of column that identifies the clusters.
#' @param div Integer. The number of observation points per sample. If provided, the data will be divided accordingly.
#' If the data has complex observational situations, please preprocess the data before using this function.
#' @param method The method employed to solve the correlation coefficient:
#' \itemize{
#'     \item Exchangeable correlation structure: \code{method = 'exchangeable'}<<Default>>
#'     \item Autoregressive(AR-1): \code{method = 'ar1'}
#'     \item k-dependent: \code{method = 'kdependent'}
#'     \item Toeplitz: \code{method = 'toeplitz'}
#'     \item Independent: \code{method = 'independent'}
#'     \item Unstructured: \code{method = 'unstructured'}
#'     }
#' @param k_value The k value only for k-dependent structure. The default value is 1.
#' @param plot_x
#' @param x_axis
#' @param y_axis
#' @param size
#' @param diagnose Diagnose option.
#' @param iteration
#'
#' @useDynLib marcox, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppEigen
#' @import Matrix
#' @import survival
#' @import ggplot2
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{coef} - The estimated regression coefficients.
#'   \item \code{exp(coef)} - The exponentiated coefficients (hazard ratios).
#'   \item \code{se(coef)} - The standard errors of the estimated coefficients.
#'   \item \code{z} - The z-statistics for testing the significance of the coefficients.
#'   \item \code{p} - The p-values associated with the coefficients.
#'   \item (hidden).correlation - Correlation coefficients of the data.
#'    }
#' @details
#' The \code{marcox()} function is specifically designed for survival data analysis using Cox proportional hazards models. It handles both clustered and time-dependent covariates effectively.
#' The survival outcome must be defined using the \code{Surv()} function in the model formula, and covariates can be included directly or by converting categorical variables with the \code{factormar()} function.
#' @examples
#'   formula <- Surv(time, cens) ~ sex + factormar('type', d_v=c(1,2,3))
#'   marcox(formula, data = kidney_data, div = 2, method = 'exchangeable', plot = TRUE, plot_x = 'sex')
#' @export
#'
marcox<-function(formula,data,method='exchangeable',iteration='newton_raphson',sep=NULL,col_id='id',div=NULL,
                 k_value=1,plot_x=NULL,x_axis='Time',y_axis='Survival Rates',size=0.5,diagnose=FALSE){

    res <- marcox.fit(formula=formula,data=data,sep=sep,col_id=col_id,div=div,method=method,k_value=k_value,diagnose=diagnose,iteration=iteration)
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
    method <- res$mtd
    methodd <- res$mtdd
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
      .method = method,
      .correlation = correlation_val,
      .phi = phi,
      .cens1=cens1,
      .cens0=cens0,
      .loglik=loglik,
      .df=df,
      .p_value=ppp
    )

    class(marcox_result) <- "marcox"


if(is.null(plot_x)==FALSE){
    coef=marcox_result$Estimation[plot_x,'coef']
    km_0 <- survfit(update.formula(formula,~1), data = data[data[,plot_x] == 0,])
    base_surv <- data.frame(time = km_0$time, surv = km_0$surv)
    pred_surv_1 <- data.frame(
      time  = base_surv$time,
      surv  = base_surv$surv^(exp(coef)),
      group = "marcox:1"
    )
    pred_surv_0 <- data.frame(
      time  = base_surv$time,
      surv  = base_surv$surv,
      group = "marcox:0"
    )
    pred_data <- rbind(
      pred_surv_1,
      pred_surv_0
    )
    plot_res <- ggplot2::ggplot()+ggplot2::geom_line(data = pred_data, ggplot2::aes(x = time, y = surv, color = group),
                                         size = size, linetype = "solid")+
      ggplot2::labs(
        x = x_axis,
        y = y_axis
      ) +
      ggplot2::theme_minimal()
    return(plot_res)
}

    return(marcox_result)
}














