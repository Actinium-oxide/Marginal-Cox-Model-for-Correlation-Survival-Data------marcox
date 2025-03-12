#' @title Generate Dummy Variables for marcox Analysis
#' @description
#' This function creates dummy variables for a specified categorical column, typically used within the \code{formula} parameter in the \code{marcox()} function.
#' If \code{d_v} is not provided, dummy variables will be automatically generated based on the unique values in the specified \code{typename} column.
#'
#' @param typename Character. The name of the categorical column for which dummy variables should be created. It must be a single string representing the column name.
#' @param d_v Optional. A vector specifying which categories should be converted into dummy variables.
#' Dummy variables will be generated for the values in \code{d_v}, while any values not included in this vector will be assigned an all-zero dummy variable.
#' If \code{d_v} is omitted, the function will automatically generate dummy variables based on all unique values in the \code{typename} column.
#' @param cluster22 temporary value.
#'
#' @return A list containing:
#' \itemize{
#'   \item A vector including the original \code{typename} and the names of the generated dummy variables.
#'   \item A data frame with the newly created dummy variables appended to the original dataset.
#' }
#'
#' @details
#' This function is primarily used to prepare categorical data for analysis within \code{marcox()}.
#' If the \code{d_v} parameter is omitted, the function automatically determines the categories from the \code{typename} column and creates dummy variables accordingly.
#' The resulting dummy variables will be added to the original data frame (\code{cluster2}).
#'
#' When providing the \code{d_v} parameter, it must explicitly specify which categories should be converted into dummy variables.
#' If categorical variables (e.g., \code{type}) are included, they should be passed as strings when using \code{factormar()}.
#'
#' Generally, binary covariates should not be converted into dummy variables using \code{factormar()}, as this may introduce unexpected errors.

factormar <- function(typename, d_v = NULL,cluster22) {
  col_num<-dim(cluster22)[2]
  #if (col_num != dim(cluster22)[2]) { cluster22 <- cluster2_backup }
  lp <- dim(cluster22)[2]
  lp_or <- lp
  if (is.null(d_v)) {
    dv_0 <- as.numeric(names(table(cluster22[, typename])))
    dv_1 <- dv_0[1]
    if (length(dv_0) > 2) {
      for (j in 2:(length(dv_0) - 1)) {
        dv_1 <- c(dv_1, dv_0[j])
      }
    }
    dv=dv_1
  } else {
    dv<-eval(str2lang(d_v))}

  typels=rep(0,length(dv))
  for (k in 1:length(dv)) {
    cluster22 <- cbind(cluster22, rep(0, dim(cluster22)[1]))
    str <- paste(typename, k, sep = '_')
    typels[k] <- str
  }
  lp <- lp + length(dv)
  diag_mat <- diag(length(dv))

  for (i in 1:dim(cluster22)[1]) {
    for (j in 1:length(dv)) {
      if (cluster22[i, typename] == dv[j]) {
        cluster22[i, (lp_or + 1):(lp_or + length(dv))] <- diag_mat[j, ]
      }
    }
  }
  colnames(cluster22)[(lp_or + 1):(lp_or + length(dv))]<-typels
  return(list(c(typename, colnames(cluster22[, (lp_or + 1):(lp_or + length(dv))])),
              cluster22[, (lp_or + 1):(lp_or + length(dv))]))
}
