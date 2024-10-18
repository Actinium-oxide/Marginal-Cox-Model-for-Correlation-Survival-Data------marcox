#' @title Generate Dummy Variables for Macrox Analysis
#'
#' @description
#' This function generates dummy variables for specified categorical columns, typically used as part of the \code{formula} parameter within the \code{marcox()} function.
#' If \code{d_v} is not provided, dummy variables will be automatically created based on the unique values within the specified \code{typename} column.
#'
#' @param typename Character. The name of the column for which dummy variables are to be created. Should be a single string representing the column name.
#' @param d_v Optional. A vector of values specifying which categories should be used to generate dummy variables.
#' If \code{d_v} is not provided, the function will automatically generate dummy variables based on the elements present in the \code{typename} column.
#'
#' @return A list containing:
#' \itemize{
#'   \item A vector with the original \code{typename} and the names of the generated dummy variables.
#'   \item A data frame of the generated dummy variables added to the original data.
#' }
#'
#' @export
#'
#' @details
#' This function is typically used to prepare categorical data for analysis within \code{marcox()}.
#' If the \code{d_v} parameter is omitted, the function will automatically determine the categories from the \code{typename} column and create dummy variables.
#' The generated dummy variables will be appended to the original data frame (\code{cluster2}).
#' If the \code{d_v} parameter is provided, it should specify which categories to convert to dummy variables.
#' If any \code{type} variables (categorical variables) are included, they must be passed as strings when using \code{factor()}.
factor <- function(typename, d_v = NULL) {
  if (col_num != dim(cluster2)[2]) { cluster2 <<- cluster2_backup }
  typename = as.character(typename)
  lp <<- dim(cluster2)[2]
  lp_or <<- lp
  if (is.null(d_v)) {
    dv_0 <<- names(table(cluster2[, typename[1]]))
    dv_1 <<- dv_0[1]
    if (length(dv_0) > 2) {
      for (j in 2:(length(dv_0) - 1)) {
        dv_1 <<- c(dv_1, dv_0[j])
      }
    }
    dv <<- dv_1
  } else { dv <<- d_v }

  for (k in 1:length(dv)) {
    cluster2 <<- cbind(cluster2, rep(0, dim(cluster2)[1]))
    str <- paste(typename[1], k, sep = '_')
    colnames(cluster2)[k + lp] <<- str
  }

  lp <<- lp + length(dv)
  diag_mat <- diag(length(dv))

  for (i in 1:dim(cluster2)[1]) {
    for (j in 1:length(dv)) {
      if (cluster2[i, typename[1]] == dv[j]) {
        cluster2[i, (lp_or + 1):(lp_or + length(dv))] <<- diag_mat[j, ]
      }
    }
  }

  return(list(c(typename, colnames(cluster2[, (lp_or + 1):(lp_or + length(dv))])),
              cluster2[, (lp_or + 1):(lp_or + length(dv))]))
}
factor<-function(typename,d_v=NULL){
  if (col_num!=dim(cluster2)[2]){cluster2<<-cluster2_backup}
  typename=as.character(typename)
  lp<<-dim(cluster2)[2]
  lp_or<<-lp
  if(is.null(d_v)==TRUE){
    dv_0<<-names(table(cluster2[,typename[1]]))
    dv_1<<-dv_0[1]
    if(length(dv_0)>2){
    for (j in 2:(length(dv_0)-1)){
      dv_1<<-c(dv_1,dv_0[j])
    }}
    dv<<-dv_1
  }
  else{dv<<-d_v}
    for(k in 1:length(dv)){
      cluster2 <<- cbind(cluster2,rep(0,dim(cluster2)[1]))
      str<-paste(typename[1],k,sep='_')
      colnames(cluster2)[k+lp] <<- str
    }
    lp<<-lp+length(dv)
    diag_mat<-diag(length(dv))
    for (i in 1:dim(cluster2)[1]){
      for (j in 1:length(dv)){
        if (cluster2[i,typename[1]]== dv[j]){
          cluster2[i,(lp_or+1):(lp_or+length(dv))] <<- diag_mat[j,]
        }
      }
    }
  #return(colnames(cluster2[,(lp_or+1):(lp_or+length(dv))]))
    return(list(c(typename,colnames(cluster2[,(lp_or+1):(lp_or+length(dv))])),cluster2[,(lp_or+1):(lp_or+length(dv))]))
}
