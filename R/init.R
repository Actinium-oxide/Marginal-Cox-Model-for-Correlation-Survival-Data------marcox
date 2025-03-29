#' @title Initialize Data for Macrox Analysis
#' @description
#' This function initializes data for further analysis, typically used as a parameter within the \code{marcox} function.
#' It loads data from a specified file or data frame, performs basic preprocessing, and assigns unique identifiers for grouping.
#' The function allows for simple clustering based on the \code{div} parameter, but if the data has more complex observational structures,
#' users should preprocess the data before running this function.
#'
#' @param ad The file path or data frame to analyze. If a file path is provided, the file will be loaded into a data frame.
#' The file should be in a tabular format (e.g., .csv, .txt).
#' @param div Integer. The number of observation points per sample. If provided, the data will be divided accordingly. Default is \code{NULL}.
#' @param sep parameter. The \code{sep} parameter specifies the character that separates
#' the fields in each line of the file. For instance, for a comma-separated file, set \code{sep = ","},
#' and for a tab-separated file, set \code{sep = "\t"}.
#' @param col_id Character. The name of column that identifies the clusters.
#' If the data has complex observational situations, please preprocess the data before using this function.
#' @return A list containing:
#' \itemize{
#'   \item \code{cluster2} - The processed data frame after applying any transformations, including the addition of unique IDs.
#'   \item \code{dv} - A list of dummy variables. If no dummy variables are created, returns an empty list.
#'   \item \code{col_name_origin} - The original column names of the data frame.
#' }
#' @export
#' @importFrom utils read.table
#' @details
#' This function is generally used as a preparatory step for the \code{marcox} function.
#' It handles loading and preprocessing data by setting up unique identifiers (`id`) based on the number of observation points specified by \code{div}.
#' If \code{div} is provided, the function will ensure consistent grouping across the data.
#' If the number of rows in the data is not divisible by \code{div}, an error will be shown.
#' For more complex observational settings, it is recommended to preprocess the data before running this function.
#' @examples
#' # Use an existing data frame without specifying div
#' sample_data_1 <- kidney_data
#' init(sample_data_1,div=2)
init <- function(ad,sep='\t',col_id='id',div=NULL){
  if(typeof(ad)=='character'){
    cluster2 <- read.table(ad,header=T,sep=sep)

}
  else {cluster2<-ad}
  #col_name_origin=colnames(cluster2)


  if(is.null(div)==FALSE){
  if (dim(cluster2)[1]%%div!=0){ stop ('ERROR IN DIVISION')}

  else{
    id<-NULL
    l<-dim(cluster2)[1]
    col_n<-dim(cluster2)[2]+1
    cluster2[,col_id] <- rep(1:l,each=div,length.out=(l%/%div)*div)
    id<-cluster2[,col_id]
    cluster2[,'original_id'] <- id
    col_num<-dim(cluster2)[2]
    cluster2_backup<-cluster2
    uid <- sort(unique(id))
    newid <- rep(0,length(id))
    for(i in 1:length(id)){
      j<-1
      repeat{
        if (id[i]!=uid[j]){
          j<-j+1
        }
        else{
          newid[i] <- j
          break
        }
      }
    }
    cluster2[,col_id] <- newid


  }
  }
  else{
    id<-cluster2[,col_id]
    cluster2[,'original_id'] <- id
    col_num<-dim(cluster2)[2]
    cluster2_backup<-cluster2
  }

  cluster2<-cluster2[order(cluster2[,col_id]),]
  colnames(cluster2)[which(colnames(cluster2)==col_id)] <- 'id'
  return(cluster2)
}
