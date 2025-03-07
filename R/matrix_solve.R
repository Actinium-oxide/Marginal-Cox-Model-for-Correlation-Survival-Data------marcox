#' Title
#'
#' @param mat d
#' @param block f
#'
#' @returns
#' @export
#'
#' @examples
mat_sol <- function(mat,block=T) {
  n_dim <- ncol(mat)
  svd_res <- svd(mat)
  min_sv <- min(svd_res$d)
  max_sv <- max(svd_res$d)
  cond_num <- ifelse(min_sv == 0, Inf, max_sv / min_sv)
  tol_cond <- 1e+3
  if(block){
  if(n_dim > 10) {
      mat_reg <- mat + diag(1e-6, nrow(mat))
      V1_inv <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
      temp3 <- 1
      for(i in 1:K) {
        ni <- n[i]  # 第 i 个聚类的实际样本数
        # 提取该聚类对应的 V1 子矩阵
        V_block <- mat[temp3:(temp3 + ni - 1), temp3:(temp3 + ni - 1)]
        # 加正则化，确保该块正定
        V_block_reg <- V_block + diag(1e-5, ni)  # 1e-3 可根据实际情况调整
        V_block_inv <- tryCatch({
          cholV <- chol(V_block_reg)
          chol2inv(cholV)
        }, error = function(e) {
          MASS::ginv(V_block_reg)
        })
        V1_inv[temp3:(temp3 + ni - 1), temp3:(temp3 + ni - 1)] <- V_block_inv
        temp3 <- temp3 + ni
      }
      result_rev=V1_inv
  }
    else {
      if(cond_num>tol_cond){
      result_rev <- tryCatch({
        chol(mat)
      }, error = function(e) {
        ginv(mat)
      })
      }
      else{
        result_rev<-solve(mat)
      }
    }

  }
  else{
    if(n_dim>10){

      if(cond_num>tol_cond){
        mat <- mat + diag(1e-6, nrow(mat))

        tol <- max(dim(mat)) * max(svd_V1$d) * .Machine$double.eps
        d_inv <- ifelse(svd_V1$d > tol, 1/svd_V1$d, 0)
        result_rev <- svd_V1$v %*% diag(d_inv) %*% t(svd_V1$u)

      }
      else{
        result_rev <- tryCatch({
            chol(mat)
        }, error = function(e) {
          ginv(mat)
        })
      }
    }



    else if(cond_num>tol_cond){
      mat <- mat + diag(1e-6, nrow(mat))
      result_rev<-ginv(mat)
    }
    else{
      result_rev<-solve(mat)
    }

  }
  return(result_rev)
  }
