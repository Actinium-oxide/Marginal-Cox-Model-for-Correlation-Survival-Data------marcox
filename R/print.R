#' Print method for marcox
#'
#' @param x A marcox object.
#' @param ... Further arguments passed to or from other methods.
#' @exportS3Method  print marcox
#' @noRd
print.marcox <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat('\n')
  print(x$Estimation)
  cat("\n")
  cat('Likelihood ratio test= ',x$.loglik,' on ',x$.df,' df, p= ',x$.p_value,'\n')
  cat("n =", x$.cens1 + x$.cens0, ", number of events =", x$.cens1, "\n")

  invisible(x)
}
