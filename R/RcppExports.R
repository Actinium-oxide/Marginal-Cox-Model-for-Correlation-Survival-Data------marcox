# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

marcox_iterCpp <- function(X1, betainit, Lambda, c1, W1, id, new_uid, n, tol = 1e-6, maxIter = 30L, maxInner = 500L, pphi = 1.0, rho = 0.0) {
    .Call(`_marcox_marcox_iterCpp`, X1, betainit, Lambda, c1, W1, id, new_uid, n, tol, maxIter, maxInner, pphi, rho)
}

matSolCpp <- function(mat, block, K, n, tol_cond = 1e3) {
    .Call(`_marcox_matSolCpp`, mat, block, K, n, tol_cond)
}

sandwich_rcpp <- function(rho, betascale, betainit, gSS, kk, covnum, K, n, new_uid, xxx, c1, t2, tt1, gg1, Lambda, id) {
    .Call(`_marcox_sandwich_rcpp`, rho, betascale, betainit, gSS, kk, covnum, K, n, new_uid, xxx, c1, t2, tt1, gg1, Lambda, id)
}

