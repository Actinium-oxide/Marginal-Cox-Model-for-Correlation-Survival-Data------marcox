// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// marcox_iter_Cpp
List marcox_iter_Cpp(Eigen::Map<Eigen::VectorXd>& betainit_origin, const Eigen::Map<Eigen::MatrixXd>& X1, Eigen::Map<Eigen::VectorXd> betainit, Eigen::Map<Eigen::VectorXd> Lambda, const NumericVector& c1, const NumericMatrix& W1, const NumericVector& id, const IntegerVector& new_uid, const IntegerVector& n, const float& tol, const short& maxIter, const short& maxInner, double pphi, double rho, Eigen::Map<Eigen::VectorXd> rho_vec_k, const int& kv, Eigen::Map<Eigen::MatrixXd> rhomat, const short& method);
RcppExport SEXP _marcox_marcox_iter_Cpp(SEXP betainit_originSEXP, SEXP X1SEXP, SEXP betainitSEXP, SEXP LambdaSEXP, SEXP c1SEXP, SEXP W1SEXP, SEXP idSEXP, SEXP new_uidSEXP, SEXP nSEXP, SEXP tolSEXP, SEXP maxIterSEXP, SEXP maxInnerSEXP, SEXP pphiSEXP, SEXP rhoSEXP, SEXP rho_vec_kSEXP, SEXP kvSEXP, SEXP rhomatSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd>& >::type betainit_origin(betainit_originSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type betainit(betainitSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type W1(W1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type id(idSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type new_uid(new_uidSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const float& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const short& >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< const short& >::type maxInner(maxInnerSEXP);
    Rcpp::traits::input_parameter< double >::type pphi(pphiSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type rho_vec_k(rho_vec_kSEXP);
    Rcpp::traits::input_parameter< const int& >::type kv(kvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type rhomat(rhomatSEXP);
    Rcpp::traits::input_parameter< const short& >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(marcox_iter_Cpp(betainit_origin, X1, betainit, Lambda, c1, W1, id, new_uid, n, tol, maxIter, maxInner, pphi, rho, rho_vec_k, kv, rhomat, method));
    return rcpp_result_gen;
END_RCPP
}
// matSolCpp
SEXP matSolCpp(const Eigen::Map<Eigen::MatrixXd>& mat, const bool block, const int K, const IntegerVector& n, const double tol_cond);
RcppExport SEXP _marcox_matSolCpp(SEXP matSEXP, SEXP blockSEXP, SEXP KSEXP, SEXP nSEXP, SEXP tol_condSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const bool >::type block(blockSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type tol_cond(tol_condSEXP);
    rcpp_result_gen = Rcpp::wrap(matSolCpp(mat, block, K, n, tol_cond));
    return rcpp_result_gen;
END_RCPP
}
// se_cpp
SEXP se_cpp(const Eigen::Map<Eigen::VectorXd>& clusteridx, double betascale, const Eigen::Map<Eigen::VectorXd>& betainit, const Eigen::Map<Eigen::VectorXd>& gSS, int kk, int covnum, int K, const IntegerVector& n, const IntegerVector& new_uid, const Eigen::Map<Eigen::MatrixXd>& xxx, const NumericVector& c1, const Eigen::Map<Eigen::VectorXd>& t2, const Eigen::Map<Eigen::VectorXd>& tt1, const Eigen::Map<Eigen::VectorXd>& gg1, const Eigen::Map<Eigen::VectorXd>& Lambda, const IntegerVector& id, const Eigen::Map<Eigen::MatrixXd>& rhomat, const Eigen::Map<Eigen::MatrixXd>& Q1R, const int kv);
RcppExport SEXP _marcox_se_cpp(SEXP clusteridxSEXP, SEXP betascaleSEXP, SEXP betainitSEXP, SEXP gSSSEXP, SEXP kkSEXP, SEXP covnumSEXP, SEXP KSEXP, SEXP nSEXP, SEXP new_uidSEXP, SEXP xxxSEXP, SEXP c1SEXP, SEXP t2SEXP, SEXP tt1SEXP, SEXP gg1SEXP, SEXP LambdaSEXP, SEXP idSEXP, SEXP rhomatSEXP, SEXP Q1RSEXP, SEXP kvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type clusteridx(clusteridxSEXP);
    Rcpp::traits::input_parameter< double >::type betascale(betascaleSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type betainit(betainitSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type gSS(gSSSEXP);
    Rcpp::traits::input_parameter< int >::type kk(kkSEXP);
    Rcpp::traits::input_parameter< int >::type covnum(covnumSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type new_uid(new_uidSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type xxx(xxxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type c1(c1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type t2(t2SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type tt1(tt1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type gg1(gg1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type id(idSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type rhomat(rhomatSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type Q1R(Q1RSEXP);
    Rcpp::traits::input_parameter< const int >::type kv(kvSEXP);
    rcpp_result_gen = Rcpp::wrap(se_cpp(clusteridx, betascale, betainit, gSS, kk, covnum, K, n, new_uid, xxx, c1, t2, tt1, gg1, Lambda, id, rhomat, Q1R, kv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_marcox_marcox_iter_Cpp", (DL_FUNC) &_marcox_marcox_iter_Cpp, 18},
    {"_marcox_matSolCpp", (DL_FUNC) &_marcox_matSolCpp, 5},
    {"_marcox_se_cpp", (DL_FUNC) &_marcox_se_cpp, 19},
    {NULL, NULL, 0}
};

RcppExport void R_init_marcox(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
