#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#include<Rcpp.h>

#include <cmath>

using namespace Rcpp;
using namespace Eigen;
static Eigen::MatrixXd blockInverseChol(const Eigen::MatrixXd &block, double reg = 1e-5) {
  Eigen::MatrixXd blockReg = block + reg * MatrixXd::Identity(block.rows(), block.cols());
  Eigen::LLT<Eigen::MatrixXd> llt(blockReg);
  if(llt.info() == Eigen::NumericalIssue) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(blockReg, ComputeThinU | ComputeThinV);
    VectorXd d = svd.singularValues();
    MatrixXd U = svd.matrixU();
    MatrixXd V = svd.matrixV();
    double tol = 1e-8;
    VectorXd inv_d(d.size());
    for(int i=0; i<d.size(); i++){
      inv_d[i] = (d[i] > tol) ? 1.0/d[i] : 0.0;
    }
    return (V * inv_d.asDiagonal() * U.transpose());
  } else {
    MatrixXd L = llt.matrixL();
    MatrixXd invBlock = llt.solve(MatrixXd::Identity(blockReg.rows(), blockReg.cols()));
    return invBlock;
  }
}

static Eigen::MatrixXd safeInverse(const Eigen::MatrixXd &M, double reg = 1e-6, double tol_cond = 1e3) {
  Eigen::JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
  double min_sv = svd.singularValues().tail(1)(0);
  double max_sv = svd.singularValues()(0);
  double cond_num = (min_sv == 0.0 ? INFINITY : max_sv / min_sv);
  if(cond_num > tol_cond) {
    MatrixXd M_reg = M + reg * MatrixXd::Identity(M.rows(), M.cols());
    Eigen::JacobiSVD<MatrixXd> svd2(M_reg, ComputeThinU | ComputeThinV);
    VectorXd d = svd2.singularValues();
    MatrixXd U = svd2.matrixU();
    MatrixXd V = svd2.matrixV();
    VectorXd inv_d(d.size());
    double eps = 1e-8;
    for(int i=0; i<d.size(); i++){
      inv_d[i] = (d[i] > eps) ? 1.0/d[i] : 0.0;
    }
    return (V * inv_d.asDiagonal() * U.transpose());
  } else {
    return M.inverse();
  }
}
// [[Rcpp::export]]
SEXP matSolCpp(const Eigen::Map<Eigen::MatrixXd> & mat,
               const bool block,
               const int K,
               const IntegerVector &n,
               const double tol_cond = 1e3)
{
  int n_dim = mat.cols();
  if(block) {
    if(n_dim > 10) {
      MatrixXd mat_reg = mat + 1e-6 * MatrixXd::Identity(mat.rows(), mat.cols());
      MatrixXd V1_inv = MatrixXd::Zero(mat.rows(), mat.cols());
      int temp3 = 0;
      for(int i=0; i < K; i++){
        int ni = n[i];
        MatrixXd V_block = mat_reg.block(temp3, temp3, ni, ni);
        MatrixXd V_block_inv = blockInverseChol(V_block, 1e-5);
        V1_inv.block(temp3, temp3, ni, ni) = V_block_inv;
        temp3 += ni;
      }
      return Rcpp::wrap(V1_inv);
    } else {

      return Rcpp::wrap(safeInverse(mat, 1e-6, tol_cond));
    }
  } else {

    if(n_dim > 10){

      Eigen::JacobiSVD<MatrixXd> svd(mat, ComputeThinU | ComputeThinV);
      VectorXd d = svd.singularValues();
      double min_sv = d[d.size()-1];
      double max_sv = d[0];
      double cond_num = (min_sv == 0 ? INFINITY : max_sv/min_sv);
      if(cond_num>tol_cond){
        MatrixXd mat_reg = mat + 1e-6 * MatrixXd::Identity(mat.rows(), mat.cols());

        Eigen::JacobiSVD<MatrixXd> svd2(mat_reg, ComputeThinU | ComputeThinV);
        VectorXd d2 = svd2.singularValues();
        MatrixXd U2 = svd2.matrixU();
        MatrixXd V2 = svd2.matrixV();
        double eps = 1e-8;
        for(int i=0; i<d2.size(); i++){
          if(d2[i] <= eps) d2[i] = 0.0;
          else d2[i] = 1.0/d2[i];
        }
        MatrixXd result = V2 * d2.asDiagonal() * U2.transpose();
        return Rcpp::wrap(result);
      } else {

        try{
          Eigen::LLT<MatrixXd> llt(mat);
          if(llt.info() == Eigen::NumericalIssue){

            return Rcpp::wrap(safeInverse(mat, 1e-6, tol_cond));
          } else {
            MatrixXd invMat = llt.solve(MatrixXd::Identity(mat.rows(), mat.cols()));
            return Rcpp::wrap(invMat);
          }
        } catch(...) {
          return Rcpp::wrap(safeInverse(mat, 1e-6, tol_cond));
        }
      }
    } else {

      return Rcpp::wrap(safeInverse(mat, 1e-6, tol_cond));
    }
  }
}
