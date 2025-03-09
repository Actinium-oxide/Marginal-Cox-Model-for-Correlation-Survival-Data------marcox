#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

SEXP matSolCpp(const Eigen::Map<Eigen::MatrixXd> & mat,
               const bool block,
               const int K,
               const IntegerVector &n,
               const double tol_cond = 1e3);

// [[Rcpp::export]]
List marcox_iterCpp(
    const Eigen::Map<Eigen::MatrixXd> &X1,
    Eigen::Map<Eigen::VectorXd>       betainit,
    Eigen::Map<Eigen::VectorXd>       Lambda,
    const NumericVector               &c1,
    const NumericMatrix               &W1,
    const NumericVector               &id,
    const IntegerVector               &new_uid,
    const IntegerVector               &n,
    const double tol   = 1e-6,
    const int maxIter  = 30,
    const int maxInner = 500,
    double pphi        = 1.0,
    double rho         = 0.0
){

  const int Kn = X1.rows();
  const int p  = X1.cols();
  const int K  = new_uid.size();



  Map<const MatrixXd> W1_eig(REAL(W1), W1.nrow(), W1.ncol());



  VectorXd mu(Kn);
  for(int i=0; i<Kn; i++){
    double val = 0.0;
    for(int j=0; j<p; j++){
      val += X1(i,j) * betainit[j];
    }
    mu[i] = std::exp(val);
  }


  int SK1 = 1;
  VectorXd beta1 = betainit;


  for(; SK1 <= maxIter; SK1++){

    bool outerBreak = false;
    VectorXd beta1_local = betainit;
    while(true){
      VectorXd newY1(Kn);
      for(int i=0; i<Kn; i++){
        newY1[i] = c1[i]/Lambda[i];
      }

      VectorXd res(Kn);
      for(int i=0; i<Kn; i++){
        double sq = (mu[i] > 0) ? std::sqrt(mu[i]) : 1e-8;
        res[i] = (newY1[i] - mu[i]) / sq;
      }

      double sumRes2 = 0.0;
      for(int i=0; i<Kn; i++){
        sumRes2 += res[i]*res[i];
      }
      pphi = sumRes2 / (Kn - p);

      double rres = 0.0;

      for(int i=0; i<K; i++){

        std::vector<double> groupRes;

        int cluster_id = new_uid[i];

        for(int row=0; row<Kn; row++){
          if( (int)id[row] == cluster_id ){
            groupRes.push_back(res[row]);
          }
        }
        int ni = (int)groupRes.size();
        if(ni==1){
          rres += groupRes[0];
        } else {
          for(int j=0; j<ni-1; j++){
            double tmp = groupRes[j];
            for(int k=j+1; k<ni; k++){
              tmp += groupRes[k];
            }
            rres += groupRes[j]* ( /*sum from j+1..ni-1*/ (tmp - groupRes[j]) );
          }
        }
      }


      double sum_n_n1 = 0.0;
      for(int i=0; i<K; i++){
        double nd = (double)n[i];
        sum_n_n1 += nd*(nd-1.0);
      }
      double denom = sum_n_n1/2.0 - p;
      rho = (1.0/pphi) * rres / denom;


      int SK=1;
      while(true){

        MatrixXd D1 = MatrixXd::Zero(Kn, p);
        {
          int temp2=0;
          for(int i=0; i<K; i++){
            int cluster_id = new_uid[i];
            int ni = n[i];
            VectorXd mu_gp(ni);

            int countRow=0;
            for(int row=0; row<Kn; row++){
              if( (int)id[row] == cluster_id ){
                mu_gp[countRow] = mu[row];
                countRow++;
              }
            }
            MatrixXd diagMu = MatrixXd::Zero(ni, ni);
            for(int r=0; r<ni; r++){
              diagMu(r,r) = mu_gp[r] + 1e-6;
            }
            MatrixXd Xsub(ni, p);
            countRow=0;
            for(int row=0; row<Kn; row++){
              if( (int)id[row] == cluster_id ){
                for(int col=0; col<p; col++){
                  Xsub(countRow, col) = X1(row, col);
                }
                countRow++;
              }
            }
            MatrixXd block = diagMu * Xsub;
            D1.block(temp2, 0, ni, p) = block;
            temp2 += ni;
          }
        }


        VectorXd S1(Kn);
        for(int i=0; i<Kn; i++){
          S1[i] = newY1[i] - mu[i];
        }


        MatrixXd V1 = MatrixXd::Zero(Kn, Kn);
        {
          int temp3=0;
          for(int i=0; i<K; i++){
            int cluster_id = new_uid[i];
            int ni = n[i];
            VectorXd mu_gp(ni);
            int countRow=0;
            for(int row=0; row<Kn; row++){
              if( (int)id[row] == cluster_id ){
                mu_gp[countRow] = mu[row];
                countRow++;
              }
            }

            MatrixXd R1 = MatrixXd::Constant(ni, ni, rho);
            for(int d=0; d<ni; d++){
              R1(d,d) = 1.0;
            }

            MatrixXd sqrtDiag = MatrixXd::Zero(ni, ni);
            for(int r=0; r<ni; r++){
              double val = (mu_gp[r] < 1e-8)? 1e-8 : mu_gp[r];
              sqrtDiag(r,r) = std::sqrt(val);
            }
            MatrixXd block = sqrtDiag * R1 * sqrtDiag * pphi;
            V1.block(temp3, temp3, ni, ni) = block;
            temp3 += ni;
          }
        }


        VectorXd Z1 = D1*betainit + S1;

        SEXP sv = matSolCpp( Map<MatrixXd>(V1.data(), V1.rows(), V1.cols()),
                             true,
                             K,
                             n,
                             1e3 );

        Map<MatrixXd> sol_V1(REAL(sv), V1.rows(), V1.cols());

        MatrixXd tD1 = D1.transpose();
        MatrixXd tmp = sol_V1 * W1_eig;
        MatrixXd M   = tD1 * tmp * D1;

        SEXP sM = matSolCpp( Map<MatrixXd>(M.data(), M.rows(), M.cols()),
                             false,
                             K,
                             n,
                             1e3 );
        Map<MatrixXd> M_inv(REAL(sM), M.rows(), M.cols());
        VectorXd right = tD1 * (sol_V1 * (W1_eig * Z1));
        VectorXd geebeta(p);
        geebeta = M_inv * right;


        VectorXd diff = geebeta - betainit;
        double maxDiff = diff.cwiseAbs().maxCoeff();
        if(maxDiff > 1e-6 && SK <= maxInner) {

          betainit = geebeta;

          for(int i=0; i<Kn; i++){
            double val=0.0;
            for(int col=0; col<p; col++){
              val += X1(i,col)*betainit[col];
            }
            mu[i] = std::exp(val);
          }
          SK++;
        } else {
          break;
        }
      }


      VectorXd diff2 = betainit - beta1_local;
      double md2 = diff2.cwiseAbs().maxCoeff();
      if(md2>1e-6 && SK1 < maxIter){
        beta1_local = betainit;
      } else {

        outerBreak = true;
        break;
      }
    }

    if(outerBreak){
      break;
    }
  }

  return List::create(
    _["betainit"] = betainit,
    _["mu"]       = mu,
    _["pphi"]     = pphi,
    _["rho"]      = rho,
    _["SK1"]      = SK1
  );
}
