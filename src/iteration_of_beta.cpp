#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#include <RcppEigen.h>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

using namespace Rcpp;
using namespace Eigen;

SEXP matSolCpp(const Eigen::Map<Eigen::MatrixXd> & mat,
               const bool block,
               const int K,
               const IntegerVector &n,
               const double tol_cond);

// [[Rcpp::export]]
List marcox_iter_Cpp(
    Eigen::Map<Eigen::VectorXd>       &betainit_origin,
    const Eigen::Map<Eigen::MatrixXd> &X1,
    Eigen::Map<Eigen::VectorXd>       betainit,
    Eigen::Map<Eigen::VectorXd>       Lambda,
    const NumericVector               &c1,
    const NumericMatrix               &W1,
    const NumericVector               &id,
    const IntegerVector               &new_uid,
    const IntegerVector               &n,
    const float &tol ,
    double pphi,
    double rho,
    Eigen::Map<Eigen::VectorXd> rho_vec_k,
    const int &kv,
    Eigen::Map<Eigen::MatrixXd> rhomat,
    const short &method,
    int SK1

){

  const int Kn = X1.rows();
  const int p  = X1.cols();
  const int K  = new_uid.size();
  Eigen::VectorXi clusteridx(K);
  clusteridx[0]=0;
  for(int i=1;i<K;i++){
    clusteridx[i]=clusteridx[i-1]+n[i-1];
  }

  Map<const MatrixXd> W1_eig(REAL(W1), W1.nrow(), W1.ncol());



  VectorXd mu(Kn);
  for(int i=0; i<Kn; i++){
    double val = 0.0;
    for(int j=0; j<p; j++){
      val += X1(i,j) * betainit[j];
    }
    mu[i] = std::exp(val);
  }


  VectorXd beta1_local = VectorXd::Zero(X1.cols());


  while(true){

    bool outerBreak = false;
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

      //indp rho iter
      if(method==5){
        rho=0;
        rhomat.diagonal().setOnes();
      }



      // exchangeable & ar1 rho iter
      else if (method==0 || method==1 ){
        for(int i=0; i<K; i++){
        //int idx=clusteridx[i];
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
            rres += groupRes[j]* ((tmp - groupRes[j]) );
          }
        }
        }
        if(method==0){
        double sum_n_n1 = 0.0;
        for(int i=0; i<K; i++){
          double nd = (double)n[i];
          sum_n_n1 += nd*(nd-1.0);
        }
        double denom = (sum_n_n1/2.0) - p;
        rho = (1.0/pphi) * rres / denom;

        rhomat.setConstant(rho);
        rhomat.diagonal().setOnes();
        }
        //ar1 rho iter
        else if(method==1){
          double rres = 0.0;
          for (int i = 0; i < K; i++) {
            std::vector<double> groupRes;
            int cluster_id = new_uid[i];
            for (int row = 0; row < Kn; row++) {
              if ((int)id[row] == cluster_id) {
                groupRes.push_back(res[row]);
              }
            }
            int ni = (int)groupRes.size();
            if (ni > 1) {
              for (int j = 0; j < ni - 1; j++) {
                rres += groupRes[j] * groupRes[j + 1];
              }
            }
          }
          int sum_n_n1 = 0.0;
          for (int i = 0; i < K; i++) {
            sum_n_n1 += n[i] - 1.0;
          }
          int denom = sum_n_n1 - p;
          rho = (1.0 / pphi) * (rres / denom);


          for(int i=0;i<K;i++){
            int ni = n[i];
            int idx=clusteridx[i];
            VectorXd indices = VectorXd::LinSpaced(ni, 0, ni - 1);
            MatrixXd absDiff = (indices.replicate(1, ni) - indices.transpose().replicate(ni, 1)).cwiseAbs();
            MatrixXd rhomattp = absDiff.unaryExpr([rho](double d) { return std::pow(rho, d); });
            rhomattp.diagonal().setOnes();
            rhomat.block(idx,idx,ni,ni) = rhomattp;
            rhomattp.resize(0,0);
          }
      }
      }

      //toep & kd rho iter
      else if(method==3||method==2){

        for (int m = 1; m <= kv; m++) {
          double rres = 0.0;
          double sum_n_n1 = 0.0;
          for (int i = 0; i < K; i++) {
            int cluster_id = new_uid[i];
            std::vector<double> groupRes;
            groupRes.reserve(n[i]);
            for (int row = 0; row < Kn; row++) {
              if ((int)id[row] == cluster_id) {
                groupRes.push_back(res[row]);
              }
            }
            int ni = groupRes.size();
            if (ni > m) {
              for (int j = 0; j < ni - m; j++) {
                rres += groupRes[j] * groupRes[j + m];
              }
            }
          }
          for (int i = 0; i < K; i++) {
            sum_n_n1 += n[i] - m;
          }
          double denom = sum_n_n1 - p;
          double rho_est;
          if (denom > 0)
            rho_est = (1.0 / pphi) * (rres / denom);
          else
            rho_est = 0.0;

          if (rho_est < -1.0) {
            rho_vec_k[m - 1] = -1.0;
          } else if (rho_est > 1.0) {
            rho_vec_k[m - 1] = 1.0;
          } else {
            double alpha_damp = 1.0;
            rho_vec_k[m - 1] = alpha_damp * rho_est + (1.0 - alpha_damp) * rho_vec_k[m - 1];
          }
        }

        for(int i=0;i<K;i++){
          int ni = n[i];
          int idx=clusteridx[i];
          VectorXd indices = VectorXd::LinSpaced(ni, 0, ni - 1);
          MatrixXd Diff = indices.replicate(1, ni) - indices.transpose().replicate(ni, 1);
          MatrixXd absDiff = Diff.cwiseAbs();
          MatrixXd rhomattp = absDiff.unaryExpr([&rho_vec_k, kv](double d) -> double {
            int lag = static_cast<int>(d + 0.5);
            return (lag == 0) ? 1.0 : ((lag > kv) ? 0.0 : rho_vec_k(lag - 1));
          });
          rhomattp.diagonal().setOnes();
          rhomat.block(idx,idx,ni,ni)=rhomattp;
          rhomattp.resize(0,0);
      }
      }


      //uns rho iter

      else if(method==4) {

        for(int ii = 0; ii < K; ii++){
          int ni  = n[ii];
          int idx = clusteridx[ii];
          std::vector<double> groupRes;
          groupRes.reserve(ni);
          int cluster_id = new_uid[ii];
          for(int row = 0; row < Kn; row++){
            if((int)id[row] == cluster_id){
              groupRes.push_back(res[row]);
            }
          }
          for(int r = 0; r < ni; r++){
            for(int q = r + 1; q < ni; q++){


              double rres    = groupRes[r] * groupRes[q];
              double denom   = K-p;
              double rho_est = 0.0;


                rho_est = (1.0 / pphi) * (rres / denom);
              if(rho_est < -1.0) rho_est = -1.0;
              if(rho_est > 1.0) rho_est = 1.0;

              rhomat(idx + r, idx + q) = rho_est;
              rhomat(idx + q, idx + r) = rho_est;
            }
          }
        }
        rhomat.diagonal().setOnes();
      }

      int SK=1;
      MatrixXd R1(0,0);
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

            // MatrixXd R1 = MatrixXd::Constant(ni, ni, rho);
            // for(int d=0; d<ni; d++){
            //   R1(d,d) = 1.0;
            // }

            R1.resize(ni,ni);
            int idx=clusteridx[i];
            R1=rhomat.block(idx,idx,ni,ni);

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
        if(maxDiff > tol && SK <= 1000) {

          betainit = geebeta;
          //std::cout <<  betainit.transpose() << std::endl;

          double divergence = betainit.cwiseAbs().maxCoeff();
          if(divergence > 5.0 || divergence < -3){
            Rcpp::stop("Ran out of iterations and did not converge");
          }

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
      if(md2>tol && SK1 < 30){
        beta1_local = betainit;
        SK1++;
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
    _["SK1"]      = SK1,
    _["rhomat"]   = rhomat,
    _["clusteridx"]= clusteridx,
    _["rho_vec_k"] = rho_vec_k
  );
}
