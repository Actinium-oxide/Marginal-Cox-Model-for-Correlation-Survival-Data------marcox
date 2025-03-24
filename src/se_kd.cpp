#ifdef __GNUC__
#pragma GCC diagnostic push
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
SEXP sandwich_kd_rcpp(const NumericVector &rho,
                       double betascale,
                       const Eigen::Map<Eigen::VectorXd> &betainit,
                       const Eigen::Map<Eigen::VectorXd> &gSS,
                       int kk,
                       int covnum,
                       int K,
                       const IntegerVector &n,
                       const IntegerVector &new_uid,
                       const Eigen::Map<Eigen::MatrixXd> &xxx,
                       const NumericVector &c1,
                       const Eigen::Map<Eigen::VectorXd> &t2,
                       const Eigen::Map<Eigen::VectorXd> &tt1,
                       const Eigen::Map<Eigen::VectorXd> &gg1,
                       const Eigen::Map<Eigen::VectorXd> &Lambda,
                       const IntegerVector &id,
                       const int &kv)
{
  VectorXd gS(kk);
  gS[0] = gSS[0];
  for(int s=1; s<kk; s++){
    gS[s] = gSS[s] - gSS[s-1];
  }
  int Kn = xxx.rows();
  VectorXd mu2(Kn);
  for(int i=0; i<Kn; i++){
    double tmp = 0.0;
    for(int col=0; col<covnum; col++){
      tmp += xxx(i, col)*betainit[col];
    }
    mu2[i] = std::exp(tmp);
  }
  VectorXd ABC1 = VectorXd::Zero(K);
  MatrixXd VA1 = MatrixXd::Zero(covnum,covnum);
  for(int v=0; v<covnum; v++){
    for(int w=0; w<covnum; w++){
      for(int i_clust=0; i_clust<K; i_clust++){
        int ni = n[i_clust];

        VectorXd indices = VectorXd::LinSpaced(ni, 0, ni - 1);
        MatrixXd Diff = indices.replicate(1, ni) - indices.transpose().replicate(ni, 1);
        MatrixXd absDiff = Diff.cwiseAbs();
        MatrixXd Q1 = absDiff.unaryExpr([rho, kv](double d) -> double {
          int lag = static_cast<int>(d + 0.5);
          return (lag == 0) ? 1.0 : ((lag > kv) ? 0.0 : rho(lag - 1));
        });
        Q1.diagonal().setOnes();

        SEXP ret = matSolCpp(Map<MatrixXd>(Q1.data(), Q1.rows(), Q1.cols()),false,K,n,1e3);
        Map<MatrixXd> IQ1(REAL(ret), Q1.rows(), Q1.cols());
        std::vector<int> idx;
        for(int r=0; r<Kn; r++){
          if(id[r]==new_uid[i_clust]) idx.push_back(r);
        }
        MatrixXd z22(ni, covnum);
        for(int r=0; r<ni; r++){
          z22.row(r) = xxx.row(idx[r]);
        }
        RowVectorXd A2(ni);
        for(int r=0; r<ni; r++){
          A2(r) = z22(r,v);
        }
        VectorXd c22(ni), Lam22(ni), mu22vec(ni);
        for(int r=0; r<ni; r++){
          c22[r] = c1[idx[r]];
          Lam22[r] = Lambda[idx[r]];
          mu22vec[r] = mu2[idx[r]];
        }
        VectorXd sqrt_mu22(ni), inv_sqrt_mu22(ni);
        for(int rr=0; rr<ni; rr++){
          sqrt_mu22[rr] = std::sqrt(mu22vec[rr]);
          inv_sqrt_mu22[rr] = 1.0/std::sqrt(mu22vec[rr]);
        }
        MatrixXd outer_mat(ni, ni);
        for(int rr=0; rr<ni; rr++){
          for(int cc=0; cc<ni; cc++){
            outer_mat(rr, cc) = sqrt_mu22[rr]*inv_sqrt_mu22[cc];
          }
        }
        MatrixXd BB1 = outer_mat.array()*IQ1.array();
        MatrixXd B2 = MatrixXd::Zero(ni, ni);
        for(int s_i=0; s_i<ni; s_i++){
          for(int l=0; l<ni; l++){
            B2(s_i,l) = 0.5*(z22(s_i,w)-z22(l,w))*BB1(s_i,l);
          }
        }
        VectorXd C2_cluster(ni);
        for(int r=0; r<ni; r++){
          C2_cluster[r] = c22[r]/Lam22[r] - mu22vec[r];
        }
        MatrixXd D2 = BB1;
        VectorXd E2(ni);
        for(int r=0; r<ni; r++){
          E2[r] = z22(r,w)*mu22vec[r];
        }
        MatrixXd G2 = Lam22.asDiagonal();
        VectorXd term1 = B2*(G2*C2_cluster);
        VectorXd term2 = D2*(G2*E2);
        double res_val = (A2*(term1-term2))(0);
        ABC1[i_clust] = res_val;
      }
      double sumVal=0.0;
      for(int i=0; i<K; i++){
        sumVal += ABC1[i];
      }
      VA1(v,w) = sumVal/betainit[0];
      VA1(v,w) = sumVal/betascale;
      ABC1.setZero();
    }
  }
  MatrixXd sdm = VA1;
  MatrixXd BBC = MatrixXd::Zero(kk, covnum);
  for(int j=0; j<covnum; j++){
    for(int s=0; s<kk; s++){
      double term1=0.0;
      double term2=0.0;
      for(int r=0; r<Kn; r++){
        if(c1[r]==1.0 && std::fabs(t2[r]-tt1[s])<1e-15){
          double dotval=0.0;
          for(int cc=0; cc<covnum; cc++){
            dotval += xxx(r,cc)*betainit[cc];
          }
          double exp_dot = std::exp(dotval);
          double BCm = gS[s]*exp_dot;
          double fac = (std::exp(-BCm)+BCm*std::exp(-BCm)-1.0)/std::pow((1-std::exp(-BCm)),2);
          term1 += exp_dot*fac*xxx(r,j);
        }
      }
      for(int r=0; r<Kn; r++){
        if(t2[r]>=tt1[s]){
          double dotval=0.0;
          for(int cc=0; cc<covnum; cc++){
            dotval += xxx(r,cc)*betainit[cc];
          }
          double exp_dot = std::exp(dotval);
          term2 += gg1[r]*exp_dot*xxx(r,j);
        }
      }
      BBC(s,j) = term1+term2;
    }
  }
  VectorXd CCC(kk);
  CCC.setZero();
  for(int s=0; s<kk; s++){
    double sum_val=0.0;
    for(int r=0; r<Kn; r++){
      if(c1[r]==1.0 && std::fabs(t2[r]-tt1[s])<1e-15){
        double dotval=0.0;
        for(int cc=0; cc<covnum; cc++){
          dotval += xxx(r,cc)*betainit[cc];
        }
        double exp_dot = std::exp(dotval);
        double BCm = gS[s]*exp_dot;
        sum_val += std::exp(2*dotval-BCm)/std::pow((1-std::exp(-BCm)),2);
      }
    }
    CCC[s] = sum_val;
  }
  MatrixXd BC = MatrixXd::Zero(covnum, kk);
  for(int r=0; r<covnum; r++){
    for(int s=0; s<kk; s++){
      double elem=0.0;
      for(int i_clust=0; i_clust<K; i_clust++){
        int ni = n[i_clust];
        std::vector<int> idx;
        for(int row=0; row<Kn; row++){
          if(id[row]==new_uid[i_clust]) idx.push_back(row);
        }
        VectorXd mu22(ni), colr(ni), t21(ni);
        for(int sub=0; sub<ni; sub++){
          mu22[sub] = mu2[idx[sub]];
          colr[sub] = xxx(idx[sub], r);
          t21[sub] = t2[idx[sub]];
        }

        VectorXd indices = VectorXd::LinSpaced(ni, 0, ni - 1);
        MatrixXd Diff = indices.replicate(1, ni) - indices.transpose().replicate(ni, 1);
        MatrixXd absDiff = Diff.cwiseAbs();
        MatrixXd Q1 = absDiff.unaryExpr([rho, kv](double d) -> double {
          int lag = static_cast<int>(d + 0.5);
          return (lag == 0) ? 1.0 : ((lag > kv) ? 0.0 : rho(lag - 1));
        });
        Q1.diagonal().setOnes();

        SEXP ret2 = matSolCpp(Map<MatrixXd>(Q1.data(),Q1.rows(),Q1.cols()),false,K,n,1e3);
        Map<MatrixXd> IQ1(REAL(ret2), ni, ni);
        for(int j_idx=0; j_idx<ni; j_idx++){
          if(t21[j_idx]>=tt1[s]){
            double sum_val=0.0;
            double inv_sqrt_j = (mu22[j_idx]>1e-15)?1.0/std::sqrt(mu22[j_idx]):0.0;
            for(int k=0; k<ni; k++){
              double sqrt_k = (mu22[k]>1e-15)?std::sqrt(mu22[k]):0.0;
              sum_val += colr[k]*sqrt_k*inv_sqrt_j*IQ1(k,j_idx);
            }
            elem += sum_val*mu22[j_idx]/betascale;
          }
        }
      }
      BC(r,s) = elem;
    }
  }
  MatrixXd M22 = -sdm;
  MatrixXd M23 = BC;
  MatrixXd M32 = BBC;
  MatrixXd M33 = CCC.asDiagonal();
  MatrixXd M(M22.rows()+M32.rows(), M22.cols()+M23.cols());
  M.block(0,0,M22.rows(),M22.cols())=M22;
  M.block(0,M22.cols(),M23.rows(),M23.cols())=M23;
  M.block(M22.rows(),0,M32.rows(),M32.cols())=M32;
  M.block(M22.rows(),M22.cols(),M33.rows(),M33.cols())=M33;
  MatrixXd fdm = MatrixXd::Zero(covnum+kk,covnum+kk);
  for(int i_clust=0; i_clust<K; i_clust++){
    std::vector<int> idx;
    for(int row=0; row<Kn; row++){
      if(id[row]==new_uid[i_clust]) idx.push_back(row);
    }
    int ni=idx.size();
    MatrixXd xxx1_mat(ni,covnum);
    VectorXd gg11(ni), c111(ni), t21(ni), g11(ni), mu22cl(ni), c22(ni);
    MatrixXd z22(ni, covnum);
    for(int k=0; k<ni; k++){
      int glob=idx[k];
      for(int cc=0; cc<covnum; cc++){
        xxx1_mat(k, cc) = xxx(glob, cc);
        z22(k, cc) = xxx(glob, cc);
      }
      gg11[k] = gg1[glob];
      c111[k] = c1[glob];
      t21[k] = t2[glob];
      g11[k] = Lambda[glob];
      mu22cl[k] = mu2[glob];
      c22[k] = c1[glob];
    }
    MatrixXd mu22m = mu22cl.asDiagonal();
    MatrixXd G2 = g11.asDiagonal();
    VectorXd C2_cluster(ni);
    for(int r=0; r<ni; r++){
      C2_cluster[r] = c22[r]/g11[r]-mu22cl[r];
    }

    VectorXd indices = VectorXd::LinSpaced(ni, 0, ni - 1);
    MatrixXd Diff = indices.replicate(1, ni) - indices.transpose().replicate(ni, 1);
    MatrixXd absDiff = Diff.cwiseAbs();
    MatrixXd Q1 = absDiff.unaryExpr([rho, kv](double d) -> double {
      int lag = static_cast<int>(d + 0.5);
      return (lag == 0) ? 1.0 : ((lag > kv) ? 0.0 : rho(lag - 1));
    });
    Q1.diagonal().setOnes();

    VectorXd sqrt_mu22cl(ni);
    for(int r=0; r<ni; r++){
      sqrt_mu22cl[r] = std::sqrt(mu22cl[r]);
    }
    MatrixXd sqrt_mu22m = sqrt_mu22cl.asDiagonal();
    MatrixXd tempMat = sqrt_mu22m*Q1*sqrt_mu22m;
    tempMat *= betascale;
    SEXP ret3 = matSolCpp(Map<MatrixXd>(tempMat.data(), ni, ni),false,K,n,1e3);
    Map<MatrixXd> sol_temp(REAL(ret3), ni, ni);
    MatrixXd Aleft = (mu22m*z22).transpose()*sol_temp;
    VectorXd AC2 = G2*C2_cluster;
    VectorXd fdv = Aleft*AC2;
    VectorXd eqalpha = VectorXd::Zero(kk);
    for(int s=0; s<kk; s++){
      std::vector<int> Aalpha, Balpha;
      for(int r=0; r<ni; r++){
        if(std::fabs(t21[r]-tt1[s])<1e-15 && c111[r]==1.0){
          Aalpha.push_back(r);
        }
        if(t21[r]>=tt1[s]){
          Balpha.push_back(r);
        }
      }
      if(Balpha.empty()){
        eqalpha[s]=0.0;
      } else if(!Balpha.empty() && Aalpha.empty()){
        double tmp=0.0;
        for(size_t b=0; b<Balpha.size(); b++){
          tmp += gg11[Balpha[b]]*mu22cl[Balpha[b]];
        }
        eqalpha[s] = -tmp;
      } else {
        double sumA=0.0;
        for(size_t a=0; a<Aalpha.size(); a++){
          int idxA=Aalpha[a];
          sumA += mu22cl[idxA]/(1-std::exp(-gS[s]*mu22cl[idxA]));
        }
        double sumB=0.0;
        for(size_t b=0; b<Balpha.size(); b++){
          sumB += gg11[Balpha[b]]*mu22cl[Balpha[b]];
        }
        eqalpha[s] = sumA - sumB;
      }
    }
    VectorXd cat_vec(fdv.size()+eqalpha.size());
    for(int i=0; i<fdv.size(); i++){
      cat_vec[i] = fdv[i];
    }
    for(int i=0; i<eqalpha.size(); i++){
      cat_vec[fdv.size()+i] = eqalpha[i];
    }
    fdm += cat_vec*cat_vec.transpose();
  }
  int rowM = M.rows();
  int colM = M.cols();
  SEXP retM = matSolCpp(Map<MatrixXd>(M.data(), rowM, colM),false,K,n,1e3);
  Map<MatrixXd> sol_M(REAL(retM), rowM, colM);
  MatrixXd vcmR = sol_M*fdm*sol_M.transpose();
  VectorXd V1(covnum);
  for(int i=0; i<covnum; i++){
    V1[i] = vcmR(i,i);
  }
  return wrap(V1);
}
