#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;
using namespace Eigen;
using namespace std;
// [[Rcpp::depends(RcppEigen)]]
//[[Rcpp::export]]
VectorXd kernal_1(VectorXd dd_c, VectorXd g11_c, VectorXd t11_c, VectorXd tt1_c,
           VectorXd t2_c, VectorXd betainit_c, MatrixXd x111_c, MatrixXd xxx_c,
           VectorXd c1_c, VectorXd n_c, VectorXd id_c,
           VectorXd Kn_ls_c, int kk_c, int Kn_c, int K_c){
  while (true){
    VectorXd gSS (kk_c);
    VectorXd gSS1 (kk_c);
    gSS1.setConstant(1);
    float rho;
    VectorXd gSSS1 (kk_c);
    MatrixXd geebeta (xxx_c.cols(),1);
    int KK1=1;
    MatrixXd beta2 (xxx_c.cols(),1);
    int min_ori=2147483647;
    MatrixXd tm_3 (1,x111_c.cols());
    for(int j=0;j<t11_c.size();j++){
      if(t11_c(j)==tt1_c(0) && Kn_ls_c(j)-1<min_ori) {min_ori=Kn_ls_c(j)-1;}
    }
    gSS(0)=dd_c(0)/exp(((g11_c.segment(min_ori,Kn_c-min_ori))*
      (x111_c.block(min_ori,0,Kn_c-min_ori,x111_c.cols())*betainit_c)).array()).sum();
    for(int p=1;p<kk_c;p++){
      for(int j=0;j<t11_c.size();j++){
        if(t11_c(j)==tt1_c[p] && Kn_ls_c(j)-1<min_ori) {min_ori=Kn_ls_c(j)-1;}
      }
      gSS[p]=gSS[p-1]+dd_c[p]/exp(((g11_c.segment(min_ori,Kn_c-min_ori))*
        ((x111_c.block(min_ori,0,Kn_c-min_ori,x111_c.cols())*betainit_c))).array()).sum();
    }
    gSS1=exp((-1*gSS).array());
    VectorXd gSS2 (Kn_c);
    VectorXd gSS3 (Kn_c);
    for(int i=0;i<Kn_c;i++){
      int kk1=0;
      if(t2_c(i)<tt1_c(i)){
        gSS2(i)=1;
        gSS3(i)=0.00000001;
      }
      else{
        if(t2_c(i)>=tt1_c(kk_c-1)){
          gSS2(i)=0;
          gSS3(i)=gSS(kk_c-1);
        }
      else{
        while(true){
          if(t2_c.coeff(i)>=tt1_c.coeff(kk1)){
            kk1+=1;
          }
          else{break;}
        }
        gSS2(i)=pow((gSS1(kk1-1)),(exp(xxx_c.row(i)*betainit_c)));
        gSS3(i)=gSS(kk1-1);
      }
    }
    }

    VectorXd Lambda=gSS3;
    MatrixXd W1=Lambda.asDiagonal();
    int SK1=1;
    MatrixXd beta1 (xxx_c.cols(),1);
    MatrixXd X1=xxx_c;
    while(true){
      VectorXd mu=(X1*betainit_c).array().exp();
      VectorXd newY1=c1_c.array()/Lambda.array();
      VectorXd res=(newY1.array()-mu.array()).array()/mu.array().sqrt();
      float rres=0;
      float pphi=(res.array()*res.array()).sum()/(n_c.sum()-X1.cols());
      MatrixXd resm (int(n_c.maxCoeff()),K_c);
      int p=0;
      for(int i=0;i<K_c;i++){
        p=0;
        for(int j=0;j<id_c.size();j++){
          if(id_c(j)==i){
            resm(p,i)=res(j);
            p++;}}}
      res=resm.transpose();
      for(int i=0;i<K_c;i++){
        if(n_c(i)==1){
          rres+=res(i,0);
        }
        else{
          for(int j=0;j<n_c(i)-1;j++){
            rres+=res(i,j)*res.row(i).segment(j+1,int(n_c(i))-j).sum();
          }
        }
      }
      VectorXd temp_2 (n_c.size()-1);
      temp_2.setConstant(1);
      rho=(pphi*(-1))*rres/(n_c.array()*(n_c-temp_2).array()).sum()/2-X1.cols();
      int SK=1;
      while(true){
        VectorXd mu_temp;
        MatrixXd D1 (xxx_c.rows(),xxx_c.cols());
        MatrixXd X1_2 (X1.rows(),X1.cols());
        for(int j=1;j<=K_c;j++){
        int temp_3=0;
        mu_temp.resize(int(n_c(j)));
        for(int i=0;i<id_c.size();i++){
          if (id_c(i)==j){
            mu_temp(temp_3)=mu(i);
            X1_2.row(temp_3)=X1.row(i);
            temp_3++;
          }
        }
        VectorXd temp_4 (int(n_c(j)));
        temp_4.setConstant(1);
        D1.row(j)=mu_temp.asDiagonal().toDenseMatrix()*temp_4.asDiagonal().toDenseMatrix()*X1_2;
        }
        VectorXd S1=newY1-mu;
        MatrixXd R1(int(n_c(0)),int(n_c(0)));
        R1.setConstant(rho);
        MatrixXd V1(xxx_c.rows(),xxx_c.rows());
        mu_temp.fill(0);
        int temp_6=0;
        int temp_3=0;

        for(int j=1;j<=K_c;j++){
          temp_3=0;
          for(int i=0;i<id_c.size();i++){
            if (id_c(i)==j){
              mu_temp(temp_3)=mu(i);
              temp_3++;
            }
          }
          MatrixXd temp_7(int(n_c(j-1)),int(n_c(j-1)));
          temp_7.setConstant(rho);
          R1=temp_7;
          R1.diagonal().setConstant(1);
          V1.block(temp_6,temp_6,xxx_c.rows()/K_c,xxx_c.rows()/K_c)=mu_temp.asDiagonal()*R1*mu_temp.asDiagonal()*pphi;
          temp_6+=xxx_c.rows()/K_c;
        }

        MatrixXd Z1 (xxx_c.rows(),1);
        Z1=D1*betainit_c+S1;
        MatrixXd b_v (V1.cols(),V1.cols());
        MatrixXd b_2 (xxx_c.cols(),xxx_c.cols());
        b_v.diagonal().setConstant(1);
        b_2.diagonal().setConstant(1);
        MatrixXd solve_V1 (V1.cols(),V1.cols());
        solve_V1=V1.colPivHouseholderQr().solve(b_v);
        geebeta=(D1.transpose()*solve_V1*W1*D1).colPivHouseholderQr().solve(b_2)*solve_V1*W1*Z1;
        if(abs((geebeta-betainit_c).maxCoeff())>0.000001 && SK<=500){
          betainit_c=geebeta;
          mu=(X1*betainit_c).array().exp();
          SK++;
        }
        else {break;}
      }
      if (abs((geebeta-beta1).maxCoeff())>0.000001 && SK1<=30){
        beta1=betainit_c;
        SK1++;
      }
      else {break;}
    }
    if(abs((geebeta-beta2).maxCoeff())>0.000001 || abs((gSS1-gSSS1).maxCoeff())>0.000001){
      beta2=betainit_c;
      gSSS1=gSS1;
      KK1++;
    }
    else {break;}
    }
  VectorXd output (betainit_c.size()+1);
  output.segment(0,betainit_c.size())=betainit_c;
  output.coeff(output.size()-1)=rho;
  return betainit_c;
  }



#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// kernal_1
VectorXd kernal_1(VectorXd dd_c, VectorXd g11_c, VectorXd t11_c, VectorXd tt1_c, VectorXd t2_c, VectorXd betainit_c, MatrixXd x111_c, MatrixXd xxx_c, VectorXd c1_c, VectorXd n_c, VectorXd id_c, VectorXd Kn_ls_c, int kk_c, int Kn_c, int K_c);
RcppExport SEXP sourceCpp_1_kernal_1(SEXP dd_cSEXP, SEXP g11_cSEXP, SEXP t11_cSEXP, SEXP tt1_cSEXP, SEXP t2_cSEXP, SEXP betainit_cSEXP, SEXP x111_cSEXP, SEXP xxx_cSEXP, SEXP c1_cSEXP, SEXP n_cSEXP, SEXP id_cSEXP, SEXP Kn_ls_cSEXP, SEXP kk_cSEXP, SEXP Kn_cSEXP, SEXP K_cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< VectorXd >::type dd_c(dd_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type g11_c(g11_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type t11_c(t11_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type tt1_c(tt1_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type t2_c(t2_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type betainit_c(betainit_cSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type x111_c(x111_cSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type xxx_c(xxx_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type c1_c(c1_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type n_c(n_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type id_c(id_cSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Kn_ls_c(Kn_ls_cSEXP);
    Rcpp::traits::input_parameter< int >::type kk_c(kk_cSEXP);
    Rcpp::traits::input_parameter< int >::type Kn_c(Kn_cSEXP);
    Rcpp::traits::input_parameter< int >::type K_c(K_cSEXP);
    rcpp_result_gen = Rcpp::wrap(kernal_1(dd_c, g11_c, t11_c, tt1_c, t2_c, betainit_c, x111_c, xxx_c, c1_c, n_c, id_c, Kn_ls_c, kk_c, Kn_c, K_c));
    return rcpp_result_gen;
END_RCPP
}
