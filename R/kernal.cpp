#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
// [[Rcpp::depends(RcppEigen)]]
//[[Rcpp::export]]
VectorXd kernal(VectorXd dd_c, VectorXd g11_c, VectorXd t11_c, VectorXd tt1_c,
           VectorXd t2_c, VectorXd betainit_h, MatrixXd x111_c, MatrixXd xxx_c,
           RowVectorXd c1_c, VectorXd n_c, VectorXd id_c,
           VectorXd Kn_ls_c, int kk_c, int Kn_c, int K_c){
  float rho;
  int s1=1;
  int s2=1;
  int s3=1;
  VectorXd gSS =Zero(kk_c);
  VectorXd gSS1 =Ones(kk_c);
  VectorXd gSSS1 =Zero(kk_c);
  MatrixXd geebeta (xxx_c.cols(),1);
  int KK1;
  min_ori=2147483646;
  VectorXd gSS2 (Kn_c);
  VectorXd gSS3 (Kn_c);
  VectorXd betainit_c (betainit_h.size());
  betainit_c=betainit_h;
  int kk1=0;
  RowVectorXd Lambda (Kn_c);
  MatrixXd W1 (Kn_c,Kn_c);
  int SK1;
  VectorXd beta1 (xxx_c.cols());
  MatrixXd X1 (xxx_c.rows(),xxx_c.cols());
  VectorXd mu (xxx_c.rows());
  RowVectorXd newY1 (xxx_c.rows());
  MatrixXd res (K_c,int(n_c[0]));
  MatrixXd resm (int(n_c[0]),K_c);
  VectorXd res_t (K_c);
  float rres=0;
  float pphi;
  int p=0;
  VectorXd temp_2 =Ones(n_c.size()-1);
  int SK=1;
  VectorXd mu_temp (int(n_c[0]));
  MatrixXd D1 (xxx_c.rows(),xxx_c.cols());
  VectorXd X1_2 =Ones(int(n_c[0]));
  int temp_3=0;
  VectorXd temp_4 =Ones(int(n_c(j)));
  VectorXd S1 (mu.size());
  MatrixXd R1(int(n_c(0)),int(n_c(0)));
  MatrixXd V1(xxx_c.rows(),xxx_c.rows());
  int temp_6=0;
  int temp_3=0;
  MatrixXd temp_7(int(n_c(j-1)),int(n_c(j-1)));
  VectorXd Z1 (xxx_c.rows());
  MatrixXd b_v (xxx_c.rows(),xxx_c.rows());
  MatrixXd b_2 (xxx_c.cols(),xxx_c.cols());
  MatrixXd solve_V1 (xxx_c.rows(),xxx_c.rows());

  //while (s1==1){
    gSS = Zero(kk_c);
    gSS1 = Ones(kk_c);
    gSSS1 = Zero(kk_c);
    KK1=1;
    min_ori=2147483646;
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
    gSS2.setConstant(0);
    gSS3.setConstant(0);
    for(int i=0;i<Kn_c;i++){
      kk1=0;
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

    Lambda=gSS3;
    W1=Lambda.asDiagonal();
    SK1=1;
    beta1.setConstant(0);
    X1=xxx_c;
    while(s2==1){
      mu=(X1*betainit_c).array().exp();
      newY1=c1_c.array()/Lambda.array();
      res_t=(newY1.transpose().array()-mu.array()).array()/mu.array().sqrt();
      rres=0;
      pphi=(res.array()*res.array()).sum()/(n_c.sum()-X1.cols());
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
      temp_2.setConstant(1);
      rho=(pphi*(-1))*rres/(n_c.array()*(n_c-temp_2).array()).sum()/2-X1.cols();
      SK=1;
      while(s3==1){
        mu_temp.setConstant(0);
        D1.setConstant(0);
        X1_2.setConstant(0);
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
        temp_4.setConstant(1);
        D1.row(j)=mu_temp.asDiagonal().toDenseMatrix()*temp_4.asDiagonal().toDenseMatrix()*X1_2;
        }
        S1=newY1.transpose()-mu;
        R1.setConstant(rho);
        V1.setConstant(0);
        mu_temp.setConstant(0);
        temp_6=0;
        for(int j=1;j<=K_c;j++){
          temp_3=0;
          for(int i=0;i<id_c.size();i++){
            if (id_c(i)==j){
              mu_temp(temp_3)=mu(i);
              temp_3++;
            }
          }
          temp_7.setConstant(rho);
          R1=temp_7;
          R1.diagonal().setConstant(1);
          V1.block(temp_6,temp_6,xxx_c.rows()/K_c,xxx_c.rows()/K_c)=mu_temp.asDiagonal()*R1*mu_temp.asDiagonal()*pphi;
          temp_6+=xxx_c.rows()/K_c;
        }

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
        else {s3=0;}
      }
      if (abs((geebeta-beta1).maxCoeff())>0.000001 && SK1<=30){
        beta1=betainit_c;
        SK1++;
      }
      else {s2=0;}
    }
    if(abs((geebeta-beta2).maxCoeff())>0.000001 || abs((gSS1-gSSS1).maxCoeff())>0.000001){
      beta2=betainit_c;
      gSSS1=gSS1;
      KK1++;
    }
    //else {s1=0;}
    //}
  VectorXd output (betainit_c.size()+1);
  output.segment(0,betainit_c.size()-1)=betainit_c;
  output(betainit_c.size())=rho;
  return output;
  }

