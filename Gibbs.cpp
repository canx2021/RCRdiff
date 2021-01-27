//--------------------------------------------------------------
// Header (header)
//--------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <cmath>
using namespace Rcpp;
using namespace arma;


//--------------------------------------------------------------
// Functions (Functions_cpp)
//--------------------------------------------------------------
arma::vec colSums(arma::mat mat) {
  int ncols = mat.n_cols;
  vec res = randn(ncols);
  for (int i=0; i<ncols; i++) {
    res(i) = sum(mat.col(i));
  }
  return(res);
}

arma::vec rowSums(arma::mat mat) {
  int nrows = mat.n_rows;
  vec res = randn(nrows);
  for (int i=0; i<nrows; i++) {
    res(i) = sum(mat.row(i));
  }
  return(res);
}

// [[Rcpp::export]]
Rcpp::List RCRdiff_cpp(
    arma::mat pos_dat,
    arma::mat neg_dat,
    arma::mat hk_dat,
    arma::mat reg_dat,
    
    arma::vec pos_conc,
    arma::vec a_itm,
    arma::vec b_itm,
    arma::vec d_pos_itm,
    double cc_itm,
    arma::vec d_neg_itm,
    arma::vec kappa_hk_itm,
    arma::vec d_hk_itm,
    arma::vec phi_itm,
    arma::vec kappa_reg_itm,
    arma::vec d_reg_itm,
    double sigma2e_neg_itm,
    double sigma2e_phr_itm,
    double sigma2a_itm,
    double sigma2_mu_a,
    double sigma2_mu_b,
    double mu_a_itm,
    double sigma2b_itm,
    double mu_b_itm,
    vec phi_L, 
    vec phi_U,
    vec lambda_hk_itm,
    double sigma2kappa_hk_itm,
    vec alpha_r_itm,
    vec beta_r_itm,
    vec invTau2_itm,
    vec lambda_tau_itm,
    double sigma2kappa_reg_itm,
    vec lambda_hk_range_L,
    vec lambda_hk_range_U,
    vec alpha_reg_range_L,
    vec alpha_reg_range_U,
    double eta,
    double nu,
    double sigma2d_neg_itm,
    double sigma2d_phr_itm,
    double mu_a_mu,
    double mu_b_mu,
    double u,
    double v,
    
    int n_hk,
    int n_reg,
    int n_neg,
    int n_pos,
    int n_patient,
    vec Z
    ){
  Rcpp::Environment truncnorm("package:truncnorm");
  Rcpp::Function rtruncnorm = truncnorm["rtruncnorm"];
  Rcpp::Environment statmod("package:statmod");
  Rcpp::Function rinvgauss = statmod["rinvgauss"];
  Rcpp::Environment stats("package:stats");
  Rcpp::Function rgamma = stats["rgamma"];
  Rcpp::Environment matrixStats("package:matrixStats");
  Rcpp::Function colMedians = matrixStats["colMedians"];
  // Update a
  mat TempA1_1=neg_dat.each_row() - (b_itm * cc_itm).t();
  mat TempA1_2=TempA1_1.each_col() - d_neg_itm;
  
  vec A1 = colSums(TempA1_2);
  
  mat TempA2_1=pos_dat - pos_conc*b_itm.t();
  vec A2 = colSums(TempA2_1.each_col() + d_pos_itm);
  
  mat TempA3_1 = reshape(kappa_hk_itm, n_hk,n_patient);
  mat TempA3_2 = TempA3_1.each_row() + phi_itm.t();
  mat TempA3_3=hk_dat - TempA3_2.each_row() % b_itm.t();
  mat TempA3_4=TempA3_3.each_col() - d_hk_itm;
  
  vec A3 = colSums(TempA3_4);
  
  mat TempA4_1 = reshape(kappa_reg_itm, n_reg,n_patient);
  mat TempA4_2 = TempA4_1.each_row() + phi_itm.t();
  mat TempA4_3=reg_dat - TempA4_2.each_row() % b_itm.t();
  mat TempA4_4=TempA4_3.each_col() - d_reg_itm;
  
  vec A4 = colSums(TempA4_4);
  
  vec mu_a = ((A1/sigma2e_neg_itm) + (A2+A3+A4)/sigma2e_phr_itm + mu_a_itm/sigma2a_itm)/
    (n_neg/sigma2e_neg_itm + (n_pos+n_hk+n_reg)/sigma2e_phr_itm + 1/sigma2a_itm);
  double sigma_a = sqrt(1/(n_neg/sigma2e_neg_itm + (n_pos+n_hk+n_reg)/sigma2e_phr_itm + 1/sigma2a_itm));
  vec new_a_itm = randn<vec>(n_patient);
  a_itm = new_a_itm * sigma_a + mu_a;
  
  // Update b
  
  mat TempB1_1 = neg_dat.each_row() - a_itm.t();
  mat TempB1_2 = TempB1_1.each_col() - d_neg_itm;
  vec B1 = cc_itm * colSums(TempB1_2);

  mat TempB2_1=pos_dat.each_row() - a_itm.t();
  mat TempB2_2=TempB2_1.each_col() - d_pos_itm;
  mat TempB2_3=TempB2_2.each_col() % pos_conc;
  vec B2 = colSums(TempB2_3);
 
  mat TempB3_1 = TempA3_1.each_row()+phi_itm.t();
  mat TempB3_2 = hk_dat.each_row()-a_itm.t();
  mat TempB3_3 = TempB3_2.each_col()-d_hk_itm;
  
  vec B3 = colSums(TempB3_1 % TempB3_3);
  
  mat TempB4_1 = TempA4_1.each_row()+phi_itm.t();
  mat TempB4_2 = reg_dat.each_row()-a_itm.t();
  mat TempB4_3 = TempB4_2.each_col()-d_reg_itm;
  
  vec B4 = colSums(TempB4_1 % TempB4_3);
  
  mat Temp_mu_b_1 = TempB3_1 % TempB3_1;
  mat Temp_mu_b_2 = TempB4_1 % TempB4_1;
  
  vec mu_b = ((B1/sigma2e_neg_itm) + (B2+B3+B4)/sigma2e_phr_itm + mu_b_itm/sigma2b_itm)/
    (n_neg*(cc_itm*cc_itm)/sigma2e_neg_itm + (sum(pos_conc % pos_conc)+
      colSums(Temp_mu_b_1)+colSums(Temp_mu_b_2))/sigma2e_phr_itm + 1/sigma2b_itm);
       
  vec sigma_b = sqrt(1/(n_neg*(cc_itm * cc_itm)/sigma2e_neg_itm + (sum(pos_conc % pos_conc)+colSums(Temp_mu_b_1)+colSums(Temp_mu_b_2))/sigma2e_phr_itm + 1/sigma2b_itm));
  vec new_b_itm = randn<vec>(n_patient);
  b_itm = new_b_itm % sigma_b + mu_b;
  // Update c
  mat Temp_c_1 = neg_dat.each_row() - a_itm.t();
  mat Temp_c_2 = Temp_c_1.each_col() - d_neg_itm;
  vec Temp_c_3 = Temp_c_2 * b_itm;
    
  double mu_c = sum(Temp_c_3)/(n_neg * sum(b_itm % b_itm));
  double sigma_c = sqrt(sigma2e_neg_itm/(n_neg * sum(b_itm % b_itm)));
  NumericVector new_cc_itm = rtruncnorm(1, -6.0, -1.0,mu_c,sigma_c);
  cc_itm = as<double>(new_cc_itm);
     
  // Update phi

  mat Temp_phi_1 = TempA3_1.each_row() % b_itm.t();
  mat Temp_phi_2 = Temp_phi_1.each_row() + a_itm.t();
  mat Temp_phi_3 = Temp_phi_2.each_col() + d_hk_itm;
  mat Temp_phi_4 = hk_dat - Temp_phi_3;
  
  mat Temp_phi_5 = TempA4_1.each_row() % b_itm.t();
  mat Temp_phi_6 = Temp_phi_5.each_row() + a_itm.t();
  mat Temp_phi_7 = Temp_phi_6.each_col() + d_reg_itm;
  mat Temp_phi_8 = reg_dat - Temp_phi_7;
  
  vec mu_phi = (colSums(Temp_phi_4)+colSums(Temp_phi_8))/((n_hk+n_reg)*b_itm);
  vec sigma_phi = sqrt(sigma2e_phr_itm/((b_itm %b_itm) * (n_hk+n_reg)));
    
  NumericVector phi_itm_n = rtruncnorm(n_patient, phi_L, phi_U,mu_phi,sigma_phi);
  phi_itm_n[n_patient-1] = -sum(phi_itm_n[Range(0, n_patient-2)]);
  phi_itm = 1.0*(phi_itm_n);
  
  // Update kappa_hk
  mat Temp_kappa_hk_0 = hk_dat.each_row() - a_itm.t();
  mat Temp_kappa_hk_1 = Temp_kappa_hk_0.each_row() - (b_itm % phi_itm).t();
  mat Temp_kappa_hk_2 = Temp_kappa_hk_1.each_col() - d_hk_itm;
  mat Temp_kappa_hk_3 = Temp_kappa_hk_2.each_row() % b_itm.t();
  mat Temp_kappa_hk_4 = Temp_kappa_hk_3/sigma2e_phr_itm;
  mat Temp_kappa_hk_5 = Temp_kappa_hk_4.each_col() + lambda_hk_itm/sigma2kappa_hk_itm;
  
  vec Temp_kappa_hk_6 = (b_itm % b_itm)/sigma2e_phr_itm + 1/sigma2kappa_hk_itm;
  mat mu_hk_mat = Temp_kappa_hk_5.each_row()/Temp_kappa_hk_6.t();
  vec mu_hk = reshape(mu_hk_mat, n_hk*n_patient,1);
  vec sigma_hk_n0 = ones(n_patient * n_hk);
  mat sigma_hk_n1 = reshape(sigma_hk_n0, n_hk,n_patient);
  vec sigma_hk_n2 = sqrt(1/((b_itm % b_itm)/sigma2e_phr_itm + 1/sigma2kappa_hk_itm));
  mat sigma_hk_n3 = sigma_hk_n1.each_row() % sigma_hk_n2.t();
  vec sigma_hk = reshape(sigma_hk_n3, n_hk*n_patient,1);
  vec new_kappa_hk_itm = randn<vec>(n_patient * n_hk);
  kappa_hk_itm = new_kappa_hk_itm % sigma_hk + mu_hk;
  
  // Update kappa_reg
  mat Temp_kappa_reg_0 = reg_dat.each_row() - a_itm.t();
  mat Temp_kappa_reg_1 = Temp_kappa_reg_0.each_row() - (b_itm % phi_itm).t();
  mat Temp_kappa_reg_2 = Temp_kappa_reg_1.each_col() - d_reg_itm;
  mat Temp_kappa_reg_3 = Temp_kappa_reg_2.each_row() % b_itm.t();
  mat Temp_kappa_reg_4 = Temp_kappa_reg_3/sigma2e_phr_itm;
  mat Temp_kappa_reg_5 = Temp_kappa_reg_4.each_col() + alpha_r_itm/sigma2kappa_reg_itm;
  Temp_kappa_reg_5 = Temp_kappa_reg_5 + beta_r_itm * Z.t();
  
  vec Temp_kappa_reg_6 = (b_itm % b_itm)/sigma2e_phr_itm + 1/sigma2kappa_reg_itm;
  mat mu_reg_mat = Temp_kappa_reg_5.each_row()/Temp_kappa_reg_6.t();
  vec mu_reg = reshape(mu_reg_mat, n_reg*n_patient,1);
  vec sigma_reg_n0 = ones(n_patient * n_reg);
  mat sigma_reg_n1 = reshape(sigma_reg_n0, n_reg,n_patient);
  vec sigma_reg_n2 = sqrt(1/((b_itm % b_itm)/sigma2e_phr_itm + 1/sigma2kappa_reg_itm));
  mat sigma_reg_n3 = sigma_reg_n1.each_row() % sigma_reg_n2.t();
  vec sigma_reg = reshape(sigma_reg_n3, n_reg*n_patient,1);
  vec new_kappa_reg_itm = randn<vec>(n_patient * n_reg);
  kappa_reg_itm = new_kappa_reg_itm % sigma_reg + mu_reg;
  
  // Update lambda_hk
  mat Temp_lambda_1 = reshape(kappa_hk_itm,n_hk,n_patient);
  vec mu_lambda = rowSums(Temp_lambda_1)/n_patient;
  double sigma_lambda = sqrt(sigma2kappa_hk_itm/n_patient);
  NumericVector lambda_hk_itm_n = rtruncnorm(n_hk, lambda_hk_range_L, lambda_hk_range_U,mu_lambda,sigma_lambda);
  lambda_hk_itm = 1.0 * lambda_hk_itm_n;
  
  // Update alpha_r
  mat Temp_alpha_1 = reshape(kappa_reg_itm,n_reg,n_patient);
  uvec id_g1 = find(Z==0);
  uvec id_g2 = find(Z==1);
  
  int n_group1 = id_g1.size();
  int n_group2 = id_g2.size();
  vec mu_alpha = rowSums(Temp_alpha_1.cols(id_g1))/n_group1;
  double sigma_alpha = sqrt(sigma2kappa_reg_itm/n_patient);
  NumericVector alpha_r_itm_n = rtruncnorm(n_reg, alpha_reg_range_L, alpha_reg_range_U,mu_alpha,sigma_alpha);
  alpha_r_itm = 1.0 * alpha_r_itm_n;
  
  mat temp_beta = zeros(100,n_reg);

  
  for (int i=0; i<100; i++) {
    vec mu_beta = (-alpha_r_itm * n_group2 + rowSums(Temp_alpha_1.cols(id_g2)))/(sum(Z%Z)+invTau2_itm);
    vec sigma_beta = sqrt(sigma2kappa_reg_itm/(sum(Z%Z)+invTau2_itm));
    vec new_beta_r_itm = randn<vec>(n_reg);
    beta_r_itm = new_beta_r_itm % sigma_beta + mu_beta;
   
    temp_beta.row(i)=beta_r_itm.t();
    
    vec mu_tau = sqrt(lambda_tau_itm * sigma2kappa_reg_itm / (beta_r_itm%beta_r_itm));
    NumericVector invTau2_itm_n = rinvgauss(n_reg,mu_tau, lambda_tau_itm);
    invTau2_itm= 1.0 * invTau2_itm_n;

    NumericVector lambda_tau_itm_n = rgamma(n_reg,1+eta,1/(invTau2_itm)/2+nu);
    lambda_tau_itm = 1.0 * lambda_tau_itm_n;

  }
  mat Temp_beta_final =  temp_beta.rows(49,99);
  
  NumericVector beta_r_itm_n = colMedians(Temp_beta_final);
  beta_r_itm = 1.0 * beta_r_itm_n;
  
  
  // Update d_neg
  mat Temp_d_neg1 = neg_dat.each_row() - (a_itm + cc_itm*b_itm).t();
  
  vec mu_d_neg = rowSums(Temp_d_neg1/sigma2e_neg_itm)/(n_patient/sigma2e_neg_itm + 1/sigma2d_neg_itm);
  double sigma_d_neg = sqrt(1/(n_patient/sigma2e_neg_itm + 1/sigma2d_neg_itm));
  vec new_d_neg_itm = randn<vec>(n_neg);
  d_neg_itm = new_d_neg_itm * sigma_d_neg + mu_d_neg;
  
  // Update d_pos
  mat Temp_d_pos1 = pos_dat - pos_conc*b_itm.t();
  mat Temp_d_pos2 = Temp_d_pos1.each_row() - a_itm.t();
  
  vec mu_d_pos = rowSums(Temp_d_pos2/sigma2e_phr_itm)/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm);
  double sigma_d_pos = sqrt(1/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm));
  vec new_d_pos_itm = randn<vec>(n_pos);
  d_pos_itm = new_d_pos_itm * sigma_d_pos + mu_d_pos;
  
  // Update d_hk
  mat Temp_d_hk_1 = reshape(kappa_hk_itm, n_hk,n_patient);
  mat Temp_d_hk_2 = Temp_d_hk_1.each_row() + phi_itm.t();
  mat Temp_d_hk_3=hk_dat - Temp_d_hk_2.each_row() % b_itm.t();
  mat Temp_d_hk_4=Temp_d_hk_3.each_row() - a_itm.t();
  
  vec mu_d_hk = rowSums(Temp_d_hk_4/sigma2e_phr_itm)/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm);
  double sigma_d_hk = sqrt(1/(n_patient/sigma2e_phr_itm + 1/sigma2d_phr_itm));
  vec new_d_hk_itm = randn<vec>(n_hk);
  d_hk_itm = new_d_hk_itm * sigma_d_hk + mu_d_hk;
  
  // Update d_reg
  
  mat Temp_d_reg_1 = reshape(kappa_reg_itm, n_reg,n_patient);
  mat Temp_d_reg_2 = Temp_d_reg_1.each_row() + phi_itm.t();
  mat Temp_d_reg_3=reg_dat - Temp_d_reg_2.each_row() % b_itm.t();
  mat Temp_d_reg_4=Temp_d_reg_3.each_row() - a_itm.t();
  mat Temp_d_reg_5 = Temp_d_reg_4.cols(id_g1);
  
  vec mu_d_reg = rowSums(Temp_d_reg_5)/(n_group1 + sigma2d_phr_itm/sigma2d_phr_itm);
  double sigma_d_reg = sqrt(1/(n_group1/sigma2e_phr_itm + 1/sigma2d_phr_itm));
  vec new_d_reg_itm = randn<vec>(n_reg);
  d_reg_itm = new_d_reg_itm * sigma_d_reg + mu_d_reg;
  //Update mu_a and mu_b
  double new_mu_a_itm = as<double>(rnorm(1));
  double mu_mu_a = (sum(a_itm)/sigma2a_itm + mu_a_mu/sigma2_mu_a)/(n_patient/sigma2a_itm + 1/sigma2_mu_a);
  double sigma_mu_a = sqrt(1/(n_patient/sigma2a_itm + 1/sigma2_mu_a));
  mu_a_itm = new_mu_a_itm* sigma_mu_a + mu_mu_a;
 
  
  double new_mu_b_itm = as<double>(rnorm(1));
  double mu_mu_b = (sum(b_itm)/sigma2b_itm + mu_b_mu/sigma2_mu_b)/(n_patient/sigma2b_itm + 1/sigma2_mu_b);
  double sigma_mu_b = sqrt(1/(n_patient/sigma2b_itm + 1/sigma2_mu_b));
  mu_b_itm = new_mu_b_itm* sigma_mu_b + mu_mu_b;   

  //Update sigma's
  mat Temp_sigma_e_neg1 = neg_dat.each_row() - (a_itm + cc_itm*b_itm).t();
  mat Temp_sigma_e_neg2 = Temp_sigma_e_neg1.each_col() - d_neg_itm;
  vec temp_e_neg_0 = colSums(Temp_sigma_e_neg2 % Temp_sigma_e_neg2);
  double temp_e_neg_1 = sum(temp_e_neg_0)/2;
  NumericVector temp_e_neg_2 = rgamma(1, u + n_patient*n_neg/2, v+temp_e_neg_1);
  NumericVector temp_e_neg_3 = 1/temp_e_neg_2;
  sigma2e_neg_itm = as<double>(temp_e_neg_3);
  
  mat Temp_sigma_e_phr0 = pos_dat - pos_conc * b_itm.t();
  mat Temp_sigma_e_phr1 = Temp_sigma_e_phr0.each_row() - a_itm.t();
  mat Temp_sigma_e_phr2 = Temp_sigma_e_phr1.each_col() -d_pos_itm;
  vec temp_e_phr_0 =colSums(Temp_sigma_e_phr2 % Temp_sigma_e_phr2);
  double temp_e_phr_1 = (sum(temp_e_phr_0)/2);
  NumericVector temp_e_phr_2 = rgamma(1, u+n_patient*(n_pos)/2,v+temp_e_phr_1);
  NumericVector temp_e_phr_3 = 1/temp_e_phr_2;
  sigma2e_phr_itm = as<double>(temp_e_phr_3);
  
  NumericVector temp_e_a1 = rgamma(1, u + n_patient/2,  v + sum((a_itm - mu_a_itm)%(a_itm - mu_a_itm))/2);
  NumericVector temp_e_a2 = 1/temp_e_a1;
  sigma2a_itm = as<double>(temp_e_a2);
  
  NumericVector temp_e_b1 = rgamma(1, u + n_patient/2,  v + sum((b_itm - mu_b_itm)%(b_itm - mu_b_itm))/2);
  NumericVector temp_e_b2 = 1/temp_e_b1;
  sigma2b_itm = as<double>(temp_e_b2);
  
  mat Temp_e_kappa_hk_1 = reshape(kappa_hk_itm, n_hk,n_patient);
  mat Temp_e_kappa_hk_2 = Temp_e_kappa_hk_1.each_col() - lambda_hk_itm;
  vec Temp_e_kappa_hk_3 = colSums(Temp_e_kappa_hk_2 % Temp_e_kappa_hk_2);
  NumericVector sigma2kappa_hk_itm_n1 = rgamma(1, u + n_patient*n_hk/2,  v + sum(Temp_e_kappa_hk_3)/2);
  NumericVector sigma2kappa_hk_itm_n2 = 1/sigma2kappa_hk_itm_n1;
  sigma2kappa_hk_itm = as<double>(sigma2kappa_hk_itm_n2);
  
  mat Temp_e_kappa_reg_1 = reshape(kappa_reg_itm, n_reg,n_patient);
  mat Temp_e_kappa_reg_2 = Temp_e_kappa_reg_1.each_col() - alpha_r_itm;
  mat Temp_e_kappa_reg_22 = Temp_e_kappa_reg_2 - beta_r_itm * Z.t();
  vec Temp_e_kappa_reg_3 = colSums(Temp_e_kappa_reg_22 % Temp_e_kappa_reg_22);
  NumericVector sigma2kappa_reg_itm_n1 = rgamma(1, u + (n_patient*n_reg+n_reg)/2,  v + (sum(Temp_e_kappa_hk_3) + sum(beta_r_itm%beta_r_itm%invTau2_itm))/2);
  NumericVector sigma2kappa_reg_itm_n2 = 1/sigma2kappa_reg_itm_n1;
  sigma2kappa_reg_itm = as<double>(sigma2kappa_reg_itm_n2);
  
    
  NumericVector sigma2d_neg_itm1 = rgamma(1, u + n_neg/2,  v + sum(d_neg_itm%d_neg_itm)/2);
  NumericVector sigma2d_neg_itm2 = 1/sigma2d_neg_itm1;
  sigma2d_neg_itm = as<double>(sigma2d_neg_itm2);
  
  NumericVector sigma2d_phr_itm1 = rgamma(1, u + n_pos/2,  v + sum(d_pos_itm%d_pos_itm)/2);
  NumericVector sigma2d_phr_itm2 = 1/sigma2d_phr_itm1;
  sigma2d_phr_itm = as<double>(sigma2d_phr_itm2);
  
  NumericVector sigmas  {sigma2e_neg_itm,sigma2e_phr_itm,sigma2a_itm,sigma2b_itm,
                         sigma2kappa_hk_itm,sigma2kappa_reg_itm,
                         sigma2d_neg_itm,sigma2d_phr_itm};
  sigmas.names() = CharacterVector({"sigma2e_neg","sigma2e_phr",
               "sigma2a","sigma2b","sigma2kappa_hk","sigma2kappa_reg",
               "sigma2d_neg","sigma2d_phr"});
 
  return List::create(
    Named("a") = a_itm,
    Named("b") = b_itm,
    Named("c") = cc_itm,
    Named("phi") = phi_itm,
    Named("kappa_hk") =kappa_hk_itm,
    Named("kappa_reg") =kappa_reg_itm,
    Named("lambda_hk") =lambda_hk_itm,
    Named("alpha_reg") =alpha_r_itm,
    Named("beta_reg") =beta_r_itm,
    Named("invTau2") =invTau2_itm,
    Named("lambda_tau") =lambda_tau_itm,
    Named("d_neg") =d_neg_itm,
    Named("d_pos") =d_pos_itm,
    Named("d_hk") =d_hk_itm,
    Named("d_reg") =d_reg_itm,
    Named("mu_a") =mu_a_itm,
    Named("mu_b") =mu_b_itm,
    Named("sigmas") = sigmas
  );
}


