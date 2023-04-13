#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

int Sign(double x)
{
  if (x < 0) {
    return -1;
  } else if (x > 0) {
    return 1;
  } else {
    return 0;
  }
}

double morethan1(double x)
{
  if (x > 1) {
    return 1;
  }else {
    return x;
  }
}

arma::mat MGF_cpp(double t, 
                  const arma::vec& Resid_unrelated_outliers, 
                  const Rcpp::List& TwoSubj_list,
                  const Rcpp::List& update_ThreeSubj_list,
                  double MAF)
{
  // Unrelated subjects.
  arma::vec lambda = arma::exp(t * Resid_unrelated_outliers);
  
  arma::vec alpha = 1 - MAF + MAF * lambda; 
  arma::vec alpha_1 = MAF * Resid_unrelated_outliers % lambda; 
  arma::vec alpha_2 = Resid_unrelated_outliers % alpha_1;
  
  arma::vec M_G0_all = alpha % alpha;
  arma::vec M_G1_all = 2 * alpha % alpha_1;
  arma::vec M_G2_all = 2 * (alpha_1 % alpha_1 + alpha % alpha_2);
  
  // Two related subjects in a family.
  int n1 = TwoSubj_list.length();
  if (n1 != 0)
  {
    for (int i = 0; i < n1; i++)
    {
      Rcpp::List TwoSubj_list_temp = TwoSubj_list[i];
      arma::vec Resid = Rcpp::as<arma::vec>(TwoSubj_list_temp["Resid"]);
      arma::vec Rho = Rcpp::as<arma::vec>(TwoSubj_list_temp["Rho"]);

      arma::vec temp = (1 - Rho) * MAF * (1 - MAF);

      double R1 = Resid[0];
      double R2 = Resid[1];
      double Rsum = R1 + R2;
      double etR1 = exp(t * R1);
      double etR2 = exp(t * R2);

      arma::vec midterm1 = etR1 * temp;
      arma::vec midterm2 = etR2 * temp;
      arma::vec midterm3 = etR1 * etR2 * (MAF - temp);

      arma::vec M_G0 = midterm1 + midterm2 + midterm3 - temp + 1 - MAF;
      arma::vec M_G1 = R1 * midterm1 + R2 * midterm2 + Rsum * midterm3;
      arma::vec M_G2 = R1*R1 * midterm1 + R2*R2 * midterm2 + Rsum*Rsum * midterm3;

      M_G0_all = arma::join_cols(M_G0_all, M_G0);
      M_G1_all = arma::join_cols(M_G1_all, M_G1);
      M_G2_all = arma::join_cols(M_G2_all, M_G2);
    }
  }

  // Three above Related Subjects.
  int n2 = update_ThreeSubj_list.length();
  if (n2 != 0)
  {
    for (int i = 0; i < n2; i++)
    {
      Rcpp::List ThreeSubj_list_temp = update_ThreeSubj_list[i];
      arma::vec stand_S = as<arma::vec>(ThreeSubj_list_temp["stand.S"]);
      arma::vec arr_prob = as<arma::vec>(ThreeSubj_list_temp["arr.prob"]);
      
      arma::vec midterm0 = exp(t * stand_S) % arr_prob;
      arma::vec midterm1 = stand_S % midterm0;
      arma::vec midterm2 = stand_S % midterm1;
      
      M_G0_all = arma::join_cols(M_G0_all, arma::vec{arma::accu(midterm0)});
      M_G1_all = arma::join_cols(M_G1_all, arma::vec{arma::accu(midterm1)});
      M_G2_all = arma::join_cols(M_G2_all, arma::vec{arma::accu(midterm2)});
    }
  }
  
  return arma::join_rows(M_G0_all, M_G1_all, M_G2_all);
}

Rcpp::List fastgetroot_cpp(double sum_R_nonOutlier,
                           double R_GRM_R_nonOutlier,
                           const arma::vec& Resid_unrelated_outliers,
                           const Rcpp::List& TwoSubj_list,
                           const Rcpp::List& update_ThreeSubj_list,
                           double Score,
                           double MAF,
                           double init_t,
                           double tol,
                           int maxiter = 100)
{
  double t = init_t;
  arma::vec MGF0; arma::vec MGF1; arma::vec MGF2;
  double CGF1 = 0; double CGF2 = 0;
  double diff_t = R_PosInf;
  bool converge = true;
  int iter;
  
  double mean = 2 * MAF * sum_R_nonOutlier;
  double var = 2 * MAF * (1 - MAF) * R_GRM_R_nonOutlier;

  for (iter = 1; iter < maxiter; iter++)
  {
    double old_t = t;
    double old_diff_t = diff_t;
    double old_CGF1 = CGF1;

    arma::mat MGF_all = MGF_cpp(t, Resid_unrelated_outliers, TwoSubj_list, update_ThreeSubj_list, MAF);
    
    MGF0 = MGF_all.col(0);
    MGF1 = MGF_all.col(1);
    MGF2 = MGF_all.col(2);
    
    arma::vec temp = MGF1 / MGF0;
    CGF1 = arma::accu(temp) + mean + var * t - Score;
    CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
    
    diff_t = - CGF1/CGF2;
    
    if (std::isnan(diff_t))
    {
      t = t / 2;
      diff_t = std::min(std::abs(t), 1.0) * Sign(Score);
      continue;
    }
    
    if (std::isnan(old_CGF1) || Sign(CGF1) != Sign(old_CGF1)) 
    {
      if (std::abs(diff_t) < tol) 
      {
        t = old_t + diff_t;
        break;
      } else {
        while (std::abs(old_diff_t) > tol && std::abs(diff_t) > std::abs(old_diff_t) - tol) 
        {
          diff_t = diff_t / 2;
        }
        t = old_t + diff_t;
        continue;
      }
    }
    
    if (Sign(old_t) != Sign(old_t + diff_t) && Sign(CGF1) == Sign(old_CGF1)) 
    {
      while (Sign(old_t) != Sign(old_t + diff_t)) 
      {
        diff_t = diff_t / 2;
      }
      t = old_t + diff_t;
      continue;
    }
    
    t = old_t + diff_t;
    if (std::abs(diff_t) < tol) break;
  }
  
  if (iter == maxiter) converge = false;
  
  return Rcpp::List::create(Named("root") = t,
                            Named("iter") = iter,
                            Named("converge") = converge);
}

Rcpp::List GetProb_cpp(double sum_R_nonOutlier,
                       double R_GRM_R_nonOutlier,
                       const arma::vec& Resid_unrelated_outliers,
                       const Rcpp::List& TwoSubj_list,
                       const Rcpp::List& update_ThreeSubj_list,
                       double Score,
                       double MAF,
                       bool lower_tail,
                       double zeta,
                       double tol)
{
  Rcpp::List out = fastgetroot_cpp(sum_R_nonOutlier, R_GRM_R_nonOutlier, Resid_unrelated_outliers,
                                   TwoSubj_list, update_ThreeSubj_list, Score, MAF, zeta, tol);
  
  zeta = as<double>(out["root"]);
  int iter = as<int>(out["iter"]);
  
  arma::mat MGF_all = MGF_cpp(zeta, Resid_unrelated_outliers, TwoSubj_list, update_ThreeSubj_list, MAF);
  
  arma::vec MGF0 = MGF_all.col(0);
  arma::vec MGF1 = MGF_all.col(1);
  arma::vec MGF2 = MGF_all.col(2);
  
  double mean = 2 * MAF * sum_R_nonOutlier;
  double var = 2 * MAF * (1 - MAF) * R_GRM_R_nonOutlier;
  
  arma::vec temp = MGF1 / MGF0;
  double CGF0 = arma::accu(log(MGF0)) + mean * zeta + 0.5 * var * zeta * zeta;
  double CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
  
  double w = Sign(zeta) * sqrt(2 * (zeta * Score - CGF0));
  double v = zeta * sqrt(CGF2);
  
  double u = w + 1/w * log(v/w);
  double pval = R::pnorm(u, 0, 1, lower_tail, false);
  
  return Rcpp::List::create(Named("pval") = pval,
                            Named("iter") = iter,
                            Named("zeta") = zeta);
}

// [[Rcpp::export]]
Rcpp::List GetTwoProb_cpp(int order_interval,
                          double MAF_ratio,
                          double R_GRM_R_TwoSubjOutlier,
                          double g_var,
                          double S_var_GRM,
                          double sum_R_nonOutlier,
                          double R_GRM_R_nonOutlier,
                          const arma::vec& Resid_unrelated_outliers,
                          const Rcpp::List& TwoSubj_list,
                          Rcpp::List ThreeSubj_list,
                          double Score,
                          double MAF,
                          double zeta,
                          double tol)
{
  double Var_ThreeOutlier = 0;
  bool lower_tail;
  double zeta1; double zeta2;
  
  int n1 = ThreeSubj_list.length();
  Rcpp::List update_ThreeSubj_list(n1);
  
  if (n1 != 0)
  {
    for (int i = 0; i < n1; i++)
    {
      Rcpp::List ThreeSubj_list_temp = ThreeSubj_list[i];
      
      Rcpp::List CLT_temp =  ThreeSubj_list_temp["CLT"];
      arma::vec CLT_temp1 = as<arma::vec>(CLT_temp[order_interval - 1]);
      arma::vec CLT_temp2 = as<arma::vec>(CLT_temp[order_interval]);
      
      arma::vec arr_prob = MAF_ratio * CLT_temp1 + (1 - MAF_ratio) * CLT_temp2;
      arma::vec stand_S = as<arma::vec>(ThreeSubj_list_temp["stand.S"]);
      
      update_ThreeSubj_list[i] = Rcpp::List::create(Named("stand.S") = stand_S,
                                                    Named("arr.prob") = arr_prob);
      
      arma::vec temp1 = stand_S % arr_prob;
      double temp2 = arma::accu(temp1);
      double Var_ThreeOutlier_temp = arma::accu(stand_S % temp1) - temp2 * temp2;
      Var_ThreeOutlier = Var_ThreeOutlier + Var_ThreeOutlier_temp;
    }
  }
  
  double Var_nonOutlier = g_var * R_GRM_R_nonOutlier;
  double Var_unrelated_outliers = g_var * arma::accu(Resid_unrelated_outliers % Resid_unrelated_outliers);
  double Var_TwoOutlier = g_var * R_GRM_R_TwoSubjOutlier;
  
  double EmpVar = Var_nonOutlier + Var_unrelated_outliers + Var_TwoOutlier + Var_ThreeOutlier;
  double Var_Ratio = S_var_GRM / EmpVar;
  double Score_adjust = Score / sqrt(Var_Ratio);
  
  zeta1 = morethan1(abs(Score_adjust) / S_var_GRM); zeta2 = - abs(zeta);
  
  Rcpp::List output1 = GetProb_cpp(sum_R_nonOutlier, R_GRM_R_nonOutlier, Resid_unrelated_outliers, TwoSubj_list, 
                                   update_ThreeSubj_list, abs(Score_adjust), MAF, lower_tail = false, zeta1, 1e-4);
  
  Rcpp::List output2 = GetProb_cpp(sum_R_nonOutlier, R_GRM_R_nonOutlier, Resid_unrelated_outliers, TwoSubj_list, 
                                   update_ThreeSubj_list, -abs(Score_adjust), MAF, lower_tail = true, zeta2, tol);
  
  double pval1 = as<double>(output1["pval"]); double pval2 = as<double>(output2["pval"]);
  
  return List::create(Named("pval1") = pval1, Named("pval2") = pval2);
}