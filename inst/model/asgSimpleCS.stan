data {
  int<lower=1> n; //total number of obs over all locations and seasons
  int<lower=1> n2; //total number of obs over all locations and seasons
  int<lower=1> num[n]; //count of total patients
  int<lower=1> num2[n2]; //count of total patients
  int<lower=1> nG; //number of regions
  int<lower=1> nS; //number of seasons
  int x[n]; //week of observation indicator
  int x2[n2]; //week of observation indicator
  int<lower = 1, upper = nG> group[n]; //indicator for each group
  int<lower = 1, upper = nG> group2[n2]; //indicator for each group
  int<lower = 0> y[n]; // observed value of flu
  int<lower = 0> z[n2]; // observed value of flu
  int<lower = 1, upper = nS> seas[n]; // season indicator
  int<lower = 1, upper = nS> seas2[n2]; // season indicator
  int<lower=1> k; // number of mean parameters
  vector[k] mu0;
  cov_matrix[k] C0;
}

parameters {
  matrix[nS, k] ctheta [nG]; // Parameters for Each Year.  This is nG x nS x K array
  vector<lower=0> [k] tau; // Standard Deviation for Parameters for S seasons
  vector<lower=0> [k] tau2;
  cholesky_factor_corr[k] Lcorr;
  vector[k] mu_theta[nG]; // Means for Parameters; This is nG x K matrix
  vector[k] mu_g;
}

transformed parameters {
  real lpsi[n]; // Logit of psi
  real lpsi2[n2]; // Logit of psi
  
  for(i in 1:n){
    if (x[i] < ctheta[group[i], seas[i], 4]){
      lpsi[i] = (ctheta[group[i], seas[i], 1] + (ctheta[group[i], seas[i], 3] - ctheta[group[i], seas[i], 1])*exp(-((x[i] - ctheta[group[i], seas[i], 4])^2)/(2*(exp(ctheta[group[i], seas[i], 5]))^2)));
    }
    else{
      lpsi[i] = (ctheta[group[i], seas[i], 2] + (ctheta[group[i], seas[i], 3] - ctheta[group[i], seas[i], 2])*exp(-((x[i] - ctheta[group[i], seas[i], 4])^2)/(2*(exp(ctheta[group[i], seas[i], 6]))^2)));
    }
    
  }

  for(i in 1:n2){
    if (x2[i] < ctheta[group2[i], seas2[i], 4]){
      lpsi2[i] = (ctheta[group2[i], seas2[i], 1] + (ctheta[group2[i], seas2[i], 3] - ctheta[group2[i], seas2[i], 1])*exp(-((x2[i] - ctheta[group2[i], seas2[i], 4])^2)/(2*(exp(ctheta[group2[i], seas2[i], 5]))^2)));
    }
    else{
      lpsi2[i] = (ctheta[group2[i], seas2[i], 2] + (ctheta[group2[i], seas2[i], 3] - ctheta[group2[i], seas2[i], 2])*exp(-((x2[i] - ctheta[group2[i], seas2[i], 4])^2)/(2*(exp(ctheta[group2[i], seas2[i], 6]))^2)));
    }
    
  }

}

model {
  //Prior for Error Terms by Year
  mu_g ~ multi_normal(mu0, C0); //prior on season mean
  tau ~ student_t(4, 0, 1); //Prior on group SD
  tau2 ~ student_t(4, 0, 1); //Prior on Season SD
  Lcorr ~ lkj_corr_cholesky(1); //prior for correlations
  
  for(g in 1:nG){
    mu_theta[g] ~ multi_normal_cholesky(mu_g, diag_pre_multiply(tau, Lcorr));
    for(s in 1:nS){
      row(ctheta[g], s) ~ multi_normal_cholesky(mu_theta[g], diag_pre_multiply(tau2, Lcorr)); 
      // c(beta1[g], beta2[g], leta[g], mu[g], lsigma1[g], lsigma2[g])
    }
  }
  
  y ~ binomial_logit(num, lpsi);
  z ~ binomial_logit(num2, lpsi2);
}

generated quantities {
  vector[n] log_lik;
  matrix[k,k] Omega;
  Omega = multiply_lower_tri_self_transpose(Lcorr);
  for(i in 1:n){
    log_lik[i] = binomial_logit_lpmf(y[i] | num[i], lpsi[i]);
  }
}

