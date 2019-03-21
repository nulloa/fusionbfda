data {
  int<lower=1> N_obs; //total number of obs
  int<lower=1> N_subj; //total number of season-region combination
  int<lower=1> N_weeks; //total number of weeks ie obs per season-region
  int<lower=1> dim_space; // number of eigen functions
  real W[N_obs]; // observed value of flu
  int week[N_obs]; // indicator for the week observed
  int subj[N_obs]; // indicator for the observed subject
  real E[N_weeks, dim_space]; // observed value of flu
}

parameters {
  real<lower=0> eps; // SD for data model
  vector[dim_space] beta[N_subj];
}

transformed parameters {
  vector[N_weeks] m[N_subj]; // Underlying smooth function
  for(i in 1:N_subj){
    for(t in 1:N_weeks){
      m[i,t] = beta[i,1]*E[t,1]+
	         beta[i,2]*E[t,2]+
	         beta[i,3]*E[t,3]+
	         beta[i,4]*E[t,4]+
	         beta[i,5]*E[t,5]+
	         beta[i,6]*E[t,6]+
	         beta[i,7]*E[t,7]+
	         beta[i,8]*E[t,8]+
	         beta[i,9]*E[t,9]+
	         beta[i,10]*E[t,10]+
	         beta[i,11]*E[t,11]+
	         beta[i,12]*E[t,12]+
	         beta[i,13]*E[t,13]+
	         beta[i,14]*E[t,14]+
	         beta[i,15]*E[t,15]+
	         beta[i,16]*E[t,16]+
	         beta[i,17]*E[t,17]+
	         beta[i,18]*E[t,18]+
	         beta[i,19]*E[t,19]+
	         beta[i,20]*E[t,20]+
	         beta[i,21]*E[t,21]+
	         beta[i,22]*E[t,22]+
	         beta[i,23]*E[t,23]+
	         beta[i,24]*E[t,24]+
	         beta[i,25]*E[t,25]+
	         beta[i,26]*E[t,26]+
	         beta[i,27]*E[t,27]+
	         beta[i,28]*E[t,28]+
	         beta[i,29]*E[t,29]+
	         beta[i,30]*E[t,20];
    }
  }
}

model {

  eps ~ cauchy(0,1); //Prior on model SD
  
  
  for(i in 1:N_obs){
    W[i]~normal(m[subj[i], week[i]], eps);
  }
  
  for(i in 1:N_subj){ 
    for (j in 1:dim_space){
      beta[i,j] ~ normal(0, 10);
	  }
  }
  
}

// generated quantities {
//   vector[N_obs*N_subj] log_lik;
//   for(i in 1:N_subj){
//     for(t in 1:N_obs){
//       log_lik[i] = normal_lpmf(W[i,t] | m[i,t], eps);
//     }
//   }
// }

