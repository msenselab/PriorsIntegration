//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//

//assuming no global prior and local prior worked independently
functions {
  real mu_lognormal(real mu, real sigma){
    return(log(mu) - 0.5 *log(1+sigma/mu^2));
  }
  
  real variance_lognormal(real mu, real sigma){
    return(log(1+sigma/mu^2));
  }
  
  real mu_normal(real mu, real sigma){
    return(exp(mu+ sigma *0.5));
  }
  
  real variance_normal(real mu, real sigma){
    return((exp(sigma)-1)*exp(2*mu + sigma));
  }
  
  real get_weight(real sig_t, real sig_p){
    return(sig_t /(sig_t + sig_p));
  }
  
  real get_mu(real wp_p, real mu_p, real x){  
    return(wp_p *mu_p +(1-wp_p) * x);
  }
  
  real get_variance(real sig_t, real sig_p){
    return(sig_t*sig_p/(sig_t + sig_p));
  }
  
  
  matrix predictor_mix_rng(real[] x, real sig2_t, real sig_pr2_s_log, real sig_pr2_l_log,real mu_p_s_log, real mu_p_l_log, real sig2_mn) {
    vector[num_elements(x)] predY;         //predication of RP generated by model
    vector[num_elements(x)] mu_r_log;
    vector[num_elements(x)] mu_r;
    vector[num_elements(x)] sig_r_log;
    vector[num_elements(x)] sig_r;
    vector[num_elements(x)] log_lik;
    vector[num_elements(x)] wp;
    vector[num_elements(x)] sig_dL_sq;
    vector[num_elements(x)] wp_local;
    
    for (m in 1:num_elements(x)) {
      //part 1 integration of local priors firstly
      if (x[m] < 1){
        wp_local[m] = get_weight(sig2_t, sig_pr2_s_log);
        mu_r_log[m] = get_mu(wp_local[m], mu_p_s_log, log(x[m]));
        sig_r_log[m] =  get_variance(sig2_t, sig_pr2_s_log);  
      }else{
        wp_local[m] = get_weight(sig2_t, sig_pr2_l_log);
        mu_r_log[m] = get_mu(wp_local[m], mu_p_l_log, log(x[m]));
        sig_r_log[m] = get_variance(sig2_t, sig_pr2_l_log);
      }
      wp[m] = 0;
      sig_dL_sq[m] = 0;
      mu_r[m] = mu_normal(mu_r_log[m], sig_r_log[m]);
      sig_r[m] =sqrt(variance_normal(mu_r_log[m], sig_r_log[m]) +sig2_mn);
      predY[m] = normal_rng(mu_r[m], sig_r[m]); 
      log_lik[m] = normal_lpdf(predY[m]|mu_r[m], sig_r[m]);
    }
    return(append_col(append_col(append_col(wp, mu_r), append_col(sig_r, predY)), append_col(append_col(wp_local, sig_dL_sq), log_lik)));
  }
  
}

data {
  int<lower=0> n_mix;  //number of data points in the mix session 
  real<lower=0> X_mix[n_mix];   //stimulus duration (mix group)
  real<lower=0> Y_mix[n_mix];   //measured reproductive duration (mixed) 
  real<lower=0> xnew[162];  //new target duration for mixed group
  real<lower=0, upper=4> sig2_mn;
  real<lower=0, upper=4> sig2_t; 
}

parameters {
  real<lower=0, upper = 1> mu_p_s;
  real<lower=0, upper =2> sig_pr2_s; 
  real<lower=1, upper =2.4> mu_p_l;
  real<lower=0, upper =2> sig_pr2_l;
}

transformed parameters { 
  real mu_p_s_log;  // mean of global prior in log scale
  real<lower=0> sig_pr2_s_log;  //variance of global prior in log scale
  real mu_p_l_log;  // mean of global prior in log scale
  real<lower=0> sig_pr2_l_log;  //variance of global prior in log scale
  mu_p_s_log = log(mu_p_s) - 0.5 *log(1+sig_pr2_s/mu_p_s^2);
  sig_pr2_s_log = log(1+sig_pr2_s/mu_p_s^2);
  mu_p_l_log = log(mu_p_l) - 0.5 *log(1+sig_pr2_l/mu_p_l^2);
  sig_pr2_l_log = log(1+sig_pr2_l/mu_p_l^2);
} 


model {
  vector[n_mix] sig_sm2;
  vector[n_mix] sig_dL_sq;
  vector[n_mix] mu_dL;
  vector[n_mix] mu_r;
  vector[n_mix] sig_r;
  vector[n_mix] wp_local;
  
  //mixed session
  for (m in 1:n_mix) {
    //integration of local priors firstly
    if (X_mix[m] < 1){
      wp_local[m] = get_weight(sig2_t, sig_pr2_s_log);
      mu_dL[m] = get_mu(wp_local[m], mu_p_s_log, log(X_mix[m]));
      sig_dL_sq[m] = get_variance(sig2_t, sig_pr2_s_log);
    }else{
      wp_local[m] = get_weight(sig2_t, sig_pr2_l_log);
      mu_dL[m] = get_mu(wp_local[m], mu_p_l_log, log(X_mix[m]));
      sig_dL_sq[m] = get_variance(sig2_t, sig_pr2_l_log);
    }
    
    mu_r[m] = mu_normal(mu_dL[m], sig_dL_sq[m]);
    sig_r[m] = sqrt(variance_normal(mu_dL[m], sig_dL_sq[m]) +sig2_mn);
    Y_mix[m] ~ normal(mu_r[m], sig_r[m]);  
  }
}

generated quantities {
  matrix[n_mix, 7] ypred_mix;
  matrix[162, 7] ynew_mix;
  vector[n_mix] log_lik;
  real log_lik_sum;
  ynew_mix = predictor_mix_rng(xnew, sig2_t, sig_pr2_s_log, sig_pr2_l_log, mu_p_s_log, mu_p_l_log, sig2_mn);
  ypred_mix = predictor_mix_rng(X_mix, sig2_t, sig_pr2_s_log, sig_pr2_l_log, mu_p_s_log, mu_p_l_log, sig2_mn);
  log_lik = col(ypred_mix,7);
  log_lik_sum = sum(log_lik);
}