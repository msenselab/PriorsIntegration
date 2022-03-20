//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// linear encoding model

functions {
real get_weight(real sig_t, real sig_p){
  return(sig_t /(sig_t + sig_p));
}

real get_mu(real wp_p, real mu_p, real x){  
  return(wp_p *mu_p +(1-wp_p) * x);
}

real get_variance(real sig_t, real sig_p){
  return(sig_t*sig_p/(sig_t + sig_p));
}

matrix predictor_mix_rng(real[] x, real sig_t, real sig_pr2_s, real sig_pr2_l, real mu_p_s, real mu_p_l, real sig2_mn) {
vector[num_elements(x)] predY;         //predication of RP generated by model
vector[num_elements(x)] wp_local;
vector[num_elements(x)] mu_r;
vector[num_elements(x)] sig_wm2;
vector[num_elements(x)] sig_r;
vector[num_elements(x)] log_lik;
vector[num_elements(x)] wp;
vector[num_elements(x)] sig_dL_sq;

for (m in 1:num_elements(x)) {
  //part 1 integration of local priors firstly
  if (x[m] < 1){
    wp_local[m] = get_weight((x[m]*sig_t)^2, sig_pr2_s); //weight of local short prior
    mu_r[m] = get_mu(wp_local[m], mu_p_s, x[m]);
    sig_wm2[m] = get_variance((x[m]*sig_t)^2, sig_pr2_s);
  }else{
    wp_local[m] = get_weight((x[m]*sig_t)^2, sig_pr2_l); //weight of local long prior
    mu_r[m] = get_mu(wp_local[m], mu_p_l, x[m]);
    sig_wm2[m] = get_variance((x[m]*sig_t)^2, sig_pr2_l);
  }
  wp[m] = 0;
  sig_dL_sq[m] = 0;
  sig_r[m] = sqrt(sig_wm2[m] +sig2_mn);
  predY[m] = normal_rng(mu_r[m], sig_r[m]); 
  log_lik[m] = normal_lpdf(predY[m]|mu_r[m], sig_r[m]);
  }
  return(append_col(append_col(append_col(wp, mu_r), append_col(sig_r, predY)), append_col(append_col(wp_local, sig_dL_sq), log_lik)));
}

matrix predictor_rng(real[] x, real sig_t, real sig_pr2, real mu_p, real sig_mn2) {
vector[num_elements(x)] predY;         
vector[num_elements(x)] sig2_sm;
vector[num_elements(x)] wp_pr;
vector[num_elements(x)] mu_r;
vector[num_elements(x)] sig_r;


for (i in 1:num_elements(x))
{
  sig2_sm[i] = get_variance((x[i]*sig_t)^2, sig_pr2);
  wp_pr[i] = get_weight(sig_t, sig_pr2); 
  mu_r[i] = get_mu(wp_pr[i], mu_p, x[i]);
  sig_r[i] =  sqrt(sig2_sm[i] +sig_mn2);
  predY[i] = normal_rng(mu_r[i], sig_r[i]); 
}
return(append_col(append_col(wp_pr, mu_r), append_col(sig_r, predY)));
}

}

data {
int<lower=0> n_s;  //number of data points in the short session 
int<lower=0> n_l;  //number of data points in the long session
int<lower=0> n_mix;  //number of data points in the mix session 
real<lower=0> Y_s[n_s];   //measured reproductive duration (short group)
real<lower=0> X_s[n_s];   //stimulus duration (short group)
real<lower=0> Y_l[n_l];   //measured reproductive duration (long group)
real<lower=0> X_l[n_l];   //stimulus duration (long group)
real<lower=0> X_mix[n_mix];   //stimulus duration (mix group)
real<lower=0> xsnew[41];  //new target duration for short group
real<lower=0> xlnew[121];  //new target duration for long group
real<lower=0> xnew[162];  //new target duration for mixed group
}

parameters {
//hyperparameters
real<lower=0.2, upper =1>  mu_p_s;   // mean of short prior of short session(linear scale)
real<lower=1, upper =2.4> mu_p_l;   // mean of long prior of long session(linear scale)
real<lower=0, upper=2> sig2_mn;
real<lower=0, upper=2> sig_t;   // variance of D_s distribution in short and long group
real<lower=0, upper=2> sig_pr2_s;  // sigma^2 of prior in linear scale
real<lower=0, upper=2> sig_pr2_l;  // sigma^2 of prior in linear scale
}


model {
vector[n_s] sig2_sm_s;
vector[n_s] wp_short;
vector[n_s] mu_r_s;
vector[n_s] sig_r_s;
vector[n_l] wp_long;
vector[n_l] mu_r_l;
vector[n_l] sig_r_l;
vector[n_l] sig2_sm_l;
sig2_mn ~ cauchy(0, 1);
sig_t ~ cauchy(0, 1);
sig_pr2_s ~ cauchy(0, 1);
sig_pr2_l ~ cauchy(0, 1);
mu_p_s ~ normal(1, 1);
mu_p_l ~ normal(1, 1);


//short session
for (i in 1:n_s)
{
  sig2_sm_s[i] = get_variance((X_s[i]*sig_t)^2, sig_pr2_s);
  wp_short[i] = get_weight((X_s[i]*sig_t)^2, sig_pr2_s); //weight of short prior
  mu_r_s[i] = get_mu(wp_short[i], mu_p_s, X_s[i]); 
  sig_r_s[i] = sqrt(sig2_sm_s[i] +sig2_mn);
  Y_s[i] ~ normal(mu_r_s[i], sig_r_s[i]);
}


//long session
for (i in 1:n_l)
{
  sig2_sm_l[i] = get_variance((X_l[i]*sig_t)^2, sig_pr2_l);
  wp_long[i] = get_weight((X_l[i]*sig_t)^2, sig_pr2_l); //weight of long prior
  mu_r_l[i] = get_mu(wp_long[i], mu_p_l, X_l[i]);
  sig_r_l[i] = sqrt(sig2_sm_l[i] +sig2_mn);
  Y_l[i] ~ normal(mu_r_l[i], sig_r_l[i]);
}

}

generated quantities {
matrix[n_s, 4] ynew_s;
matrix[n_l, 4] ynew_l;
matrix[41, 4] ypred_s;
matrix[121, 4] ypred_l;
matrix[n_mix, 7] ynew_mix;
matrix[162, 7] ypred_mix;
vector[n_mix] log_lik;
real log_lik_sum;
ynew_mix = predictor_mix_rng(X_mix, sig_t, sig_pr2_s, sig_pr2_l, mu_p_s, mu_p_l, sig2_mn);
ypred_mix = predictor_mix_rng(xnew, sig_t, sig_pr2_s, sig_pr2_l, mu_p_s, mu_p_l, sig2_mn);
log_lik =  col(ynew_mix, 7);
log_lik_sum = sum(log_lik);
ynew_s = predictor_rng(X_s, sig_t, sig_pr2_s, mu_p_s, sig2_mn);
ynew_l = predictor_rng(X_l, sig_t, sig_pr2_l, mu_p_l, sig2_mn);
ypred_s = predictor_rng(xsnew, sig_t, sig_pr2_s, mu_p_s, sig2_mn);
ypred_l = predictor_rng(xlnew, sig_t, sig_pr2_l, mu_p_l, sig2_mn);
}
