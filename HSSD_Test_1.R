library(rstan)
load("c:/users/sam/desktop/rstan-hierarchical-state-space-dynamics/testdata.rdata")



model="
data {
  int<lower=0> P;
  int<lower=0> N;
  int<lower=0> T;
  int<lower=0> L;
  matrix[P,N] x_data[T];
  matrix[T,N] y_data;

}

parameters {
  matrix[N,P*L+L-1] beta;
  matrix[T,N] theta_innovation;
  corr_matrix[P*L+L-1] Omega;
  vector<lower=0>[P*L+L-1] tau;
  vector[N] mu_theta;
  vector[P*L+L-1] mu_beta; 
  real<lower=0> sigmasq_theta[N];
  real<lower=0> sigmasq_y[N];
  vector<lower=0>[P*L+L-1] sigmasq_beta;

}
transformed parameters {
  matrix[N,T] theta;
  matrix[N,T] explanvars;  


  for(n in 1:N)
    for(t in 1:T) explanvars[n,t]<-0;

  for(n in 1:N)
    for(t in (1+L):T)
      for(p in 1:P)
        for(l in 1:L)
          explanvars[n,t]<-beta[n,p +(l-1)*P]*x_data[t-(l-1),p,n] + explanvars[n,t];



  for(n in 1:N) for(t in 1:(L+1)) theta[n,t]<-y_data[t,n];


  for (n in 1:N)
    for (t in (2+L):T){
      theta[n,t]<- explanvars[n,t] +sigmasq_theta[n]*theta_innovation[t,n];
    }
  

  for (n in 1:N)
    for (t in (2+L):T)
      for(l in 1:(L-1))
        theta[n,t]<-beta[n,P*L+l]*theta[n,t-l]+theta[n,t];
    
}

model {
  matrix[L*P+L-1,L*P+L-1] Sigma_beta;
  mu_beta ~ normal(0, 1);
  Sigma_beta <- quad_form_diag(Omega,tau);
  tau ~ cauchy(0,2.5);
  Omega ~ lkj_corr(1);
  for(n in 1:N) beta[n]~multi_normal(mu_beta,Sigma_beta);

  sigmasq_y ~ inv_gamma(0.001, 0.001);
  sigmasq_theta ~inv_gamma(0.001, 0.001);

  for(t in 1:T) theta_innovation[t]~normal(0,1);

  for (n in 1:N)
    for (t in (2+L):T){
      y_data[t,n]~normal(theta[n,t] , sigmasq_y);
  }

}


"
fit <- stan(model_code = model, data = stan.dat, iter = 1000, chains = 4, verbose = TRUE)
mat=as.matrix(fit)
