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
matrix[T,N] y_data_indicator_missing;

}
parameters {
matrix[N,P*L+L-1] beta;
matrix[T,N] y_data_missing;
matrix[T,N] theta_innovation;
matrix[T,N] error_innovation;
corr_matrix[P*L+L-1] Omega; // prior correlation
vector<lower=0>[P*L+L-1] tau; // prior scale
vector[N] mu_theta;

vector[P*L+L-1] mu_beta; 
real<lower=0> sigmasq_theta[N];
real<lower=0> sigmasq_y[N];
vector<lower=0>[P*L+L-1] sigmasq_beta;

}

transformed parameters {
matrix[N,T] theta;
matrix[N,T] lagged_theta;
real<lower=0> sigma_y[N];      
real<lower=0> sigma_theta[N];    
vector<lower=0>[P*L+L-1] sigma_beta;
matrix[N,T] explanvars;  
matrix[T,N] y_data_observed;



for (n in 1:N)
for (t in 1:T)
y_data_observed[t,n] <-y_data[t,n];


for (n in 1:N)
for (t in 1:T){
if(y_data_indicator_missing[t,n])
y_data_observed[t,n] <-y_data_missing[t,n];
}


for(n in 1:N)
for(t in 1:T) explanvars[n,t]<-0;

for(n in 1:N)
for(t in (1+L):T)
for(p in 1:P)
for(l in 1:L)
explanvars[n,t]<-beta[n,p +(l-1)*P]*x_data[t-(l-1),p,n] + explanvars[n,t];


for(n in 1:N)
for(t in (1+L):T)
for(l in 1:(L-1))
lagged_theta[n,t]<-beta[n,P*L+l]*theta[n,t-l] + lagged_theta[n,t];



for(n in 1:N) sigma_y[n] <- sqrt(sigmasq_y[n]);
for(p in 1:(P*L+L-1))sigma_beta[p] <- (sigmasq_beta[p]);
for(n in 1:N) sigma_theta[n] <- sqrt(sigmasq_theta[n]);


for(n in 1:N) for(t in 1:(L+1)) theta[n,t]<-y_data[t,n];

for (n in 1:N)
for (t in (2+L):T){
theta[n,t]<-lagged_theta[n,t] + explanvars[n,t] +sigma_theta[n]*theta_innovation[t,n];
}

}

model {
matrix[L*P+L-1,L*P+L-1] Sigma_beta;

mu_beta ~ normal(0, 1);
Sigma_beta <- quad_form_diag(Omega,tau);
tau ~ cauchy(0,2.5);
\
Omega ~ lkj_corr(5);
for(n in 1:N) beta[n]~multi_normal(mu_beta,Sigma_beta);
sigmasq_y ~ inv_gamma(0.001, 0.001);
sigmasq_beta ~ inv_gamma(0.001, 0.001);
sigmasq_theta ~inv_gamma(0.001, 0.001);
for(n in 1:N) theta_innovation[n]~normal(0,1);

for (n in 1:N)
for (t in (2+L):T){
y_data_observed[t,n]~normal(theta[n,t] , sigma_y);
}

}


"
fit <- stan(model_code = model, data = stan.dat, iter = 100, chains = 1, verbose = TRUE)
