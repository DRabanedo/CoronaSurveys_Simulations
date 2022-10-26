### Zheng et al. (2005) overdispersed model
###########################################
library(rstan)

overdispersed_model = "
data {             
  int<lower=0> I;                        // respondents
  int<lower=0> K;                        // subpopulations
  vector[K] mu_beta;                     // prior mean of beta
  vector<lower=0>[K] sigma_beta;         // prior variance of beta
  int  y[I,K];                           // known by respondent i in subpopulation k
  }

parameters {
  vector[I] alpha;                       // log degree
  vector[K] beta;                        // log prevalence of group in population
  vector<lower = 0 , upper = 1>[K] inv_omega;  // ineverse overdispersion; implies the uniform prior 
  real mu_alpha;                         // prior mean for alpha
  real<lower=0> sigma_alpha;             // prior scale for alpha
  }

model {
// priors
  alpha ~ normal(mu_alpha, sigma_alpha);  
  beta ~ normal(mu_beta, sigma_beta);     // informative prior on beta: location and scale are identified             

// hyperpriors
  //mu_alpha ~ normal(0,25);                // weakly informative (no prior in paper)
  //sigma_alpha ~ normal(0,5);              // weakly informative (no prior in paper)


  for (k in 1:K) {
    real omega_k_m1;
    omega_k_m1 = inv(inv(inv_omega[k]) - 1) ;
    for (i in 1:I) {
      real xi_i_k;
      xi_i_k = omega_k_m1 * exp(alpha[i] + beta[k])  ;
      y[i,k] ~ neg_binomial(xi_i_k, omega_k_m1);             
      }
    }
  }"

overdispersed_model2 = "
data {             
  int<lower=0> I;                        // respondents
  int<lower=0> K;                        // subpopulations
  vector[K] mu_beta;                     // prior mean of beta
  vector[I] mu_alpha;                    // prior mean of alpha
  vector<lower=0>[K] sigma_beta;         // prior variance of beta
  int  y[I,K];                           // known by respondent i in subpopulation k
  }

parameters {
  vector[I] alpha;                       // log degree
  vector[K] beta;                        // log prevalence of group in population
  vector<lower = 0 , upper = 1>[K] inv_omega;  // ineverse overdispersion; implies the uniform prior 
  //real mu_alpha;                         // prior mean for alpha
  real<lower=0> sigma_alpha;             // prior scale for alpha
  }

model {
// priors
  alpha ~ normal(mu_alpha, sigma_alpha);  
  beta ~ normal(mu_beta, sigma_beta);     // informative prior on beta: location and scale are identified             

// hyperpriors
  //mu_alpha ~ normal(0,25);                // weakly informative (no prior in paper)
  sigma_alpha ~ normal(0,5);              // weakly informative (no prior in paper)


  for (k in 1:K) {
    real omega_k_m1;
    omega_k_m1 = inv(inv(inv_omega[k]) - 1) ;
    for (i in 1:I) {
      real xi_i_k;
      xi_i_k = omega_k_m1 * exp(alpha[i] + beta[k])  ;
      y[i,k] ~ neg_binomial(xi_i_k, omega_k_m1);             
      }
    }
  }"


overdispersed = function(survey, v_pop_total,N, warmup,iterations,chains=1){
  y0 = survey %>% dplyr::select(starts_with("KP")| HP_total_apvis )
  y <- array(dim = c(nrow(y0), ncol(y0)))
  for (i in 1:nrow(y)) {
    for (k in 1:ncol(y)) {
      y[i,k] <- y0[i,k]
    }
  }
  #Inizialization
  mu_beta = rep(NA,ncol(y))
  sigma_beta =rep(NA,ncol(y))
  mu_beta[1:(ncol(y)-1)] = log(v_pop_total)
  sigma_beta[-ncol(y)] = log(sd(y[,-ncol(y)]/survey$Reach_memory*N))
  mu_beta[ncol(y)] =  log(mean(survey$HP_total_apvis/survey$Reach_memory*N))
  sigma_beta[ncol(y)] = log(sd(survey$HP_total_apvis/survey$Reach_memory*N))
  
  data <- list(I = nrow(y), K = ncol(y), mu_beta = mu_beta, sigma_beta = sigma_beta, y = y)
  
  fit <- stan(model_code = overdispersed_model, #file='NB_norecall.stan', 
              data = data, 
              # warmup = 1000, iter = 2000,  # takes ~ 40min/chain
              warmup = warmup, iter=iterations,
              chains = chains)
  out <- extract(fit)
  return(out)
}

overdispersed2 = function(survey, v_pop_total,N, warmup,iterations,chains=1){
  y0 = survey %>% dplyr::select(starts_with("KP")| HP_total_apvis )
  y <- array(dim = c(nrow(y0), ncol(y0)))
  for (i in 1:nrow(y)) {
    for (k in 1:ncol(y)) {
      y[i,k] <- y0[i,k]
    }
  }
  #Inizialization
  mu_beta = rep(NA,ncol(y))
  sigma_beta =rep(NA,ncol(y))
  mu_alpha = rep(NA,nrow(y))
  mu_beta[1:(ncol(y)-1)] = log(v_pop_total)
  sigma_beta[-ncol(y)] = sd(log(y[,-ncol(y)]/survey$Reach_memory*N))
  mu_beta[ncol(y)] =  mean(log(survey$HP_total_apvis/survey$Reach_memory*N))
  sigma_beta[ncol(y)] = sd(log(survey$HP_total_apvis/survey$Reach_memory*N))
  for (i in 1:nrow(y)) {
    mu_alpha[i] = glm(y[i,]~mu_beta-1,family = "poisson")$coefficients
  }
  
  data <- list(I = nrow(y), K = ncol(y), mu_beta = mu_beta, mu_alpha= mu_alpha, sigma_beta = sigma_beta, y = y)
  
  fit <- stan(model_code = overdispersed_model2, #file='NB_norecall.stan', 
              data = data, 
              # warmup = 1000, iter = 2000,  # takes ~ 40min/chain
              warmup = warmup, iter=iterations,
              chains = chains)
  out <- extract(fit)
  return(out)
}


