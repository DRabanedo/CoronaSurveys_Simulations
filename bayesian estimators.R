library(asbio)
library(rjags)

getNh_Zheng = function(survey, v_pop_prob,N, iterations = 5000,burnins =1000){
  data = survey %>% dplyr::select(starts_with("kp_r") | hp_total)
  y = matrix(data = NA,nrow = nrow(data),ncol = ncol(data))
  for (i in 1:ncol(data)) {
    y[,i]=data[,i]
  }
  known = N*v_pop_prob
  
  N.mc = iterations
  N.i = nrow(y)
  N.k = ncol(y)
  
  
  
  prevalences = v_pop_prob
  pg1.ind = length(v_pop_prob)
  Pg1 = sum(prevalences[pg1.ind])
  
  ## Declare parameters
  alphas = matrix(NA, nrow = N.mc, ncol = N.i)
  betas = matrix(NA, nrow = N.mc, ncol = N.k)
  omegas = matrix(NA, nrow = N.mc, ncol = N.k)
  mu.alpha = mu.beta = sigma.sq.alpha = sigma.sq.beta = rep(NA, N.mc)
  C1 = C2 = C = NA
  
  alphas[1,] = rnorm(N.i, sd = 2)
  betas[1,] = rnorm(N.k, sd = 2)
  omegas[1,] = 20
  mu.alpha[1] = mu.beta[1] = sigma.sq.alpha[1] = sigma.sq.beta[1] = 5
  
  
  for(ind in 2:N.mc){
    ## Step 1
    for(i in 1:N.i){
      alpha.prop = alphas[ind - 1,i] + rnorm(1, 0, 0.4)
      zeta.prop = exp(alpha.prop + betas[ind - 1,]) / (omegas[ind - 1,] - 1)
      zeta.old = exp(alphas[ind - 1, i] + betas[ind - 1,]) / (omegas[ind - 1,] - 1)
      sum1 = sum(lgamma(y[i,] + zeta.prop) - lgamma(zeta.prop) - zeta.prop * log(omegas[ind - 1,])) +
        dnorm(alpha.prop, mu.alpha[ind - 1], sqrt(sigma.sq.alpha[ind - 1]), log = T)
      sum2 = sum(lgamma(y[i,] + zeta.old) - lgamma(zeta.old) - zeta.old * log(omegas[ind - 1,])) +
        dnorm(alphas[ind - 1,i], mu.alpha[ind - 1], sqrt(sigma.sq.alpha[ind - 1]), log = T)
      prob.acc = exp(sum1 - sum2)
      
      if(prob.acc > runif(1)){
        alphas[ind, i] = alpha.prop
      }else{
        alphas[ind, i] = alphas[ind - 1, i]
      }
    }
    
    ## Step 2
    for(k in 1:N.k){
      beta.prop = betas[ind - 1,k] + rnorm(1, 0, 0.2)
      zeta.prop = exp(alphas[ind, ] + beta.prop) / (omegas[ind - 1,k] - 1)
      zeta.old = exp(alphas[ind, ] + betas[ind - 1,k]) / (omegas[ind - 1,k] - 1)
      sum1 = sum(lgamma(y[,k] + zeta.prop) - lgamma(zeta.prop) - zeta.prop * log(omegas[ind - 1,k])) +
        dnorm(beta.prop, mu.beta[ind - 1], sqrt(sigma.sq.beta[ind - 1]), log = T)
      sum2 = sum(lgamma(y[,k] + zeta.old) - lgamma(zeta.old) - zeta.old * log(omegas[ind - 1,k])) +
        dnorm(betas[ind - 1,k], mu.beta[ind - 1], sqrt(sigma.sq.beta[ind - 1]), log = T)
      prob.acc = exp(sum1 - sum2)
      
      if(prob.acc > runif(1)){
        betas[ind, k] = beta.prop
      }else{
        betas[ind, k] = betas[ind - 1, k]
      }
    }
    
    ## Step 3
    mu.alpha.hat = mean(alphas[ind,])
    mu.alpha[ind] = rnorm(1, mu.alpha.hat, sqrt(sigma.sq.alpha[ind - 1] / 2))
    
    ## Step 4
    sigma.alpha.hat = mean((alphas[ind,] - mu.alpha[ind])^2)
    sigma.sq.alpha[ind] = rinvchisq(1, N.i - 1, sigma.alpha.hat)
    
    ## Step 5
    mu.beta.hat = mean(betas[ind,])
    mu.beta[ind] = rnorm(1, mu.beta.hat, sqrt(sigma.sq.beta[ind - 1] / 2))
    
    ## Step 6
    sigma.beta.hat = mean((betas[ind,] - mu.beta[ind])^2)
    sigma.sq.beta[ind] = rinvchisq(1, N.k - 1, sigma.beta.hat)
    
    
    ## Step 7
    for(k in 1:N.k){
      omega.prop = omegas[ind - 1,k] + rnorm(1, 0, 0.2)
      if(omega.prop > 1){
        zeta.prop = exp(alphas[ind, ] + betas[ind,k]) / (omega.prop - 1)
        zeta.old = exp(alphas[ind, ] + betas[ind,k]) / (omegas[ind - 1,k] - 1)
        sum1 = sum(lgamma(y[,k] + zeta.prop) - lgamma(zeta.prop) - zeta.prop * log(omega.prop) +
                     y[,k] * log((omega.prop - 1)/omega.prop))
        sum2 = sum(lgamma(y[,k] + zeta.old) - lgamma(zeta.old) - zeta.old * log(omegas[ind - 1,k]) +
                     y[,k] * log((omegas[ind - 1, k] - 1)/omegas[ind - 1,k]))
        prob.acc = exp(sum1 - sum2)
        
        if(prob.acc > runif(1)){
          omegas[ind, k] = omega.prop
        }else{
          omegas[ind, k] = omegas[ind - 1, k]
        }
      }else{
        omegas[ind,k] = omegas[ind - 1, k]
      }
    }
    
    ## Step 8
    C1 = log(sum(exp(betas[ind, pg1.ind]) / Pg1))
    C = C1
    alphas[ind,] = alphas[ind,] + C
    mu.alpha[ind] = mu.alpha[ind] + C
    betas[ind,] = betas[ind,] - C
    mu.beta[ind] = mu.beta[ind] - C
    
    if(ind %% 1 == 0){
      print(ind)
    }
  }
  
  ## Burn-in and thin
  alphas = alphas[-c(1:burnins),]
  betas = betas[-c(1:burnins),]
  alphas = alphas[seq(1, nrow(alphas), by = 10),]
  betas = betas[seq(1, nrow(betas), by = 10),]
  
  return(mean(exp(betas[,length(prevalences)+1]) * N))
}

### Teo et al. model

model1 = 'model {
for(i in 1:N)
{
  
  for(k in 1:Ku)
  
  {
  
  nu[i,k] ~ dpois(lambda*alpha[i]*Su[k])
  
  }
  
  for(k in 1:Kk)
  
  {
  
  nk[i,k] ~ dpois(lambda*alpha[i]*Sk[k])
  
  }
  
  alpha[i]~dlnorm(0,tau)
  
}
for(k in 1:Ku)
{
  
  Su[k]~dunif(0,2500000)
  
}
for(k in 1:Kk)
{
  
  Sk[k]~dunif(0,2500000)
  
}
lambda ~ dunif(0,10)
tau ~ dunif(0,10)
}
'

runModel = function(indexu,indexk,NITERATION,y,v_pop_prob)
{
  dataset=list(
    N=dim(y)[1],
    Kk=length(indexk),
    nk=y[,indexk],
    Ku=length(indexu),
    nu=as.matrix(y[,indexu], ncol = length(indexu)),
    #Sk=Curitiba$known[indexk],
    Sk = N* v_pop_prob, 
    Su=rep(NA,length(indexu)))
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(model1),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda","tau","Su"),n.iter=NITERATION,progress.bar="text",thin=10)
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior)
  return(results)
}

getNh_TEO = function(survey,v_pop_prob,N,iter){
  data = survey %>% dplyr::select(starts_with("kp_r") | hp_total)
  y = matrix(data = NA,nrow = nrow(data),ncol = ncol(data))
  for (i in 1:ncol(data)) {
    y[,i]=data[,i]
  }
  #known = N*v_pop_prob
  indexu = length(v_pop_prob)+1
  indexk = 1:length(v_pop_prob)
  teo.basic.res = runModel(indexu, indexk, iter,y,v_pop_prob)
  teo.basic.post = teo.basic.res$posterior
  teo.basic.su = c(teo.basic.post[[1]][,1], teo.basic.post[[2]][,1])
  
  return(mean(teo.basic.su))
}