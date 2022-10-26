# Teo et al. (2019) model
#########################
library(rjags)

### Models
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

model2 ='model {

for(i in 1:N)

{
  
  for(k in 1:Ku)
  
  {
  
  nu[i,k] ~ dpois(lambda*alpha[i]*exp(beta[k]*x[i,k])*Su[k])
  
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

for(k in 1:Ku)

{
  
  beta[k]~dnorm(0,(1/sigma)^2)
  
}

lambda ~ dunif(0,10)

tau ~ dunif(0,10)

sigma ~ dgamma(1,0.01)
}
'

runModel = function(survey,knowpopulation_data,NITERATION)
{
  #knownpopulation_data contains the number of individuals from each subpopulation
  popindex = colnames(survey %>% dplyr::select(starts_with("KP")| HP_total_apvis ))
  data0 = survey[,popindex]
  indexk = grep("K", colnames(data0))
  indexu = grep("H",colnames(data0))
  dataset=list(
    N=dim(data0)[1],
    Kk=length(indexk),
    nk=data0[,indexk],
    Ku=length(indexu),
    nu=as.data.frame(x=data0[,indexu],col.names=indexu),
    Sk=knowpopulation_data, 
    #Sk=as.data.frame(knowpopulation_data),
    Su=rep(NA,length(indexu)))
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(model1),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda","tau","Su"),n.iter=NITERATION,progress.bar="text",thin=10)
  dicsamples = dic.samples(jagmod,type = "pD",n.iter=20000,thin=10)
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior,dicsamples=dicsamples)
  return(results)
}
