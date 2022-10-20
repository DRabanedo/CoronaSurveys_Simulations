
########################################################################
# Functions to create populations and surveys, and the NSUM estimators #
########################################################################

#######################
# Libraries used #    #
library(igraph)       #
library(igraph)       #
library(tidyverse)    #
library(stringr)      #
library(ggplot2)      #
library(sampler)      #
library(dplyr)        #
library(truncnorm)    #
#
#######################

################################################################################

#############################
# Generation of populations #
#############################

# This function uniformly assings if the individuals belong to the hidden population

getHP <- function(n, prob_hp) {
  # Hidden population generator: returns a binary vector, the ones represent the Hidden Population
  
  # prob_hp: probability of the occurrence of the hidden population
  # n:    people in the general population
  
  vect = sample(1:n, prob_hp*n, replace = FALSE, p = NULL)
  vect_hp = rep(NA, n)
  for (i in 1:n){
    if (as.logical(sum(i %in% vect))){
      vect_hp[i] = 1
    }
    else{
      vect_hp[i] = 0
    }
    
  }
  return(vect_hp)
}


# This function generates the population with not disjoint populations

genPopulation <- function(n, prob_vect, prob_hp) {
  # Generates a data frame with the population and the belonging to the hidden population
  
  # n:         the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # prob_hp:   probability of the occurrence of the Hidden Population
  
  subpop_vect = n * prob_vect
  population_buc = data.frame(hidden_population = getHP(n,prob_hp))
  rownames(population_buc) <- c(1:n)
  
  for (k in 1:length(subpop_vect)) {
    # Index belonging to the subpopulation 
    subpop_ind = sample(1:n, subpop_vect[k], replace = FALSE)
    
    # Index transformed into a 0 & 1 vector
    subpop = rep(NA, n)
    for (i in 1:n){
      if (as.logical(sum(i %in% subpop_ind))){
        subpop[i] = 1
      }
      else{
        subpop[i] = 0
      }
      
    }
    
    #Dataframe append
    population_buc = cbind(population_buc, Subpopulation = subpop)
    names(population_buc)[dim(population_buc)[2]] = str_c("subpopulation_",k)
  }
  
  return(population_buc)
}


# This function generates a basic population with disjoint populations (Network script)

genPopulation_Disjoint_basic <- function(n, prob_vect,HP) {
  # Generates a data frame with the population and the belonging to the Hidden Population
  
  # n: the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # HP:  Hidden Population vector
  
  population_buc = data.frame("hidden_population" = HP)
  
  # Subpopulation size vector
  subpop_vect = round(n*prob_vect)
  
  # Variables for the loop
  sampling_vect = 1:n
  gen_subpop = rep(0, n)
  n_vect = 1:n
  
  for (k in 1:length(subpop_vect)) {
    # Index belonging to the subpopulation k
    subpop_ind = sample(sampling_vect, subpop_vect[k], replace = FALSE)
    
    # Index transformed into a 0 & 1 vector to represent the populations
    subpop = rep(NA, n)
    
    for (i in 1:n){
      if (as.logical(sum(i %in% subpop_ind))){
        subpop[i] = 1
        # for k in 1:n appends 1 if a population is assigned
        gen_subpop[i] = 1
      }
      else{
        subpop[i] = 0
      }
      
    }
    
    # Creates a vector with the people who does not have a population 
    sampling_vect = n_vect[gen_subpop == 0]
    
    #Dataframe append population k
    population_buc = cbind(population_buc, Subpopulation = subpop)
    names(population_buc)[dim(population_buc)[2]] = str_c("subpopulation_",k)
  }
  
  return(population_buc)
}


# This function generates the population with disjoint populations

genPopulation_Disjoint <- function(n, prob_vect,HP, M_vis, sub_mem_factor, r, r_mem, hp_t, hp_s) {
  # Generates a data frame with the population and the belonging to the Hidden Population
  
  # n: the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # HP:  Hidden Population vector
  
  population_buc = data.frame("hidden_population" = HP)
  
  # Subpopulation size vector
  subpop_vect = round(n*prob_vect)
  
  # Variables for the loop
  sampling_vect = 1:n
  gen_subpop = rep(0, n)
  n_vect = 1:n
  
  for (k in 1:length(subpop_vect)) {
    # Index belonging to the subpopulation k
    subpop_ind = sample(sampling_vect, subpop_vect[k], replace = FALSE)
    
    # Index transformed into a 0 & 1 vector to represent the populations
    subpop = rep(NA, n)
    
    for (i in 1:n){
      if (as.logical(sum(i %in% subpop_ind))){
        subpop[i] = 1
        # for k in 1:n appends 1 if a population is assigned
        gen_subpop[i] = 1
      }
      else{
        subpop[i] = 0
      }
      
    }
    
    # Creates a vector with the people who does not have a population 
    sampling_vect = n_vect[gen_subpop == 0]
    
    #Dataframe append population k
    population_buc = cbind(population_buc, Subpopulation = subpop)
    names(population_buc)[dim(population_buc)[2]] = str_c("subpopulation_",k)
  }
  
  population_buc = cbind(population_buc, reach = r)
  population_buc = cbind(population_buc, reach_memory = r_mem)
  population_buc = cbind(population_buc, hp_total = hp_t)
  population_buc = cbind(population_buc, hp_survey = hp_s)
  
  for(j in 1:(length(prob_vect)-1)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(population_buc[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(as.character(j)))) 
      vis_yij = sum(population_buc[net_sw[[i]][[1]],]["hidden_population"][as.logical(dplyr::select(population_buc[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(as.character(j)))[,1]),]) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_mem_factor*vis_pob)))
    }
    
    population_buc = cbind(population_buc,Subpoblacion_total = v_1)
    names(population_buc)[dim(population_buc)[2]] = str_c("kp_reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:(length(prob_vect)-1)) {
    for (j in ind1){
      ind2 = dplyr::select(population_buc, starts_with("subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(M_vis[ind2,j]) + 0.5, mean = sum(M_vis[ind2,j]), sd = sum(M_vis[ind2,j])*sub_mem_factor)) 
    }
    population_buc = cbind(population_buc, Subpoblacion_total = i_hp_vis)
    names(population_buc)[dim(population_buc)[2]] = str_c("kp_alters_",i)
  }
  
  
  return(population_buc)
}

# This function generates a basic population using a poisson distribution (network scripts)


getSurvey = function(n_enc, dataframe){
  # This function makes a sample of size n_enc from dataframe
  # n_enc = number of individuals interviewed
  # dataframe = general population 
  
  sur = dataframe[sample(nrow(dataframe), n_enc, replace = FALSE),]
  
  return(sur)
}


getSurvey_VF = function(n_enc, pop, vis_matrix, memory_fact){
  # This function makes a sample of size n_enc from dataframe to make an estimate 
  # of the visibility factor
  # n_enc: number of people surveyed
  # pop: Population dataframe
  # vis_matrix: Visibility matrix
  # memory_fact: Memoty factor
  
  
  enc_hp = pop[sample(nrow(pop[pop$hidden_population==1,]), n_enc, replace = FALSE),]
  
  ind_survey = as.numeric(rownames(pop[pop$hidden_population==1,]))[as.numeric(rownames(enc_hp))]
  vect_reach_hp = colSums(vis_matrix[,ind_survey])
  
  vect_reach = pop$reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double truncation + double truncation calculate
  for (i in 1:length(vect_reach_hp)) {
    
    if (vect_reach_hp[i] == vect_reach[i]){
      mem_vect_reach[i] = vect_reach_hp[i]
    } else {
      mem_vect_reach[i] = max(0,round(rtruncnorm(1, a = -0.5 + vect_reach_hp[i] , b = 0.5 + 2 * vect_reach[i] - vect_reach_hp[i], mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    }
    
    if (vect_reach_hp[i] == mem_vect_reach[i]){
      mem_vect_reach_hp[i] = mem_vect_reach[i]
    } else {
      mem_vect_reach_hp[i] = max(0,round(rtruncnorm(1, a = max( vect_reach_hp[i] - (mem_vect_reach[i]-vect_reach_hp[i]) - 0.5, -0.5), b = min(mem_vect_reach[i] + 0.5, 2 * vect_reach_hp[i] + 0.5),  mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))   
    }
    
  }
  enc_pop = pop[ind_survey,]
  enc_pop$reach_memory =  mem_vect_reach
  enc_pop = cbind(enc_pop, reach_hp = vect_reach_hp)
  enc_pop = cbind(enc_pop, reach_hp_memory = mem_vect_reach_hp)
  return(enc_pop)
  
}


getV_pop = function(n_pop, Population){
  # Number of people on each subpopulation
  # n_pop: number of subpopulations
  # Population: Population dataframe
  
  v_pop_total = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total[k] = sum(dplyr::select(Population, starts_with("subpop") & ends_with(as.character(k)) ) ) # N_k
  }
  return(v_pop_total)
}

################################################################################

########################
# Matrix for the GNSUM #
########################

matrixHP = function(grafo,Pob){
  # Adjacency matrix of the directed graph of connections with the hidden population
  
  # grafo: population graph
  # Pob: Population obtained by getData()
  
  ad = as_adj(grafo) # adjacency matrix
  for (j in 1:ncol(ad)) {
    #if (V(grafo)$label[j] == 0){
    if (Pob$hidden_population[j]==0){
      ad[,j] = 0
    }
  }
  return(ad)
}


# Visibility factor calculate for a matrix

berHP = function(x,p){
  # Binomial general function for the visibility matrix (element by element)
  
  # x: matrix 
  # p: binomial probability
  if(x!=0){
    return(x*rbinom(1,1,p))
  }
  else {
    return(0)
  }
}

################################################################################

#################################
# General Population generation #
#################################

# General population generation #

getData = function(N, prob_vect, prob_hp, dim, nei, p, vis_factor, mem_factor, sub_mem_factor){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  # N:  Population size
  # prob_vect:  Vector with the population's probabilities
  # prob_hp: Hidden Population proportion
  # dim: Integer constant, the dimension of the starting lattice.
  # nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  # p: Real constant between zero and one, the rewiring probability.
  # vis_factor: Visibility factor
  # mem_factor: Memory factor
  # sub_mem_factor: Subpopulation memory factor
  
  
  Population = genPopulation(N, prob_vect,prob_hp)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(prob_vect)
  # initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = vis_factor)
  
  
  for (i in 1:N) {
    # net_sw[[i]], list with one element, the list of the adjacent vertices to i
    
    vect_hp[i] = sum(Mhp[i,])
    vect_reach[i] = length(net_sw[[i]][[1]])
    vect_hp_vis[i] = round(rtruncnorm(1, a = max(-0.5,  2 * sum(Mhp_vis[i,]) - vect_reach[i] + 0.5 ) , b = min(2 * sum(Mhp_vis[i,]) + 0.5, vect_reach[i]-0.5), mean = sum(Mhp_vis[i,]), sd = mem_factor*sum(Mhp_vis[i,])))
    
    vect_reach_re[i] = max(1,round(rtruncnorm(1, a = vect_hp_vis[i] - 0.5, b = 2*vect_reach[i] - vect_hp_vis[i] + 0.5, mean = vect_reach[i], sd = mem_factor*vect_reach[i])))
  }
  
  Population = cbind(Population, reach = vect_reach)
  Population = cbind(Population, reach_memory = vect_reach_re)
  Population = cbind(Population, hp_total = vect_hp) 
  Population = cbind(Population, hp_survey = vect_hp_vis)
  
  for(j in 1:length(prob_vect)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(as.character(j)))) 
      vis_yij = sum(Population[net_sw[[i]][[1]],]["hidden_population"][as.logical(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(as.character(j)))[,1]),]) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_mem_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("kp_reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:length(prob_vect)) {
    for (j in ind1){
      ind2 = dplyr::select(Population, starts_with("subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_mem_factor)) 
    }
    Population = cbind(Population, Subpoblacion_total = i_hp_vis)
    names(Population)[dim(Population)[2]] = str_c("kp_alters_",i)
  }
  
  
  return(list(net_sw, Population, Mhp_vis))
}

################################################################################
# Visibility factor estimate #

VF_Estimate = function(enc_hp){
  return(sum(enc_hp$reach_hp_memory))/sum(enc_hp$reach_memory)
}

################################################################################

###################
# Basic estimator #
###################

getNh_basic_sum = function(survey,N) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  
  Nh_f =  N*sum(survey$hp_survey)/sum(survey$reach_memory)
  
  return(Nh_f)
}

getNh_basic_mean = function(survey,N) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  
  Nh_f =  N*mean(survey$hp_survey/survey$reach_memory)
  
  return(Nh_f)
}

getNh_basicvis_sum = function(survey,N,vis) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  #vis: estimation of the visibility factor
  
  Nh_f =  N*sum(survey$hp_survey)/sum(survey$reach_memory) * (1/vis)
  
  return(Nh_f)
}


getNh_basicvis_mean = function(survey,N,vis) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  #vis: estimation of the visibility factor
  
  Nh_f =  N*sum(survey$hp_survey)/sum(survey$reach_memory) * (1/vis)
  
  return(Nh_f)
}


#################
# MLE estimator #
#################

getNh_MLE = function(enc,v_pob) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  
  suma_KP = sum( dplyr::select(enc, starts_with("kp_reach_")) ) # Known Population sum
  # (\sum y_{iu})/(\frac{\sum N_k}{\sum \sum y_{ik}} )
  Nh_f = sum(enc$hp_survey)*(sum(v_pob)/suma_KP)
  
  return(Nh_f)
}


getNh_MLEvis = function(enc,v_pob,vis) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  #vis: estimation of the visibility factor  
  
  suma_KP = sum( dplyr::select(enc, starts_with("kp_reach_")) )
  Nh_f = (sum(enc$hp_survey))*(sum(v_pob)/suma_KP)*(1/vis)
  Nh_f
}



###################
# PIMLE estimator #
###################


getNh_PIMLE = function(enc,v_pob,N) {
  #NSUM Plug-in Maximum Likelihood Estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  #N: population's size
  #vis: estimation of the visibility factor
  
  #reach estimate
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum( dplyr::select(enc, starts_with("kp_reach_"))[i,] )/sum(v_pob)
  }
  Nh_f = N * mean(enc$hp_survey/d_iest)
  Nh_f
}


getNh_PIMLEvis = function(enc,v_pob,N,vis) {
  #NSUM Plug-in Maximum Likelihood Estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  #N: population's size
  #vis: estimation of the visibility factor
  
  #reach estimate
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum( dplyr::select(enc, starts_with("kp_reach_"))[i,] )/sum(v_pob)
  }
  Nh_f = N * mean(enc$hp_survey/d_iest) * (1/vis) # \frac{y_{iu}}{\hat{d_i}}
  Nh_f
}


#################
# MoS estimator #
#################


getNh_MoS = function(enc, v_pob, N){
  #NSUM Mean of Sums(MoS) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  #N: population's size
  
  # \hat{d_i} = N/L \sum_k (y_{ik}/N_k)
  
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum( dplyr::select(enc, starts_with("kp_reach_"))[i,] /v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$hp_survey/d_i_est)
  Nh_f
}

# Using dplyr in this case increases a lot the complexity of the function

getNh_MoSvis = function(enc, v_pob, N, vis){
  #NSUM Mean of Sums(MoS) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  #N: population's size
  #vis: estimation of the visibility factor 
  
  # \hat{d_i} = N/L \sum_k (y_{ik}/N_k)
  
  # reach estimate
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum(dplyr::select(enc, starts_with("kp_reach_"))[i,]/v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$hp_survey/d_i_est) * (1/vis)
  Nh_f
} 


#########
# GNSUM #
#########


getNh_GNSUM  = function(enc, enc_hp, v_pob, N){
  #General NSUM (GNSUM) (Formula from "GENERALIZING THE NETWORK SCALE-UP METHOD")
  #enc:     survey
  #enc_hp:  hidden population's survey
  #v_pob:   vector with the number of people in each Subpopulation
  #N:       population's size
  
  #Numerator estimate
  n_enc = nrow(enc)
  prob_inc = n_enc/N  #Same inclusion probability for all samples
  numerador = (1/prob_inc) * sum(enc$hp_total) #Numerator estimate
  
  #Denominator estimate
  ind1 = as.numeric(rownames(enc_hp))
  suma = sum(dplyr::select( enc_hp, starts_with("kp_alters_") ))
  denominador = N/sum(v_pob) * suma/nrow(enc_hp)      #Denominator estimate
  
  Nh = numerador/denominador
  return(Nh)
}


####################
# Direct estimator #
####################

getNh_Direct = function(survey,N){
  #Direct estimation
  #survey: survey
  #N: Population size
  
  Nh = sum(survey$hidden_population)/nrow(survey) * N
  return(Nh)
}
