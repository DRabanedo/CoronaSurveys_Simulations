
########################################################################
# Functions to create populations and surveys, and the NSUM estimators #
########################################################################

#######################
# Libraries used #    #
library(igraph)       #
library(tidyverse)    #
library(stringr)      #
library(ggplot2)      #
library(sampler)      #
library(dplyr)        #
library(truncnorm)    #
library(rjags)        #
library(rstan)        #
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
# This function assings using a SIR structure the individuals who belong to the hidden population
gen_SIRpop = function(N, net ,beta, gamma, chosen_nodes, infected_people){
  # SIR method for determination of the hidden population
  
  # beta: infection rate
  # gamma: removal rate
  # infected_people: number of infected people
  # chosen_nodes: nodes from which we start
  
  # Loop variables
  infected_nodes = c()
  infected_recovered = c()
  
  infected_nodes = c(infected_nodes,chosen_nodes)
  
  while (length(infected_nodes) < infected_people) {
    
    #Infection stage
    for (node in infected_nodes) {
      for (neighbor in net[[node]][[1]]){
        if (runif(1) < beta & !(neighbor %in% infected_nodes) & !(neighbor %in% infected_recovered )) {
          infected_nodes = c(infected_nodes,neighbor)
          #print(infected_nodes)
        }
      }
    }
    
    # Removal stage
    infected_survivors = c()
    
    for (node in infected_nodes) {
      if (runif(1) < gamma){
        infected_recovered = c(infected_recovered,node)
      }
      else{
        infected_survivors = c(infected_survivors,node)
      }
    }
    infected_nodes = infected_survivors
  }
  
  final_infected_nodes = infected_nodes[1:infected_people] # CUIDADO (cambiar esto)
  
  hp_vector = rep(NA, N)
  
  for (i in 1:N){
    if (as.logical(sum(i %in% final_infected_nodes))){
      hp_vector[i] = 1
    }
    else{
      hp_vector[i] = 0
    }
  }
  return(hp_vector)
}



# This function generates the population with not disjoint populations

genPopulation <- function(n, prob_vect, prob_hp) {
  # Generates a data frame with the population and the belonging to the hidden population
  
  # n:         the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # prob_hp:   probability of the occurrence of the Hidden Population
  
  subpop_vect = round(n * prob_vect)
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
  
  for(j in 1:(length(prob_vect))){
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
  for (i in 1:(length(prob_vect))) {
    for (j in ind1){
      ind2 = dplyr::select(population_buc, starts_with("subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(M_vis[ind2,j]) + 0.5, mean = sum(M_vis[ind2,j]), sd = sum(M_vis[ind2,j])*sub_mem_factor)) 
    }
    population_buc = cbind(population_buc, Subpoblacion_total = i_hp_vis)
    names(population_buc)[dim(population_buc)[2]] = str_c("kp_alters_",i)
  }
  
  
  return(population_buc)
}


getSurvey = function(n_enc, dataframe){
  # This function makes a sample of size n_enc from dataframe
  # n_enc = number of individuals interviewed
  # dataframe = general population 
  
  sur = dataframe[sample(nrow(dataframe), n_enc, replace = FALSE),]
  
  return(sur)
}


getSurvey_VF = function(n_enc, pop, vis_matrix, memory_fact){
  # This function makes a ordered sample of size n_enc from dataframe to make an estimate 
  # of the visibility factor
  
  # n_enc: number of people surveyed
  # pop: Population dataframe
  # vis_matrix: Visibility matrix
  # memory_fact: Memoty factor
  
  # Sample from the hidden population
  enc_hp = pop[sample(nrow(pop[pop$hidden_population==1,]), n_enc, replace = FALSE),]
  
  # Survey index
  ind_survey = as.numeric(rownames(pop[pop$hidden_population==1,]))[as.numeric(rownames(enc_hp))]
  
  # Known variables 
  vect_reach_hp = colSums(vis_matrix[,ind_survey])
  vect_reach = pop$reach[ind_survey]
  
  # New variables
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
  # Output dataframe construction
  enc_pop = pop[ind_survey,]
  
  # New reach_memory variable
  enc_pop$reach_memory =  mem_vect_reach
  
  # New variables
  enc_pop = cbind(enc_pop, reach_hp = vect_reach_hp)
  enc_pop = cbind(enc_pop, reach_hp_memory = mem_vect_reach_hp)
  
  # Ordering the dataframe by index (future needs)
  enc_pop = enc_pop[order(as.numeric(row.names(enc_pop))), ]
  
  
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

to_matrix = function(x){
  # Converts to matrix the ad. matix
  if(x!=0){
    return(1)
  }
  else {
    return(0)
  }
}

to_matrix_SIR = function(x) {
  ifelse(x %in% c(0,2,3), 0, 1)
}

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
  ad = apply(ad, c(1,2), to_matrix)
  return(ad)
}


# Visibility factor calculate for an adjacent matrix (apply method)

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

getData = function(N, prob_vect, prob_hp, dim, nei, p, vis_factor, mem_factor, sub_mem_factor, beta = 0.3, gamma = 0.1, n_chosen_nodes = 2){
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
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  # SIR generation of the hidden population 
  infected_people   = round(prob_hp*N)
  chosen_nodes      = round(seq(1,N,N/n_chosen_nodes))
  
  hidden_population = gen_SIRpop(N, net_sw ,beta, gamma, chosen_nodes, infected_people)
  
  # First we generate the population, then we remplace the apropiate column
  Population = genPopulation(N, prob_vect,prob_hp) # hidden population without SIR
  Population$hidden_population = hidden_population # hidden population using SIR
  
  n_populations = length(prob_vect)
  # initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = vis_factor)
  
  # Mhp_hp calculation for making two different Bernouillis
  
  # Mhp_vis_comp = 2*Mhp_vis - 1
  # Mhp_hp = 1*(t(Mhp_vis) ==  Mhp_vis_comp)
  # Mhp_hp_vis = apply(Mhp_hp, c(1,2), berHP, p = (0.5 +  0.5*vis_factor))
  
  # Matrix for knowing with bernoilli should be applied
  # Mhp_sir = Mhp_vis + 2*Mhp_hp + 2*Mhp_hp_vis
  
  # Final visibility matrix
  # Mhp_vis = apply(Mhp_sir, c(1,2), to_matrix_SIR)
  
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


##################################
# Teo et al. (2019) bayesian model
##################################

modelTeo = 'model {

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

getNh_Teo = function(survey,knowpopulation_data,NITERATION)
{
  #knownpopulation_data contains the number of individuals from each subpopulation
  popindex = colnames(survey %>% dplyr::select(starts_with("kp_reach_")| hp_survey ))
  data0 = survey[,popindex]
  indexk = grep("k", colnames(data0))
  indexu = grep("h",colnames(data0))
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
  jagmod=jags.model(textConnection(modelTeo),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda","tau","Su"),n.iter=NITERATION,progress.bar="text",thin=1)
  dicsamples = dic.samples(jagmod,type = "pD",n.iter=20000,thin=1)
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior,dicsamples=dicsamples)
  # uses only the first chain to obtain the hidden population estimation
  Nh = mean(as.matrix(results$posterior[[1]])[,1])
  return(list(results,Nh))
}

#########################################
# Zheng et al. (2005) overdispersed model
#########################################

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

getNh_overdispersed = function(survey, v_pop_total,N, warmup,iterations,chains=1){
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
  beta_post <- out $ beta
  beta_hat <- apply(beta_post, 2, mean)
  Nh = exp(beta_hat[length(beta_hat)])
  return(list(out,Nh))
}

################################################################################
# Functions for the graphs

data_analysis = function(Nh_df, Nh_ref_df){
  # Estimation dataframe analysis
  
  df_analysis = data.frame( abserror = rowMeans(as.matrix(abs(Nh_df-Nh_ref_df))),
                            mse      = rowMeans(as.matrix((Nh_df-Nh_ref_df)^2)),
                            bias     = rowMeans(as.matrix(Nh_df)),
                            sd       = rowSds(as.matrix(Nh_df)),
                            median   = rowMedians(as.matrix(Nh_df)) )
  
  return(df_analysis)
}




