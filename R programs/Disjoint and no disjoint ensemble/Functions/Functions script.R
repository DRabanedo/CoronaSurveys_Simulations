
########################################################################
# Functions to create populations and surveys, and the NSUM estimators #
########################################################################

#######################
                      #
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

  sample(c(0,1), n, replace = TRUE, p = c(1-prob_hp,prob_hp))
}


# This function generates the population with not disjoint populations

genPopulation <- function(n, prob_vect, prob_hp) {
  # Generates a data frame with the population and the belonging to the hidden population
  
  # n:         the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # prob_hp:   probability of the occurrence of the Hidden Population
  
  pop_vect = 1:length(prob_vect)
  population_buc = data.frame(Hidden_Population = getHP(n,prob_hp))
  rownames(population_buc) <- c(1:n)
  
  for (i in 1:length(pop_vect)) {
    population_buc = cbind(population_buc, Subpopulation = sample(c(0,1), n, replace = TRUE, p = c(1-prob_vect[i],prob_vect[i])))
    names(population_buc)[dim(population_buc)[2]] = str_c("Subpopulation_",i)
  }

  return(population_buc)
}


# This function generates the population with disjoint populations

genPopulation_Disjoint <- function(n, prob_vect,HP) {
  # Generates a data frame with the population and the belonging to the Hidden Population
  
  # n: the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # HP:  Hidden Population vector
  
  population_buc = data.frame("Hidden_Population" = HP)
  
  #Population 0 introduction
  prob_vect = c(1-sum(prob_vect), prob_vect)
  Subpop_vector = sample(0:(length(prob_vect)-1), n, replace = TRUE, p = prob_vect)
  
  for (i in 1:(length(prob_vect)-1)) {
    population_buc = cbind(population_buc, Subpopulation_buc = as.integer(Subpop_vector == i))
    names(population_buc)[dim(population_buc)[2]] = str_c("Subpopulation_",i)
  }
  
  return(population_buc)
}


# This function generates the population using a poisson distribution

genPopulation_poisson <- function(n, prob_vect, prob_hp) {
  # Generates a data frame with the population and the belonging to the hidden population
  
  # n: the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # HP:  Hidden Population vector
  
  population_buc = data.frame(Hidden_Population = getHP(n, prob_hp))
  
  pop_vect_num = rep(NA, length(prob_vect))
  for (j in 1:length(prob_vect)) {
    pop_vect_num[j] = rpois(1,prob_vect[j]*n)
  }
  
  for (i in 1:(length(prob_vect))) {
    sam_pop = sample(1:n, size = pop_vect_num[i])
    Subpop_vector = rep(0,n)
    for (k in sam_pop){
      Subpop_vector[k] = 1
    }
    population_buc = cbind(population_buc, Subpopulation_buc = Subpop_vector)
    names(population_buc)[dim(population_buc)[2]] = str_c("Subpopulation_",i)
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
  # This function makes a sample of size n_enc from dataframe to make an estimate 
  # of the visibility factor
  enc_hp = dataframe[sample(nrow(pop[pop$Hidden_Population==1,]), n_enc, replace = FALSE),]
  
  ind_survey = as.numeric(rownames(enc_hp))
  vect_reach_hp = colSums(vis_matrix[,ind_survey])
  
  vect_reach = Population$Reach[ind_survey]
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
      mem_vect_reach_hp[i] = max(0,round(rtruncnorm(1, a = vect_reach_hp[i] - (mem_vect_reach[i]-vect_reach_hp[i]), b = mem_vect_reach[i], mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))   
    }
    
  }
  
  pop$Reach_Memory =  mem_vect_reach
  pop = cbind(pop, Reach_HP = vect_reach_hp)
  pop = cbind(pop, Reach_HP_Memory = mem_vect_reach_hp)
  
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
    if (Pob$Hidden_Population[j]==0){
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

getData = function(N, prob_vector, prob_hp, dim, nei, p, vis_factor, mem_factor, sub_mem_factor){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  # N:  Population size
  # prob_vector:  Vector with the population's probabilities
  # prob_hp: Hidden Population proportion
  # dim: Integer constant, the dimension of the starting lattice.
  # nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  # p: Real constant between zero and one, the rewiring probability.
  # vis_factor: Visibility factor
  # mem_factor: numeric value, reach divided by the standard deviation of the normal we use to correct the reach
  # sub_mem_factor:  the Subpopulations visibility divided by the standard deviation of the normal distributions we use to correct the Subpopulations visibility
  #                  it is applied to each Subpopulation
  
  
  Population = genPopulation(N, prob_vector,prob_hp)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(prob_vector)
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
    vect_hp_vis[i] = max(0,round(rtruncnorm(1, a = -0.5, b = 2 * sum(Mhp_vis[i,]) + 0.5, mean = sum(Mhp_vis[i,]), sd = mem_factor*sum(Mhp_vis[i,]))))
    
    vect_reach_re[i] = max(1,round(rtruncnorm(1, a = vect_hp_vis[i] - 0.5, b = Inf, mean = vect_reach[i], sd = mem_factor*vect_reach[i])))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_Memory = vect_reach_re)
  Population = cbind(Population, HP_Total = vect_hp) 
  Population = cbind(Population, HP_Survey = vect_hp_vis)
  
  for(j in 1:length(prob_vector)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("Subpop") & ends_with(as.character(j)))) 
      vis_yij = sum(Population[net_sw[[i]][[1]],]["Hidden_Population"][as.logical(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("Subpop") & ends_with(as.character(j)))[,1]),]) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_mem_factor*vis_pob)))
    }

    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_Reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:length(prob_vector)) {
    for (j in ind1){
      ind2 = dplyr::select(Population, starts_with("Subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_mem_factor)) 
    }
    Population = cbind(Population, Subpoblacion_total = i_hp_vis)
    names(Population)[dim(Population)[2]] = str_c("KP_Alters_",i)
  }
  
  
  return(list(net_sw, Population, Mhp_vis))
}




getData_poisson = function(N, prob_vector, prob_hp, dim, nei, p, vis_factor, mem_factor, sub_mem_factor){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  #N: population size
  #dis_populations: vector with the populations
  #prob_vector: vector with the population's probabilities
  #prob_hp: Hidden Population proportion
  #dim: Integer constant, the dimension of the starting lattice.
  #nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p: Real constant between zero and one, the rewiring probability.
  #vis_factor: the visibility factor
  #memory factor: numeric value, reach divided by the standard deviation of the normal we use to correct the reach
  # sub_mem_factor: the Subpopulations visibility divided by the standard deviation of the normal distributions we use to correct the Subpopulations visibility
  #                    it is applied to each Subpopulation
  # sub_mem_factor: the Subpopulations' visibility divided by the standard deviation of the normals we use to correct the Subpopulations' visibility
  
  Population = genPopulation_poisson(N, prob_vector, prob_hp)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(dis_populations)
  # initializes the vectors
  vect_hp = rep(NA,N)        # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)    # vect_hp applying visibility
  vect_reach = rep(NA,N)     # the degrees of each individual
  vect_reach_re = rep(NA,N)  # reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = vis_factor)
  
  
  for (i in 1:N) {
    # net_sw[[i]], list with one element, the list of the adjacent vertices to i
    
    vect_hp[i] = sum(Mhp[i,])
    vect_reach[i] = length(net_sw[[i]][[1]])
    vect_hp_vis[i] = max(0,round(rnorm(1, mean = sum(Mhp_vis[i,]), sd = mem_factor*sum(Mhp_vis[i,]))))
    vect_reach_re[i] = max(1,round(rtruncnorm(1, a = vect_hp_vis[i] - 0.5 , b = 2*vect_reach[i] - vect_hp_vis[i] + 0.5, mean = vect_reach[i], sd = mem_factor*vect_reach[i])))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_Memory = vect_reach_re)
  Population = cbind(Population, HP_Total = vect_hp) 
  Population = cbind(Population, HP_Survey = vect_hp_vis)
  
  for(j in 1:length(prob_vector)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("Subpop") & ends_with(as.character(j)))) 
      vis_yij = sum(Population[net_sw[[i]][[1]],]["Hidden_Population"][as.logical(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("Subpop") & ends_with(as.character(j)))[,1]),]) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_mem_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_Reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:length(prob_vector)) {
    for (j in ind1){
      ind2 = dplyr::select(Population, starts_with("Subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_mem_factor)) 
    }
    Population = cbind(Population, Subpoblacion_total = i_hp_vis)
    names(Population)[dim(Population)[2]] = str_c("KP_Alters_",i)
  }
  
  return(list(net_sw, Population, Mhp_vis))
}

################################################################################

# Visibility factor estimate #

vf_dt_dt = function(enc_hp){
  return(sum(enc_hp$Reach_HP_Memory))/sum(enc_hp$Reach_Memory)
}


################################################################################

###################
# Basic estimator #
###################

getNh_basic_sum = function(survey,N) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  
  Nh_f =  N*sum(survey$HP_Survey)/sum(survey$Reach_Memory)
  
  return(Nh_f)
}

getNh_basic_mean = function(survey,N) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  
  Nh_f =  N*mean(survey$HP_Survey/survey$Reach_Memory)
  
  return(Nh_f)
}

getNh_basicvis_sum = function(survey,N,vis) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  #vis: estimation of the visibility factor
  
  Nh_f =  N*sum(survey$HP_Survey)/sum(survey$Reach_Memory) * (1/vis)
  
  return(Nh_f)
}


getNh_basicvis_mean = function(survey,N,vis) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  #vis: estimation of the visibility factor
  
  Nh_f =  N*sum(survey$HP_Survey)/sum(survey$Reach_Memory) * (1/vis)
  
  return(Nh_f)
}


#################
# MLE estimator #
#################

getNh_MLE = function(enc,v_pob) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  
  suma_KP = sum( dplyr::select(enc, starts_with("KP_Reach_")) ) # Known Population sum
  # (\sum y_{iu})/(\frac{\sum N_k}{\sum \sum y_{ik}} )
  Nh_f = sum(enc$HP_Survey)*(sum(v_pob)/suma_KP)
  
  return(Nh_f)
}


getNh_MLEvis = function(enc,v_pob,vis) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each Subpopulation
  #vis: estimation of the visibility factor  
  
  suma_KP = sum( dplyr::select(enc, starts_with("KP_Reach_")) )
  Nh_f = (sum(enc$HP_Survey))*(sum(v_pob)/suma_KP)*(1/vis)
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
    d_iest[i] = N * sum( dplyr::select(enc, starts_with("KP_Reach_"))[i,] )/sum(v_pob)
  }
  Nh_f = N * mean(enc$HP_Survey/d_iest)
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
    d_iest[i] = N * sum( dplyr::select(enc, starts_with("KP_Reach_"))[i,] )/sum(v_pob)
  }
  Nh_f = N * mean(enc$HP_Survey/d_iest) * (1/vis) # \frac{y_{iu}}{\hat{d_i}}
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
    d_i_est[i] = (sum( dplyr::select(enc, starts_with("KP_Reach_"))[i,] /v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$HP_Survey/d_i_est)
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
    d_i_est[i] = (sum(dplyr::select(enc, starts_with("KP_Reach_"))[i,]/v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$HP_Survey/d_i_est) * (1/vis)
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
  numerador = (1/prob_inc) * sum(enc$HP_Survey) #Numerator estimate
  
  #Denominator estimate
  ind1 = as.numeric(rownames(enc_hp))
  suma = 0
  for (i in 1:length(v_pob)){
    suma = sum(dplyr::select(enc_hp, starts_with("KP_Alters_") )[dplyr::select(enc_hp, starts_with("Subpop") & ends_with(as.character(i)))[,1] == 1,] ) + suma
  }
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
  
  Nh = sum(survey$Hidden_Population)/nrow(survey) * N
  return(Nh)
}
