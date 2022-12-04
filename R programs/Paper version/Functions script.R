
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
library(gridExtra)    #
library(cowplot)      #
library(rjags)        #
#######################

################################################################################

#############################
# Generation of populations #
#############################

####################################
## Hidden population distribution ##

# This function uniformly assings if the individuals belong to the hidden population
gen_HP <- function(n, prob_hp) {
  # Hidden population generator: returns a dataframe with a single column,  
  # the ones represent belonging to the Hidden Population
  
  # prob_hp: probability of the occurrence of the hidden population (bernoulli parameter)
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
  population_buc = data.frame(hidden_population = vect_hp)
  return(population_buc)
}


# This function assings using a SIR structure the individuals who belong to the hidden population
gen_SIRpop = function(n, net ,beta, gamma, chosen_nodes, n_iter){
  # SIR method to determine the distribution of the hidden population
  
  # n: number of individuals in the population
  # beta: infection rate
  # gamma: removal rate
  # n_iter: number of infected people
  # chosen_nodes: nodes from which we start the infection
  
  # Loop variables
  infected_nodes = c()        # Infected nodes
  infected_recovered = c()    # Nodes recovered from the infection
  
  infected_nodes = c(infected_nodes,chosen_nodes)
  for (w in 1:n_iter) {
    
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
  # Final nodes
  final_infected_nodes = infected_survivors
  
  # 0 & 1 vector conversion
  hp_vector = rep(NA, n)
  
  for (i in 1:n){
    if (as.logical(sum(i %in% final_infected_nodes))){
      hp_vector[i] = 1
    }
    else{
      hp_vector[i] = 0
    }
  }
  
  # Dataframe output
  population_buc = data.frame(hidden_population = hp_vector)
  
  return(population_buc) 
}

######################################
## Reach & Hidden population number ##

# Gets the Reach value of each individual in the network
gen_Reach = function(net){
  # Reach calculation of each individual from the network
  # net: network to be analysed
  
  # Net size 
  n = net$size
  
  # Loop variables
  vect_reach = rep(NA,n)    # the degrees of each individual
  
  for (i in 1:n) {
    # net[[i]], list with one element, the list of the adjacent vertices to i
    vect_reach[i] = length(net[[i]][[1]])
  }
  
  # Dataframe output
  population_buc = data.frame(reach = vect_reach)
  
  return(population_buc)
}

# Calculates the people each individual knows of the hidden population
gen_Reach_hp = function(M_hp){
  # Number of neighbors from the hidden population that each individual knows
  # M_hp: Visibility matrix (before applying the visibility factor)
  
  # Population size
  n = nrow(M_hp)
  
  # Loop variables
  vect_hp = rep(NA,n)      # number of hidden population individuals known by each person
  
  for (i in 1:n) {
    vect_hp[i] = sum(M_hp[i,])
  }
  
  # Dataframe output
  population_buc = data.frame(hp_total = vect_hp)
  
  return(population_buc)
  
}

gen_Reach_hp_memory = function(population_buc, M_vis, mem_factor){
  # Number of neighbors from the hidden population that each individual knows after
  # applying the visibility factor (transmission error) and the memory factor(recall error)
  
  # population_buc: Population dataframe with at least the "hp_total" variable
  # M_vis: Visibility matrix (after applying the visibility factor)
  # mem_factor: memory factor (recall error mesure)
  
  # Population size
  n = nrow(M_vis)
  
  # initializes the vectors
  vect_reach = population_buc$reach      # number of hidden population individuals known for each person
  
  vect_hp_vis = rep(NA,n)                # vect_hp applying visibility

  for (i in 1:n) {
    vect_hp_vis[i] = round(rtruncnorm(1, a = max(-0.5,  2 * sum(M_vis[i,]) - vect_reach[i] + 0.5 ) , b = min(2 * sum(M_vis[i,]) + 0.5, vect_reach[i]-0.5), mean = sum(M_vis[i,]), sd = mem_factor*sum(M_vis[i,])))
  }
  
  # Dataframe output
  population_buc = data.frame(hp_survey = vect_hp_vis)
  
  return(population_buc)
  
}


# Calculates the reach applied the memory factor
gen_Reach_memory = function(population_buc, mem_factor){
  # Number of neighbors that each individual knows after applying the memory factor (recall error)
  
  # population_buc: Population dataframe with at least the "hp_survey" & the "reach" variable
  # mem_factor: memory factor (recall error mesure)
  
  # Population size
  n = nrow(population_buc)
  
  # Loop variables
  vect_hp_vis   = population_buc$hp_survey   # vect_hp applying visibility
  vect_reach    = population_buc$reach       # the degrees of each individual
  vect_reach_re = rep(NA,n)                # reach vector applying memory error
  
  for (i in 1:n) {
    vect_reach_re[i] = max(1,round(rtruncnorm(1, a = vect_hp_vis[i] - 0.5, b = 2*vect_reach[i] - vect_hp_vis[i] + 0.5, mean = vect_reach[i], sd = mem_factor*vect_reach[i])))
  }
  
  # Dataframe output
  population_buc = data.frame(reach_memory = vect_reach_re)
  return(population_buc)
}

###############################
## Subpopulations generation ##

# This function generates the subpopulations uniformly (not disjoint)
gen_Subpopulation = function(n, prob_vect){ 
  # Function that generates a dataframe with the subpopulations uniformly distributed 
  # in a not disjoint way, so each individual can belong to several subpopulations.
  
  # n:         number of people in the subpopulation
  # prob_vect: vector with the subpopulations size 
  
  # Loop variables
  subpop_vect = round(n * prob_vect)       # Number of individuals on each subpopulations
  population_buc = data.frame(data = 1:n)  # Dataframe initialisation
  rownames(population_buc) <- c(1:n)       # Dataframe rownames
  
  for (k in 1:length(subpop_vect)) {
    # Index to determine the belonging to each subpopulation (SAMPLE)
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
  
  return(dplyr::select(population_buc, -starts_with("data") ))
}



# This function generates the subpopulations uniformly (disjoint)
gen_Subpopulation_disjoint = function(n, prob_vect){
  # Function that generates a dataframe with the subpopulations uniformly distributed 
  # in a disjoint way, so each individual belongs to a unique subpopulations.
  
  # n:         number of people in the subpopulation
  # prob_vect: vector with the subpopulations size 
  
  population_buc = data.frame(data = 1:n)
  
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
  
  return(dplyr::select(population_buc, -starts_with("data") ))
  
  
}  




# Subpopulation memory factor function
gen_Subpopulation_memoryfactor = function(population_buc, M_vis, sub_mem_factor, net){
  # Function that appends to a dataframe the memory factor (recall error) applied
  # to the subpopulation present
  
  # population_buc: dataframe with at least the belonging to the subpopulations 
  # M_vis: Visibility matrix (after applying the visibility factor)
  # sub_mem_factor: memory factor for the subpopulation
  
  # Loop variables
  n     = nrow(population_buc)
  n_pop = ncol(dplyr::select(population_buc, starts_with("subpop")))
  ind1  = 1:n
  
  for(j in 1:n_pop){
    v_1  = rep(NA,n)
    # Index of the people that belongs to the subpopulation j
    ind2 = dplyr::select(population_buc, starts_with("subpop") & ends_with( str_c("_", as.character(j) )))[,1] != 0
    for(i in ind1) {
      # Number of people that belongs to the hidden population in the subpopulations
      vis_yij = sum(M_vis[i,ind2]) 
      vis_pob = sum(dplyr::select(population_buc[net[[i]][[1]],],starts_with("subpop") & ends_with( str_c("_", as.character(j) )))) 
      
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_mem_factor*vis_pob)))
    }
    
    population_buc = cbind(population_buc,Subpoblacion_total = v_1)
    names(population_buc)[dim(population_buc)[2]] = str_c("kp_reach_",j)
  }
  
  return(population_buc)
}




# Subpopulation alters memory factor
gen_Subpopulation_alters_memoryfactor = function(population_buc, M_vis, sub_mem_factor){
  # Function that appends to a dataframe the memory factor (recall error) applied
  # to the subpopulation alters
  
  # population_buc: dataframe with at least the belonging to the subpopulations 
  # M_vis: Visibility matrix (after applying the visibility factor)
  # sub_mem_factor: memory factor for the subpopulation
  
  # Loop variables
  n     = nrow(population_buc)
  n_pop = ncol(dplyr::select(population_buc, starts_with("subpop")))
  ind1  = 1:n

  for (i in 1:n_pop) {
    i_hp_vis = rep(NA,n)
    ind2 = dplyr::select(population_buc, starts_with("subpop") & ends_with( str_c("_", as.character(i) )))[,1] != 0
    for (j in ind1){
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(M_vis[ind2,j]) + 0.5, mean = sum(M_vis[ind2,j]), sd = sum(M_vis[ind2,j])*sub_mem_factor)) 
    }
    population_buc = cbind(population_buc, Subpoblacion_total = i_hp_vis)
    names(population_buc)[dim(population_buc)[2]] = str_c("kp_alters_",i)
  }  
  return(population_buc)
  
}

# Vector with the subpopulation number
getV_pop = function(n_pop, population_buc){
  # Number of people on each subpopulation
  # n_pop: number of subpopulations
  # Population: Population dataframe
  
  v_pop_total = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total[k] = sum(dplyr::select(population_buc, starts_with("subpop") & ends_with( str_c("_", as.character(k) )) ) ) # N_k
  }
  return(v_pop_total)
}

#############
## Surveys ##

# Uniform survey
gen_Survey = function(n_enc, dataframe){
  # This function makes a sample of size n_enc from dataframe
  
  # n_enc = number of individuals interviewed
  # dataframe = general population 
  
  sur = sample(nrow(dataframe), n_enc, replace = FALSE)
  
  # Ordering the dataframe by index
  sur = sur[order(sur)]
  
  return(sur)
}

# Survey for the visibility factor estimate
gen_Survey_VF = function(n_enc, pop, vis_matrix, memory_fact){
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

matrixHP_visibility = function(M_hp, vis_factor){
  M_vis = apply(M_hp,c(1,2), berHP,p = vis_factor)
  return(M_vis)
}

################################################################################

#################################
# General Population generation #
#################################

# General population generation #

# Uniform hidden populatio distribution
gen_Data_uniform = function(n, prob_vect, prob_hp, vis_factor, mem_factor, sub_mem_factor, beta = 0.115, gamma = 0.115/1.5, chosen_nodes = 1, n_iter = 5, net){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  # N:  Population size
  # prob_vect:  Vector with the population's probabilities
  # prob_hp: Hidden Population proportion
  # vis_factor: Visibility factor
  # mem_factor: Memory factor
  # sub_mem_factor: Subpopulation memory factor
  
  # Subpopulation dataframe
  subpop_df = gen_Subpopulation(n, prob_vect)
  
  # Hidden population distribution dataframe
  hp_df = gen_HP(n, prob_hp)
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  M_hp     =  matrixHP(net, hp_df)
  M_vis =  matrixHP_visibility(M_hp, vis_factor)
  
  # Dataframe of the population generation
  population_buc  = hp_df #Hidden population
  population_buc  = cbind(population_buc, subpop_df) #Subpopulations
  population_buc  = cbind(population_buc, gen_Reach(net)) #Reach variable
  population_buc  = cbind(population_buc, gen_Reach_hp(M_hp)) # HP reach variable
  population_buc  = cbind(population_buc, gen_Reach_hp_memory(population_buc, M_vis, mem_factor)) # HP reach recall error variable
  population_buc  = cbind(population_buc, gen_Reach_memory(population_buc, mem_factor)) #Reach recall error variable
  population_buc  = gen_Subpopulation_memoryfactor(population_buc, M_vis, sub_mem_factor, net)
  population_buc  = gen_Subpopulation_alters_memoryfactor(population_buc, M_vis, sub_mem_factor)
  
  # Returns the netwpork, the dataframe with the population data and the visibility matrix
  return(list(net, population_buc, M_vis))
}

# Epidemic data distribution
gen_Data_SIR = function(n, prob_vect, prob_hp, vis_factor, mem_factor, sub_mem_factor, beta = 0.115, gamma = 0.115/1.5, chosen_nodes = 1, n_iter = 5, net){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  # N:  Population size
  # prob_vect:  Vector with the population's probabilities
  # prob_hp: Hidden Population proportion
  # vis_factor: Visibility factor
  # mem_factor: Memory factor
  # sub_mem_factor: Subpopulation memory factor
  
  # Subpopulation dataframe
  subpop_df = gen_Subpopulation(n, prob_vect)
  
  # Hidden population distribution dataframe
  hp_df = gen_SIRpop(n, net ,beta, gamma, chosen_nodes, n_iter)
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  M_hp     =  matrixHP(net, hp_df)
  M_vis =  matrixHP_visibility(M_hp, vis_factor)
  
  # Dataframe of the population generation
  population_buc  = hp_df #Hidden population
  population_buc  = cbind(population_buc, subpop_df) #Subpopulations
  population_buc  = cbind(population_buc, gen_Reach(net)) #Reach variable
  population_buc  = cbind(population_buc, gen_Reach_hp(M_hp)) # HP reach variable
  population_buc  = cbind(population_buc, gen_Reach_hp_memory(population_buc, M_vis, mem_factor)) # HP reach recall error variable
  population_buc  = cbind(population_buc, gen_Reach_memory(population_buc, mem_factor)) #Reach recall error variable
  population_buc  = gen_Subpopulation_memoryfactor(population_buc, M_vis, sub_mem_factor, net)
  population_buc  = gen_Subpopulation_alters_memoryfactor(population_buc, M_vis, sub_mem_factor)
  
  # Returns the netwpork, the dataframe with the population data and the visibility matrix
  return(list(net, population_buc, M_vis))
}

# This function generates the population with disjoint populations
gen_Population_disjoint <- function(n, net, prob_vect, HP, M_vis, sub_mem_factor, r, r_mem, hp_t, hp_s) {
  # Generates the entire data for the population
  
  # n: the number of individuals
  # prob_vect: vector with the Subpopulations probabilities
  # HP:  Hidden Population vector
  
  population_buc  = data.frame("hidden_population" = HP)
  population_buc  = cbind(population_buc, gen_Subpopulation_disjoint(n, prob_vect))
  population_buc  = cbind(population_buc, reach = r)
  population_buc  = cbind(population_buc, reach_memory = r_mem)
  population_buc  = cbind(population_buc, hp_total = hp_t)
  population_buc  = cbind(population_buc, hp_survey = hp_s)
  population_buc  = gen_Subpopulation_memoryfactor(population_buc, M_vis, sub_mem_factor, net)
  population_buc  = gen_Subpopulation_alters_memoryfactor(population_buc, M_vis, sub_mem_factor)
  
  return(population_buc)
}


################################################################################
# Visibility factor estimate #

VF_Estimate = function(enc_hp){
  return(sum(enc_hp$reach_hp_memory))/sum(enc_hp$reach_memory)
}
################################################################################

# Code ideas

# Mhp_hp calculation for making two different Bernouillis, one for the people in the
# hidden population that knows people in the hidden population and for those who are
# not in the hidden population and know people from there.
################################################################################
# M_vis_comp = 2*M_vis - 1
# Mhp_hp = 1*(t(M_vis) ==  M_vis_comp)
# Mhp_hp_vis = apply(Mhp_hp, c(1,2), berHP, p = (0.5 +  0.5*vis_factor))
# Matrix for knowing with bernoilli should be applied
# Mhp_sir = M_vis + 2*Mhp_hp + 2*Mhp_hp_vis
# Final visibility matrix
# M_vis = apply(Mhp_sir, c(1,2), to_matrix_SIR)
############################################################################
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
  update(jagmod, n.iter=5000, progress.bar="none")
  posterior = coda.samples(jagmod, c("alpha","lambda","tau","Su"),n.iter=NITERATION,progress.bar="none",thin=1)
  dicsamples = dic.samples(jagmod,type = "pD",n.iter=20000,thin=1)
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior,dicsamples=dicsamples)
  # uses only the first chain to obtain the hidden population estimation
  Nh = mean(as.matrix(results$posterior[[1]])[,1])
  return(Nh)
}

#########################################
# Zheng et al. (2005) overdispersed model
#########################################

overdispersed_model = "
data {             
  int<lower=0> I;                        // respondents
  int<lower=0> K;                        // subpopulations
  vector[K] mu_beta;                     // prior mean of beta
  //vector<lower=0>[K] sigma_beta;         // prior variance of beta
  int  y[I,K];                           // known by respondent i in subpopulation k
  }

parameters {
  vector[I] alpha;                       // log degree
  vector[K] beta;                        // log prevalence of group in population
  vector<lower = 0 , upper = 1>[K] inv_omega;  // ineverse overdispersion; implies the uniform prior 
  real mu_alpha;                         // prior mean for alpha
  real<lower=0> sigma_alpha;             // prior scale for alpha
  vector<lower=0>[K] sigma_beta;         // prior variance of beta
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



getNh_overdispersed = function(survey, v_pop_total,N, warmup = 2000, iterations = 5000,chains=1){
  y0 = survey %>% dplyr::select(starts_with("kp_reach_")| hp_survey )
  y <- array(dim = c(nrow(y0), ncol(y0)))
  for (i in 1:nrow(y)) {
    for (k in 1:ncol(y)) {
      y[i,k] <- y0[i,k]
    }
  }
  #Inizialization
  b= rep(NA,ncol(y))
  b[1:(ncol(y)-1)] = v_pop_total
  a=rep(NA,nrow(y))
  for (i in 1:length(a)) {
    a[i] = glm(y[i,1:(ncol(y)-1)]~ b[1:(ncol(y)-1)]-1,family = "poisson")$coefficients
  }
  b[ncol(y)] = glm(y[,ncol(y)]~a-1,family = "poisson")$coefficients
  mu_beta = log(b)
  
  data <- list(I = nrow(y), K = ncol(y), mu_beta = mu_beta, y = y)
  
  fit <- stan(model_code = overdispersed_model, #file='NB_norecall.stan', 
              data = data, 
              # warmup = 1000, iter = 2000,  # takes ~ 40min/chain
              warmup = warmup, iter=iterations,
              chains = chains)
  out <- extract(fit)
  beta_post <- out $ beta
  beta_hat <- apply(beta_post, 2, mean)
  Nh = exp(beta_hat[length(beta_hat)])
  return(Nh)
}

################################################################################
####################### Functions for the graphs ###############################

data_analysis = function(Nh_df, Nh_ref_df){
  # Estimation dataframe analysis
  
  df_analysis = data.frame( abserror = rowMeans(as.matrix(abs(Nh_df-Nh_ref_df))),
                            mse      = rowMeans(as.matrix((Nh_df-Nh_ref_df)^2)),
                            bias     = rowMeans(as.matrix(Nh_df)),
                            sd       = rowSds(as.matrix(Nh_df)),
                            median   = rowMedians(as.matrix(Nh_df)) )
  
  return(df_analysis)
}

################################################################################
################# Functions for the network analysis ###########################

net_degree_distribution = function(net, p, nei){
  # Variables #
  size = net$size 
  degree_vect = c()
  degree_loop = rep(NA, size)
  for (j in 1:size){
    degree_loop[j] = length(net[[j]][[1]])
  }
  
  degree_vect = append(degree_vect, degree_loop)
  degree_df = data.frame(degrees = degree_vect)
  
  # Variables of interest #
  degree_var    = round(var(degree_vect), digits = 2)
  degree_mean   = round(mean(degree_vect), digits = 2)
  degree_median = median(degree_vect)
  degree_max    = max(degree_vect)
  degree_min    = min(degree_vect)
  
  # Graph representation #
  sub_title = str_c("Mean = ", degree_mean, ", median = ", degree_median, ", var = ", degree_var,", min = ", degree_min, ", max = ", degree_max, ". Small World model with p = ", p, " and nei = ", nei, ".")
  
  degree_graph = ggplot(degree_df) +
    geom_histogram( aes(x = degrees, y = ..count../sum(..count..)), binwidth = 1, color = "black", fill = "grey", alpha = 0.4) +
    labs(title = "Network degree distribution",
         subtitle = sub_title,
         x = "Number of neighbors",
         y = "Proportion")
  
  return(degree_graph)
}

net_hplinks_distribution = function(net, pop){
  
  final_infected_nodes = as.integer(row.names(pop)[pop$hidden_population == 1])
  links_hp = rep(NA, N)
  for (j in 1:N){
    count = 0
    for (l in net[[j]][[1]]){
      if (as.logical(sum(l %in% final_infected_nodes))){
        count = count + 1
      }
    }
    links_hp[j] = count
  }
  
  # Variables of interest #
  link_var    = round(var(links_hp), digits = 2)
  link_mean   = round(mean(links_hp), digits = 2)
  link_median = median(links_hp)
  link_max    = max(links_hp)
  link_min    = min(links_hp)
  
  sub_title = str_c("Mean = ", link_mean, ", median = ", link_median, ", var = ", link_var,", min = ", link_min, ", max = ", link_max, ". SIR model with beta = 0.115, gamma = 0.0767, 1 hotspot & 5 iterations. HP = ", sum(pop$hidden_population))
  links_graph = ggplot() + 
    geom_line(aes(x = 1:N , y = links_hp)) +
    scale_color_discrete("Legend") + 
    labs(title = "Hidden population distribution",
         subtitle = sub_title,
         x = "People",
         y = "Hidden population links")
  
  return(links_graph)
  
}

net_analysis = function(net, pop, p, nei){
  # Double plot
  plot1 =  net_degree_distribution(net, p, nei)
  plot2 = net_hplinks_distribution(net, pop)
  plt   = grid.arrange(plot1,plot2)
  
  #Variables analysis
  Global_cluster_coefficent = transitivity(net, type = "global")
  Mean_distance = mean_distance(net, weights = NULL, directed = F, unconnected = TRUE, details = FALSE)
  Diameter = diameter(net, directed = F, unconnected = TRUE, weights = NULL)
  Radius = radius(net, mode = "all")
  
  sub_title = str_c("Network parameters: Cluster coefficient = ", round(Global_cluster_coefficent, 2), ", Mean distance = ", round(Mean_distance, 2), ", Diameter = ", Diameter,", Radius = ", Radius) 
  
  title <- ggdraw() + 
    draw_label(sub_title, x = 0.05, fontfamily = "bold", hjust = 0, size = 14)
  
  
  result_graph = plot_grid(title, plt, ncol = 1, rel_heights = c(0.1, 1)
  )
  return(result_graph)
}
