
# Functions to create populations and surveys, and the NSUM estimators
############################
library(igraph)
library(igraph)
library(tidyverse)
library(stringr)
library(ggplot2)
library(sampler)
library(dplyr)


#######################################
# Generation of surveys and populations
#######################################

# This function uniformly assings if the individuals belong to the Hidden Population
getHiddenPop <- function(k, prob) {
  # Hidden Population generator
  # returns a binary vector, the ones represent the Hidden Population
  # prob: probability of the occurrence of the Hidden Population

  sample(c(0,1), k, replace = TRUE, p = c(1-prob,prob))
}


# This function generates the Population
genPopulation <- function(n, dis_pop, pop_vect,prob_hidden) {
  # Generates a data frame with the population and the belonging to the Hidden Population
  # n: the number of individuals
  # dis_pop: integer vector representing the different subpopulations
  # pop_vect: vector with the subpopulations probabilities
  # prob_hidden: probability of the occurrence of the Hidden Population
  enc = data.frame(Hidden_Population = getHiddenPop(n,prob_hidden))
  rownames(enc) <- c(1:n)
  for (i in 1:length(pop_vect)) {
    enc = cbind(enc, Subpopulation = sample(c(0,1), n, replace = TRUE, p = c(1-pop_vect[i],pop_vect[i])))
    names(enc)[dim(enc)[2]] = str_c("Subpopulation_",i)
  }

  return(enc)
}

genPopulation_disjoint <- function(n, pop_vect,HP) {
  # Generates a data frame with the population and the belonging to the Hidden Population
  # n: the number of individuals
  # dis_pop: integer vector representing the different subpopulations
  # pop_vect: vector with the subpopulations probabilities
  # HP:  Hidden Population vector
  enc = data.frame("Hidden_Population" = HP)
  
  #Population 0 introduction
  pop_vect = c(1-sum(pop_vect), pop_vect)
  subpop_vector = sample(0:(length(pop_vect)-1), n, replace = TRUE, p = pop_vect)
  
  for (i in 1:(length(pop_vect)-1)) {
    enc = cbind(enc, Subpopulation = as.integer(subpop_vector == i))
    names(enc)[dim(enc)[2]] = str_c("Subpopulation_",i)
  }
  
  return(enc)
}


genPopulation_poisson <- function(n, pop_vect,prob_hidden) {
  # Generates a data frame with the population and the belonging to the Hidden Population
  # n: the number of individuals
  # dis_pop: integer vector representing the different subpopulations
  # pop_vect: vector with the subpopulations probabilities
  # HP:  Hidden Population vector
  enc = data.frame(Hidden_Population = getHiddenPop(n,prob_hidden))
  
  #Population 0 introduction
  pop_vect_num = rep(NA, length(pop_vect))
  for (j in 1:length(pop_vect)) {
    pop_vect_num[j] = rpois(1,pop_vect[j]*n)
  }
  
  for (i in 1:(length(pop_vect))) {
    sam_pop = sample(1:n, size = pop_vect_num[i])
    subpop_vector = rep(0,n)
    for (k in sam_pop){
      subpop_vector[k] = 1
    }
    enc = cbind(enc, Subpopulation = subpop_vector)
    names(enc)[dim(enc)[2]] = str_c("Subpopulation_",i)
  }
  
  return(enc)
}


######################
# Matrix for the GNSUM
######################

matrixHP = function(grafo,Pob){
  # adjacency matrix of the directed graph of connections with the Hidden Population
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

#######################

berHP = function(x,p){
  # Binomial general function for the visibility matrix (element by element)
  # x: matrix element
  # p: binomial probability
  if(x!=0){
    return(x*rbinom(1,1,p))
  }
  else {
    return(0)
  }
}

######################

###############################
# Population generation
###############################

getData = function(N, dis_populations,prob_vector,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor, sub_memory_factor){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  #N: population size
  #dis_populations: vector with the populations
  #prob_vector: vector with the population's probabilities
  #PropHiddenPop: Hidden Population proportion
  #dim: Integer constant, the dimension of the starting lattice.
  #nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p: Real constant between zero and one, the rewiring probability.
  #visibility_factor: the visibility factor
  #memory factor: numeric value, Reach divided by the standard deviation of the normal we use to correct the Reach
  # sub_memory_factor: the subpopulations visibility divided by the standard deviation of the normal distributions we use to correct the subpopulations visibility
  #                    it is applied to each subpopulation
  # sub_memory_factor: the subpopulations' visibility divided by the standard deviation of the normals we use to correct the subpopulations' visibility
  
  Population = genPopulation(N, dis_populations, prob_vector,PropHiddenPop)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(dis_populations)
  # initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # Reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = visibility_factor)
  
  
  for (i in 1:N) {
    # net_sw[[i]], list with one element, the list of the adjacent vertices to i
    
    vect_hp[i] = sum(Mhp[i,])
    vect_reach[i] = length(net_sw[[i]][[1]])
    vect_hp_vis[i] = max(0,round(rnorm(1, mean = sum(Mhp_vis[i,]), sd = memory_factor*sum(Mhp_vis[i,]))))
    
    vect_reach_re[i] = max(1,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_memory = vect_reach_re)
  Population = cbind(Population, HP_total = vect_hp) 
  Population = cbind(Population, HP_total_apvis = vect_hp_vis)
  
  for(j in 1:length(prob_vector)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("Subpop") & ends_with(as.character(j)))) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }

    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
  
  
  return(list(net_sw, Population, Mhp_vis))
}


getData_poisson = function(N, dis_populations,prob_vector,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor, sub_memory_factor){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  #N: population size
  #dis_populations: vector with the populations
  #prob_vector: vector with the population's probabilities
  #PropHiddenPop: Hidden Population proportion
  #dim: Integer constant, the dimension of the starting lattice.
  #nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p: Real constant between zero and one, the rewiring probability.
  #visibility_factor: the visibility factor
  #memory factor: numeric value, Reach divided by the standard deviation of the normal we use to correct the Reach
  # sub_memory_factor: the subpopulations visibility divided by the standard deviation of the normal distributions we use to correct the subpopulations visibility
  #                    it is applied to each subpopulation
  # sub_memory_factor: the subpopulations' visibility divided by the standard deviation of the normals we use to correct the subpopulations' visibility
  
  Population = genPopulation_poisson(N, dis_populations, prob_vector,PropHiddenPop)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(dis_populations)
  # initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # Reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = visibility_factor)
  
  
  for (i in 1:N) {
    # net_sw[[i]], list with one element, the list of the adjacent vertices to i
    
    vect_hp[i] = sum(Mhp[i,])
    vect_reach[i] = length(net_sw[[i]][[1]])
    vect_hp_vis[i] = max(0,round(rnorm(1, mean = sum(Mhp_vis[i,]), sd = memory_factor*sum(Mhp_vis[i,]))))
    
    vect_reach_re[i] = max(1,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_memory = vect_reach_re)
  Population = cbind(Population, HP_total = vect_hp) 
  Population = cbind(Population, HP_total_apvis = vect_hp_vis)
  
  for(j in 1:length(prob_vector)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("Subpop") & ends_with(as.character(j)))) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
  
  
  return(list(net_sw, Population, Mhp_vis))
}


getSurvey = function(n_enc, dataframe){
  #This function makes a sample of size n_enc from dataframe
  # n_enc = number of individuals interviewed
  # dataframe = general population 

  sur = dataframe[sample(nrow(dataframe), n_enc, replace = FALSE),]
  
  return(sur)
}


################################################################################

#####################################################
# Visibility factor prediction using subpopulations #
#####################################################

vf_subpop_es = function(survey_hp,Population, Mhp_vis,sub_memory_factor){
  # Total sums variables creation
  sum_pop = 0
  sum_pop_hp = 0
  
  # Number of subpopulations calculus 
  n_pop = length(names(dplyr::select(Population, starts_with("KP_"))))
  
  # Index of the people who has been surveyed
  ind_survey = as.numeric(rownames(survey_hp))
  
  # Loop for studying each answered value for each subpopulation on the survey 
  for (k in 1:n_pop) {
    
    #People who belong to the subpopulation k
    ind_subpop = as.numeric(rownames(Population[as.logical(dplyr::select(Population, starts_with("Subpop") & ends_with(as.character(k)))[,1]),]))
    
    # Vector representing how many people of subpopulation k knows that the surveyed
    # person belongs to the hidden population
    vect_pop_hp = colSums(Mhp_vis[ind_subpop,ind_survey])
    
    # Application of a memory factor to that answer
    mem_vect_pop_hp = rep(NA,length(vect_pop_hp))
    for (i in 1:length(vect_pop_hp)) {
      mem_vect_pop_hp[i] = max(0,round(rnorm(1,mean = vect_pop_hp[i],sub_memory_factor*vect_pop_hp[i])))
    }
    
    # Vector who represent how many people each surveyed person knows from subpopulation k 
    # (with memory factor applied)
    mem_vect_pop = dplyr::select(Population, starts_with("KP_") & ends_with(as.character(k)))[,1][ind_survey]
    mem_vect_pop
    for (j in 1:length(mem_vect_pop)) {
      
      # As the people a person knows is always bigger that the people that a person knows
      # AND knows its belonging to the hidden population, we make a loop that stops
      # when the assumption commented is verified.
      while (mem_vect_pop_hp[j] > mem_vect_pop[j]) {
        
        mem_vect_pop_hp[j] = max(0,round(rnorm(1,mean = vect_pop_hp[j],sub_memory_factor*vect_pop_hp[j])))
        
        # It makes it converge easier (reduces computation time)
        if (mem_vect_pop_hp[j]<vect_pop_hp[j])
          vect_pop_hp[j] = mem_vect_pop_hp[j]
      }
      
    }
    
    # Sum of the results on subpopulation k to the corresponding general variables
    sum_pop_hp = sum_pop_hp + sum(mem_vect_pop_hp)
    sum_pop = sum_pop + sum(mem_vect_pop)
    
  }
  
  
  # Visibility factor estimate
  vf_subpop = sum_pop_hp/sum_pop
  
  return(vf_subpop)
}


#Outliers detection
#####################################################
# Visibility factor prediction using subpopulations #
#####################################################

vf_subpop_es_out = function(survey_hp,Population, Mhp_vis,sub_memory_factor){
  # Total sums variables creation
  sum_pop = 0
  sum_pop_hp = 0
  
  # Number of subpopulations calculus 
  n_pop = length(names(dplyr::select(Population, starts_with("KP_"))))
  
  # Index of the people who has been surveyed
  ind_survey = as.numeric(rownames(survey_hp))
  
  # Loop for studying each answered value for each subpopulation on the survey 
  for (k in 1:n_pop) {
    
    #People who belong to the subpopulation k
    ind_subpop = as.numeric(rownames(Population[as.logical(dplyr::select(Population, starts_with("Subpop") & ends_with(as.character(k)))[,1]),]))
    
    # Vector representing how many people of subpopulation k knows that the surveyed
    # person belongs to the hidden population
    vect_pop_hp = colSums(Mhp_vis[ind_subpop,ind_survey])
    
    # Application of a memory factor to that answer
    mem_vect_pop_hp = rep(NA,length(vect_pop_hp))
    for (i in 1:length(vect_pop_hp)) {
      mem_vect_pop_hp[i] = max(0,round(rnorm(1,mean = vect_pop_hp[i],sub_memory_factor*vect_pop_hp[i])))
    }
    
    # Vector who represent how many people each surveyed person knows from subpopulation k 
    # (with memory factor applied)
    mem_vect_pop = dplyr::select(Population, starts_with("KP_") & ends_with(as.character(k)))[,1][ind_survey]
    mem_vect_pop
    for (j in 1:length(mem_vect_pop)) {
      
      # As the people a person knows is always bigger that the people that a person knows
      # AND knows its belonging to the hidden population, we eliminate that value
      if (mem_vect_pop_hp[j] > mem_vect_pop[j]){
        mem_vect_pop_hp[j] =  0
        mem_vect_pop[j] = 0
      }
    }
    
    # Sum of the results on subpopulation k to the corresponding general variables
    sum_pop_hp = sum_pop_hp + sum(mem_vect_pop_hp)
    sum_pop = sum_pop + sum(mem_vect_pop)
    
  }
  
  
  # Visibility factor estimate
  vf_subpop = sum_pop_hp/sum_pop
  
  return(vf_subpop)
}


#Outliers detection
#############################################
# Visibility factor predictions using reach #
#############################################

vf_reach_es_out = function(survey_hp,Population, Mhp_vis, memory_factor) {
  # People from the general population who know that the people from the survey_hp belong to the hidden population
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach_hp[i] = max(0,round(rnorm(1,mean = vect_reach_hp[i],memory_factor*vect_reach_hp[i])))
  }
  
  #People from the subpopulations known by the hidden population in survey_hp  
  mem_vect_reach = Population$Reach_memory[ind_survey]
  for (j in 1:length(mem_vect_reach)) {
    
    # As the people a person knows is always bigger that the people that a person knows
    # AND knows its belonging to the hidden population, we eliminate that value
    if (mem_vect_reach_hp[j] > mem_vect_reach[j]){
      mem_vect_reach_hp[j] =  0
      mem_vect_reach[j] = 0
    }
  }
  # Final estimate
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(list(vf_reach))
}


#############################################
# Visibility factor predictions using reach #
#############################################

vf_reach_es = function(survey_hp,Population, Mhp_vis, memory_factor) {
  # People from the general population who know that the people from the survey_hp belong to the hidden population
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach_hp[i] = max(0,round(rnorm(1,mean = vect_reach_hp[i],memory_factor*vect_reach_hp[i])))
  }
  
  #People from the subpopulations known by the hidden population in survey_hp  
  mem_vect_reach = Population$Reach_memory[ind_survey]
  for (i in 1:length(mem_vect_reach)) {
    while (mem_vect_reach[i] < mem_vect_reach_hp[i]){
      mem_vect_reach_hp[i] = max(0,round(rnorm(1,mean = vect_reach_hp[i],memory_factor*vect_reach_hp[i])))
      
      if (mem_vect_reach_hp[i]<vect_reach_hp[i])
        vect_reach_hp[i] = mem_vect_reach_hp[i]
    }
  }
  # Final estimate
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(list(vf_reach))
}

################################################################################

#################
# Basic estimator
#################

getNh_basic_sum = function(survey,N) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  
  Nh_f =  N*sum(survey$HP_total_apvis)/sum(survey$Reach_memory)
  
  return(Nh_f)
}

getNh_basic_mean = function(survey,N) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  
  Nh_f =  N*mean(survey$HP_total_apvis/survey$Reach_memory)
  
  return(Nh_f)
}

getNh_basicvis_sum = function(survey,N,vis) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  #vis: estimation of the visibility factor
  
  Nh_f =  N*sum(survey$HP_total_apvis)/sum(survey$Reach_memory) * (1/vis)
  
  return(Nh_f)
}


getNh_basicvis_mean = function(survey,N,vis) {
  #NSUM Basic estimator  
  #survey: survey
  #N: Population size
  #vis: estimation of the visibility factor
  
  Nh_f =  N*sum(survey$HP_total_apvis)/sum(survey$Reach_memory) * (1/vis)
  
  return(Nh_f)
}


###############
# MLE estimator
###############

getNh_MLE = function(enc,v_pob) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  
  suma_KP = sum(enc[tail(names(enc),length(v_pob))]) # Known Population sum
  # (\sum y_{iu})/(\frac{\sum N_k}{\sum \sum y_{ik}} )
  Nh_f = sum(enc$HP_total_apvis)*(sum(v_pob)/suma_KP)
  Nh_f
}


getNh_MLEvis = function(enc,v_pob,vis) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  #vis: estimation of the visibility factor  
  
  suma_KP = sum(enc[tail(names(enc),length(v_pob))])
  Nh_f = (sum(enc$HP_total_apvis))*(sum(v_pob)/suma_KP)*(1/vis)
  Nh_f
}



#################
# PIMLE estimator
#################


getNh_PIMLE = function(enc,v_pob,N) {
  #NSUM Plug-in Maximum Likelihood Estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  #N: population's size
  #vis: estimation of the visibility factor
  
  #Reach estimate
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  Nh_f = N * mean(enc$HP_total_apvis/d_iest)
  Nh_f
}


getNh_PIMLEvis = function(enc,v_pob,N,vis) {
  #NSUM Plug-in Maximum Likelihood Estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  #N: population's size
  #vis: estimation of the visibility factor
  
  #Reach estimate
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  Nh_f = N * mean(enc$HP_total_apvis/d_iest) * (1/vis) # \frac{y_{iu}}{\hat{d_i}}
  Nh_f
}


###############
# MoS estimator
###############


getNh_MoS = function(enc, v_pob, N){
  #NSUM Mean of Sums(MoS) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  #N: population's size
  
  # \hat{d_i} = N/L \sum_k (y_{ik}/N_k)
  
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$HP_total_apvis/d_i_est)
  Nh_f
}

# Using dplyr in this case increases a lot the complexity of the function

getNh_MoSvis = function(enc, v_pob, N, vis){
  #NSUM Mean of Sums(MoS) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  #N: population's size
  #vis: estimation of the visibility factor 
  
  # \hat{d_i} = N/L \sum_k (y_{ik}/N_k)
  
  # Reach estimate
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$HP_total_apvis/d_i_est) * (1/vis)
  Nh_f
} 


#######
# GNSUM
#######


getNh_GNSUM  = function(Pob, enc, enc_hp, Mhp_vis, v_pob, N, sub_memory_factor){
  #General NSUM (GNSUM) (Formula from "GENERALIZING THE NETWORK SCALE-UP METHOD")
  #Pob:     Population
  #enc:     survey
  #enc_hp:  hidden population's survey
  #Mhp_vis: matrix of the graph connections visibility (see getDatos)
  #v_pob:   vector with the number of people in each subpopulation
  #N:       population's size
  #vis:     estimation of the visibility factor
  
  #Numerator estimate
  n_enc = nrow(enc)
  prob_inc = n_enc/N  #Same inclusion probability for all samples
  numerador = (1/prob_inc) * sum(enc$HP_total_apvis) #Numerator estimate
  
  #Denominator estimate
  ind1 = as.numeric(rownames(enc_hp))
  suma = 0
  for (j in ind1){
    for (i in 1:length(v_pob)) {
      ind2 = Pob[,i+1] != 0
      suma = rnorm(1,sum(Mhp_vis[ind2,j]), sum(Mhp_vis[ind2,j])*sub_memory_factor) + suma
    }
  }
  denominador = N/sum(v_pob)*suma/nrow(enc_hp)      #Denominator estimate
  
  Nh = numerador/denominador
}

##################
# Direct estimator
##################

getNh_Direct = function(survey,N){
  #Direct estimation
  #survey: survey
  #N: Population size
  
  Nh = sum(survey$Hidden_Population)/nrow(survey) * N
  return(Nh)
}
