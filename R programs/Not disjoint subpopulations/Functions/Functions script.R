
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
    vect_hp_vis[i] = sum(Mhp_vis[i,])
    
    vect_reach_re[i] = max(1,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_memory = vect_reach_re)
  Population = cbind(Population, HP_total = vect_hp) 
  Population = cbind(Population, HP_total_apvis = vect_hp_vis)
  
  for(j in 1:length(prob_vector)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population[net_sw[[i]][[1]],][,j+1]) # Subpopulation_j column is in (j+1)th dataframe column
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


############################

##############################
# Visibility factor estimate #
##############################

##### Visibility factor estimate (using subpopulations) #####


vf_subpop_es = function(survey_hp,Population, Mhp_vis){
  
  # People from the subpopulations who know that the people from the survey_hp belong to the hidden population
  ind_survey = as.numeric(rownames(survey_hp))
  sum_pop_hp = 0
  for (i in 1:length(v_pop_total)) {
    ind_subpop = Population[,i+1] != 0
    sum_pop_hp = sum(Mhp_vis[ind_subpop,ind_survey]) + sum_pop_hp
  }
  
  #People from the subpopulations known by the hidden population in survey_hp  
  sum_pop = sum(select(Population, starts_with("KP_"))[ind_survey,])
  
  # Final estimate
  vf_subpop = sum_pop_hp/sum_pop
  
  return(vf_subpop)
}

#vf_subpop = vf_subpop_es(survey_hp,Population, Mhp_vis)




##### Visibility factor estimate (using Reach) #####

vf_reach_es = function(survey_hp,Population, Mhp_vis) {
  
  # People from the general population who know that the people from the survey_hp belong to the hidden population
  ind_survey = as.numeric(rownames(survey_hp))
  sum_reach_hp = sum(Mhp_vis[,ind_survey])
  
  #People from the subpopulations known by the hidden population in survey_hp  
  sum_reach = sum(Population$Reach_memory[ind_survey])
  
  # Final estimate
  vf_reach = sum_reach_hp/sum_reach
  return(vf_reach)
}


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


getNh_GNSUM  = function(Pob, enc, enc_hp, Mhp_vis, v_pob, N){
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
  for (i in 1:length(v_pob)) {
    ind2 = Pob[,i+1] != 0
    suma = sum(Mhp_vis[ind2,ind1]) + suma
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




