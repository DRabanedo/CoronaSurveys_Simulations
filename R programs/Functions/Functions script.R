
# Functions to create populations and surveys, and the NSUM estimators
############################
library(igraph)
library(igraph)
library(tidyverse)
library(stringr)
library(ggplot2)
library(sampler)


#######################################
# Generation of surveys and populations
#######################################

# This function asigns each individual to a subpopulation
getPopDis <- function(k, pop, p) {
  # Generator of disjoint populations
  # Returns a vector with the population corresponding to each individual
  # k: the number of individuals
  # pop: integer vector representing the different subpopulations
  # p: integer vector with the subpopulations probabilities
  sample(pop, k, replace = TRUE, p = p)
}

# This function uniformly assings if the individuals belong to the Hidden Population
getHiddenPop <- function(k, prob) {
  # Hidden Population generator
  # returns a binary vector, the ones represent the Hidden Population
  # prob: probability of the occurrence of the Hidden Population

  sample(c(0,1), k, replace = TRUE, p = c(1-prob,prob))
}



getStratum <- function(n,StratumProp){
  sample(1:length(StratumProp),n,replace = TRUE, p= StratumProp)
}

# This function generates the Population
genPopulation <- function(n, dis_pop, pop_vect,prob_hidden) {
  # Generates a data frame with the population and the belonging to the Hidden Population
  enc = data.frame(Population = getPopDis(n,dis_pop,pop_vect))
  enc = cbind(enc, Hidden_Population = getHiddenPop(n,prob_hidden))
  rownames(enc) <- c(1:n)
  return(enc)
}

######################
# Matrix for the GNSUM
######################

matrixHP = function(grafo,Pob){
  # adjacency matrix of the directed graph of connections with the Hidden Population
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
  
  n_populations = length(dis_populations)-1
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
    
    vect_reach_re[i] = round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i]))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_memory = vect_reach_re)
  Population = cbind(Population, HP_total = vect_hp) 
  Population = cbind(Population, HP_total_apvis = vect_hp_vis)
  
  for(j in 0:n_populations){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population[net_sw[[i]][[1]],]$Population == j)
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
  
  
  return(list(net_sw, Population, Mhp_vis))
}

getDatos_conEstrat = function(N, poblaciones_dis,v_probabilidades,PropHiddenPop,StratumProp, dim, nei, p, visibility_factor, memory_factor, sub_memory_factor){
  # Lista que contine el grafo con v?rtices etiquetados, y los datos de la poblaci?n
  
  #N es el tamaño de la población que queremos
  #poblaciones_dis vector con todas las poblaciones
  #v_probabilidades vector con las probabilidades de las poblaciones
  #PropHiddenPop proporción de la población oculta
  #StratumProp , vector con las proporciones totales del estrato
  #dim Integer constant, the dimension of the starting lattice.
  #nei Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p Real constant between zero and one, the rewiring probability.
  #memory_factor es el factor de memoria
  #visibility_factor es el factor de visibilidad
  #memory_factor es la proporción que vamos a tomar la varianza de la normal que vamos a aplicar al Reach
  #sub_memory_factor es la proporción que vamos a tomar como varianza de la normal que vamos a aplicar a las subpoblaciones
  
  Pob_general = genPopulation(N, poblaciones_dis,v_probabilidades,PropHiddenPop)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(poblaciones_dis)-1
  
  # Vector con la pertenencia a los diferentes estratos
  Estratos = getEstrato(N,StratumProp)
  # Muestra con reemplazamiento de los factores de correcci?n con probabilidades vect_prob_memory_subpob
  vect_hp = rep(NA,N)       #Vector sin visibility
  vect_hp_vis = rep(NA,N)   #Vector con visibility aplicado
  vect_reach = rep(NA,N)    #Vector del número de nodos que se conecta cada nodo
  vect_reach_re = rep(NA,N) #Vector Reach con error de memoria
  
  # Creamos la matriz del grafo dirigido de personas que conocen a la poblaci?n oculta
  Mhp = matrixHP(net_sw,Pob_general)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = visibility_factor)
  
  
  for (i in 1:N) {
    # En la clase de igraph, [[]] busca vértices adyacentes
    # net_sw[[i]] es un lista con un elemento, la lista con los vértices adyacentes
    
    vect_hp[i] = sum(Mhp[i,])
    vect_reach[i] = length(net_sw[[i]][[1]])
    #  binomial de tamaño las conexiones con la población oculta y probabilidad la visibilidad
    #vect_hp_vis[i] = rbinom(1,vect_hp[i],prob = visibility_factor)
    vect_hp_vis[i] = sum(Mhp_vis[i,])
    
    vect_reach_re[i] = round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i]))
  }
  
  Pob_general = cbind(Pob_general, Estrato = Estratos)
  Pob_general = cbind(Pob_general, Reach = vect_reach)
  Pob_general = cbind(Pob_general, Reach_memory = vect_reach_re)
  Pob_general = cbind(Pob_general, HP_total_conocida = vect_hp) 
  Pob_general = cbind(Pob_general, HP_total_apvis = vect_hp_vis)
  
  for(j in 0:n_populations){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Pob_general[net_sw[[i]][[1]],]$Poblacion == j)
      # Visibilidad de la poblaci?n j por i aplicando un factor de visibilidad para las subpoblaciones
      v_1[i] = round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob))
    }
    
    Pob_general = cbind(Pob_general,Subpoblacion_total = v_1)
    names(Pob_general)[dim(Pob_general)[2]] = str_c("KP_total_apvis_",j)
  }
  
  
  
  
  
  # A?adimos las etiquetas de la pertenencia a la Hidden Population al grafo
  #V(net_sw)$label = Pob_general$Poblacion_Oculta
  
  return(list(net_sw, Pob_general, Mhp_vis))
}


getSurvey = function(n_enc, dataframe){
  #This function makes a sample of size n_enc from dataframe

  sur = dataframe[sample(nrow(dataframe), n_enc, replace = FALSE),]
  
  return(sur)
}

getMuestraEstrat = function(n_enc,dataframe){
  # Muestro estratificado por asignaci?n proporcional
  ssamp(dataframe,n_enc,Estrato)
}
############################


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
  ind2 = as.numeric(rownames(Pob[Pob$Population != 0,]))
  suma = sum(Mhp_vis[ind2,ind1])
  
  denominador = N/sum(v_pob)*suma/nrow(enc_hp)      #Denominator estimate
  
  Nh = numerador/denominador
  return(Nh)
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




