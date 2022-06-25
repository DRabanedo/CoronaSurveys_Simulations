########################
# Population generator #
########################

############################
library(igraph)
library(igraph)
library(tidyverse)
library(stringr)
library(ggplot2)
library(sampler)
############################

N = 1000                  # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1     # Visibility factor (Binomial's probability)
seed = 207                # Seed
set.seed(seed)

#Graph
dim = 1    # Graph dimension 
nei = 75   # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


################################################################################


#############################
# Generation of populations #
#############################

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

########################
# Matrix for the GNSUM #
########################

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

################################################################################

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix
