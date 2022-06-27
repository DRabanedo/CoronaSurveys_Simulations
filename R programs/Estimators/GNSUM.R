###########################
# General estimator of NSUM
###########################

############################
library(igraph)
library(igraph)
library(tidyverse)
library(stringr)
library(ggplot2)

N = 1000                  # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
sub_memory_factor = 0     # Subpopulation's memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1     # Visibility factor (Binomial's probability)
seed = 207                # Seed
set.seed(seed)

#Graph
dim = 1    # Graph dimension 
nei = 75   # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
           # before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


#Population and Survey

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

survey = getSurvey(n_survey,Population)    #Survey
survey_hp = getSurvey(n_survey_hp,Population[Population$Hidden_Population==1,]) #Hidden population survey


#Vector with the number of people in each subpopulation

v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population$Population == k) # N_k
  
}


################################

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


#################################


# Value of estimates
t = Sys.time()
Nh_GNSUM = getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
Nh_GNSUM
Sys.time() - t

# Real value
sum(Population$Hidden_Population)

#################### COMPUTATION TIME ANALYSIS ###########################
# Computation time (N=1000)  (my PC)
#timer ->  0.017416 secs 
###########################################################################
