#######################
# MLE estimator of NSUM
#######################

N = 1000                 # Population size
v_pop = c(1:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)   # Number of subpopulations
v_pop_prob = c(0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey

memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
sub_memory_factor = 0        # Subpopulation's memory factor (parameter to change variance of the perturbations' normal)
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

survey = getSurvey(n_survey,Population)    # Survey


#Vector with the number of people in each subpopulation
v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}

##### DPLYR syntax (same speed as R syntax) #####
# t <- Sys.time()
# v_pop_total = rep(NA, n_pop)
# v_pop_total = as.vector((Population %>% group_by(Population) %>% summarise(n = n()))[,2])
# Sys.time() - t

################################################################################

getNh_MLE = function(enc,v_pob) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  
  suma_KP = sum(enc[tail(names(enc),length(v_pob))]) # Known Population sum
  # (\sum y_{iu})/(\frac{\sum N_k}{\sum \sum y_{ik}} )
  Nh_f = sum(enc$HP_total_apvis)*(sum(v_pob)/suma_KP)
  Nh_f
}


##### DPLYR syntax (less velocity than R syntax) #####

#getNh_MLEdplyr = function(enc,v_pob) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulationKnown Population
  #Nh_f = sum(select(enc,HP_total_apvis))*sum(v_pob)/sum(select(enc, starts_with("KP_total_apvis") & -KP_total_apvis0))
  
  #return(Nh_f)
#}



getNh_MLEvis = function(enc,v_pob,vis) {
  #NSUM maximum likelihood estimator(MLE) (Formula from "Thirty Years of the NSUM method")
  #enc: survey
  #v_pob: vector with the number of people in each subpopulation
  #vis: estimation of the visibility factor  
  
  suma_KP = sum(enc[tail(names(enc),length(v_pob))])
  Nh_f = (sum(enc$HP_total_apvis))*(sum(v_pob)/suma_KP)*(1/vis)
  Nh_f
}


################################################################################


# Value of estimates 

t = Sys.time()
Nh_MLE = getNh_MLE(survey, v_pop_total)
Nh_MLE
Sys.time() - t

Nh_MLEvis = getNh_MLEvis(survey, v_pop_total, visibility_factor)
Nh_MLEvis


# Real value

sum(Population$Hidden_Population) 

#################### COMPUTATION TIME ANALYSIS ###########################
# Computation time (N=1000)  (my PC)
#timer ->  0.003318787 secs 
###########################################################################
