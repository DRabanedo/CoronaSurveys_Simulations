#######################
# MoS estimator of NSUM
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

survey = getSurvey(n_survey,Population)


#Vector with the number of people in each subpopulation
v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}

##### DPLYR sintax (same speed as R sintax) #####
# t <- Sys.time()
# v_pop_total = rep(NA, n_pop)
# v_pop_total = as.vector((Population %>% group_by(Population) %>% summarise(n = n()))[,2])
# Sys.time() - t

################################################################################



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

########################################################################


# Value of estimates
t = Sys.time()
Nh_MoS = getNh_MoS(survey, v_pop_total, N)
Nh_MoS
Sys.time() - t

Nh_MoSvis = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
Nh_MoSvis



# Real value
sum(Population$Hidden_Population)

#################### COMPUTATION TIME ANALYSIS ###########################
# Computation time (N=1000)  (my PC)
#timer ->  0.4987271 secs 
###########################################################################

