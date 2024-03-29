#########################
# Basic estimator of NSUM
#########################


N = 1000                 # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 30            # Number of individuals we draw in the survey

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

################################################################################

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
  
  Nh_f =  N*mean(survey$HP_total_apvis/survey$Reach_memory) * (1/vis)
  
  return(Nh_f)
}


################################################################################


# Value of estimates
t = Sys.time()
Nh_basic = getNh_basic_mean(survey,N) 
Nh_basic
Sys.time() - t

t = Sys.time()
Nh_basic =getNh_basic_sum(survey,N) 
Nh_basic
Sys.time() - t

t = Sys.time()
Nh_basicvis_mean =getNh_basicvis_mean(survey,N,visibility_factor) 
Nh_basicvis_mean
Sys.time() - t

t = Sys.time()
Nh_basicvis_sum =getNh_basicvis_sum(survey,N,visibility_factor) 
Nh_basicvis_sum
Sys.time() - t


# Real value
sum(Population$Hidden_Population) 


#################### COMPUTATION TIME ANALYSIS ###########################
# Computation time (N=1000)  (my PC)
#timer ->  0.003814936 secs 
###########################################################################

