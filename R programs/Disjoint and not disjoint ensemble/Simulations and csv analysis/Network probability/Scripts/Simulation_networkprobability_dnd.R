#########################################################################################################################
# Simulation based on the value of the probability of aleatorize a graph connection, leaving the rest of parameters fixed
#########################################################################################################################

t = Sys.time()

N = 10000                  # Population size
v_pop = c(1:5)              # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)       # Number of subpopulations
v_pop_prob = rep(1/10, 5)   # Probability of each subpopulation. As we are working with disjoint and no disjoint subpopulations
# sum(v_pop_prob) < 1.
hp_prob = 0.1               # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500              # Number of individuals we draw in the survey
n_survey_hp = 50            # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1     # Visibility factor (Binomial's probability)

seed = 207                # Seed
set.seed(seed)

#Graph
dim = 1     # Graph dimension 
nei = 75    # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1   # Probability of randomize a connection. It is applied to all connections



# Study parameters
parameters = seq(from = 0.05, to = 1, length.out = 50)


################################################################################
# AUXILIARY DATA FOR THE SIMULATION

Population_ref = genPopulation(N, v_pop, v_pop_prob,hp_prob)
Population_disjoint_ref = genPopulation_disjoint(N, v_pop_prob,Population_ref$Hidden_Population)

b = 100 #Number of iterations for the simulation

lista_simulacion = list()
lista_simulacion_disjoint =list()

################################################################################

## Surveys ##

# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population_ref), n_survey, replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population_ref[Population_ref$Hidden_Population == 1,]), n_survey_hp, replace = FALSE)
}

################################################################################

#Simulation
for (w in 1:length(parameters)) {
  p = parameters[w]
  
  # Not disjoint population #
  Population = Population_ref
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_pop = length(v_pop)
  
  # Initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # Reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP, p = visibility_factor)
  
  
  for (h in 1:N) {
    # net_sw[[h]], list with one element, the list of the adjacent vertices to h
    
    vect_hp[h] = sum(Mhp[h,])
    vect_reach[h] = length(net_sw[[h]][[1]])
    vect_hp_vis[h] = sum(Mhp_vis[h,])
    
    vect_reach_re[h] = round(rnorm(1, mean = vect_reach[h], sd = memory_factor*vect_reach[h]))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_memory = vect_reach_re)
  Population = cbind(Population, HP_total = vect_hp) 
  Population = cbind(Population, HP_total_apvis = vect_hp_vis)
  
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population[net_sw[[i]][[1]],][,j+1])
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
  
  #Vector with the number of people in each subpopulation
  
  v_pop_total = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total[k] = sum(Population[,k+1]) # N_k
  }
  
  
  
  # Disjoint population #
  Population_disjoint = Population_disjoint_ref
  
  Population_disjoint = cbind(Population_disjoint, Reach = vect_reach)
  Population_disjoint = cbind(Population_disjoint, Reach_memory = vect_reach_re)
  Population_disjoint = cbind(Population_disjoint, HP_total = vect_hp) 
  Population_disjoint = cbind(Population_disjoint, HP_total_apvis = vect_hp_vis)
  
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population_disjoint[net_sw[[i]][[1]],][,j+1])
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population_disjoint = cbind(Population_disjoint,Subpoblacion_total = v_1)
    names(Population_disjoint)[dim(Population_disjoint)[2]] = str_c("KP_total_apvis_",j)
  }
  
  #Vector with the number of people in each subpopulation
  
  v_pop_total_disjoint = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total_disjoint[k] = sum(Population_disjoint[,k+1]) # N_k
  }
  
  # Not disjoint population analysis #
  
  
  #Vector with the number of people in each subpopulation
  
  v_pop_total = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total[k] = sum(Population[,k+1]) # N_k
  }
  
  #Variable reset
  Nh_real =  rep(NA,b) 
  
  Nh_basic_sum = rep(NA,b) 
  #Nh_basicvis_sum = rep(NA,b) 
  Nh_basic_mean = rep(NA,b) 
  #Nh_basicvis_mean = rep(NA,b)                                      
  
  Nh_PIMLE = rep(NA,b) 
  #Nh_PIMLEvis = rep(NA,b) 
  
  Nh_MLE = rep(NA,b) 
  #Nh_MLEvis = rep(NA,b) 
  
  Nh_MoS = rep(NA,b) 
  #Nh_MoSvis = rep(NA,b) 
  
  Nh_GNSUM = rep(NA,b) 
  
  lista_sim = list()
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_subpop = vf_subpop_es(survey_hp, Population, Mhp_vis, sub_memory_factor)
    
    #Hidden population estimates
    Nh_real = sum(Population$Hidden_Population) 
    
    Nh_basic_sum    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum = getNh_basicvis_sum(survey,N,vf_subpop) 
    Nh_basic_mean    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean = getNh_basicvis_mean(survey,N,vf_subpop) 
    
    Nh_PIMLE    = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, vf_subpop)
    
    Nh_MLE     = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, vf_subpop)
    
    Nh_MoS     = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, vf_subpop)
    
    Nh_GNSUM   =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
    
    
    #Dataframe for saving the estimates
    sim = data.frame(Nh_real = Nh_real)
    names(sim)[dim(sim)[2]] = str_c("Nh_real_",l)
    
    sim = cbind(sim,Nh_basic_sum = Nh_basic_sum)
    names(sim)[dim(sim)[2]] = str_c("Nh_basic_sum_",l)
    
    #sim = cbind(sim,Nh_basicvis_sum = Nh_basicvis_sum)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basicvis_sum_",l)
    
    sim = cbind(sim,Nh_basic_mean = Nh_basic_mean)
    names(sim)[dim(sim)[2]] = str_c("Nh_basic_mean_",l)
    
    #sim = cbind(sim,Nh_basicvis_mean = Nh_basicvis_mean)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basicvis_mean_",l)
    
    sim = cbind(sim,Nh_PIMLE = Nh_PIMLE)
    names(sim)[dim(sim)[2]] = str_c("Nh_PIMLE_",l)
    
    #sim = cbind(sim,Nh_PIMLEvis = Nh_PIMLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_PIMLEvis_",l)
    
    sim = cbind(sim,Nh_MLE = Nh_MLE)
    names(sim)[dim(sim)[2]] = str_c("Nh_MLE_",l)
    
    #sim = cbind(sim,Nh_MLEvis = Nh_MLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MLEvis_",l)
    
    sim = cbind(sim,Nh_MoS = Nh_MoS)
    names(sim)[dim(sim)[2]] = str_c("Nh_MoS_",l)
    
    #sim = cbind(sim,Nh_MoSvis = Nh_MoSvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MoSvis_",l)
    
    sim = cbind(sim,Nh_GNSUM = Nh_GNSUM)
    names(sim)[dim(sim)[2]] = str_c("Nh_GNSUM_",l)
    
    lista_sim[[l]] = sim
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[w]] = simulacion
  
  ######################################
  ## Disjoint subpopulations analysis ##
  
  ## Variable reset ##
  
  Nh_real_disjoint =  rep(NA,b) 
  
  Nh_basic_sum_disjoint = rep(NA,b) 
  #Nh_basicvis_sum_disjoint = rep(NA,b) 
  Nh_basic_mean_disjoint = rep(NA,b) 
  #Nh_basicvis_mean_disjoint = rep(NA,b)                                      
  
  Nh_PIMLE_disjoint = rep(NA,b) 
  #Nh_PIMLEvis_disjoint = rep(NA,b) 
  
  Nh_MLE_disjoint = rep(NA,b) 
  #Nh_MLEvis_disjoint = rep(NA,b) 
  
  Nh_MoS_disjoint = rep(NA,b) 
  #Nh_MoSvis_disjoint = rep(NA,b) 
  
  Nh_GNSUM_disjoint = rep(NA,b) 
  
  lista_sim_disjoint = list()
  
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population_disjoint[list_surveys[[l]],]
    survey_hp = Population_disjoint[Population_disjoint$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_subpop = visibility_factor
    
    #Hidden population estimates
    Nh_real_disjoint = sum(Population_disjoint$Hidden_Population) 
    
    Nh_basic_sum_disjoint    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum_disjoint = getNh_basicvis_sum(survey,N,vf_subpop) 
    Nh_basic_mean_disjoint    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean_disjoint = getNh_basicvis_mean(survey,N,vf_subpop) 
    
    Nh_PIMLE_disjoint    = getNh_PIMLE(survey, v_pop_total_disjoint, N)
    #Nh_PIMLEvis_disjoint = getNh_PIMLEvis(survey, v_pop_total_disjoint, N, vf_subpop)
    
    Nh_MLE_disjoint     = getNh_MLE(survey, v_pop_total_disjoint)
    #Nh_MLEvis_disjoint  = getNh_MLEvis(survey, v_pop_total_disjoint, vf_subpop)
    
    Nh_MoS_disjoint     = getNh_MoS(survey, v_pop_total_disjoint, N)
    #Nh_MoSvis_disjoint  = getNh_MoSvis(survey, v_pop_total_disjoint, N, vf_subpop)
    
    Nh_GNSUM_disjoint   =  getNh_GNSUM(Population_disjoint, survey, survey_hp, Mhp_vis, v_pop_total_disjoint, N)
    
    
    #Dataframe for saving the estimates
    sim_disjoint = data.frame(Nh_real = Nh_real_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_real_",l)
    
    sim_disjoint = cbind(sim_disjoint,Nh_basic_sum = Nh_basic_sum_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_basic_sum_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_basicvis_sum = Nh_basicvis_sum_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_basicvis_sum_",l)
    
    sim_disjoint = cbind(sim_disjoint,Nh_basic_mean = Nh_basic_mean_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_basic_mean_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_basicvis_mean = Nh_basicvis_mean_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_basicvis_mean_",l)
    
    sim_disjoint = cbind(sim_disjoint,Nh_PIMLE = Nh_PIMLE_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_PIMLE_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_PIMLEvis = Nh_PIMLEvis_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_PIMLEvis_",l)
    
    sim_disjoint = cbind(sim_disjoint,Nh_MLE = Nh_MLE_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_MLE_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_MLEvis = Nh_MLEvis_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_MLEvis_",l)
    
    sim_disjoint = cbind(sim_disjoint,Nh_MoS = Nh_MoS_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_MoS_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_MoSvis = Nh_MoSvis_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_MoSvis_",l)
    
    sim_disjoint = cbind(sim_disjoint, Nh_GNSUM = Nh_GNSUM_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_GNSUM_",l)
    
    lista_sim_disjoint[[l]] = sim_disjoint
  }
  simulacion_disjoint = bind_cols(lista_sim_disjoint)
  lista_simulacion_disjoint[[w]] = simulacion_disjoint
  print(w)
  
}

simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = parameters)

simulaciones_disjoint = bind_rows(lista_simulacion_disjoint)
simulaciones_disjoint = cbind(simulaciones_disjoint, data = parameters)



################################################################################

write.csv(simulaciones,                                           # Data frame 
          file = "Simulation_networkprobability_notdisjoint.csv", # Csv name
          row.names = TRUE )                          # Row names: TRUE or FALSE 

################################################################################


################################################################################

write.csv(simulaciones_disjoint,                                # Data frame 
          file = "Simulation_networkprobability_disjoint.csv",  # Csv name
          row.names = TRUE )                       # Row names: TRUE or FALSE 

################################################################################



timer = Sys.time() - t
timer
#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC)
#timer -> 17.22806 mins

# Computation time (N=10000) (virtual machine)
#timer ->  6.721422 hours

###########################################################################



