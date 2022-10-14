#################################################################################################
# Simulation based on the number of neighbours of each node, leaving the rest of parameters fixed
#################################################################################################

t = Sys.time()


N = 10000                   # Population size
v_pop_prob = rep(1/10, 5)   # Probability of each subpopulation. As we are working with disjoint and no disjoint subpopulations
# sum(v_pop_prob) < 1.
n_pop = length(v_pop_prob)  # Number of subpopulations

hp_prob = 0.1               # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500              # Number of individuals we draw in the survey
n_survey_hp = 50            # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0       # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0           # reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1       # Visibility factor (Binomial's probability)

seed = 207                # Seed
set.seed(seed)

#Graph
dim = 1     # Graph dimension 
nei = 75    # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1   # Probability of randomize a connection. It is applied to all connections


# Study parameters
parameters = round(seq(from = 10, to = 150, length.out =50))


################################################################################
# AUXILIARY DATA FOR THE SIMULATION

Population_ref = genPopulation(N, v_pop_prob, hp_prob)
Population_Poisson_ref = genPopulation_Poisson_basic(N, v_pop_prob, Population_ref$hidden_population)

b = 100 #Number of iterations for the simulation

lista_simulacion = list()
lista_simulacion_Poisson =list()

################################################################################

## Surveys ##

# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population_ref), n_survey, replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population_ref[Population_ref$hidden_population == 1,]), n_survey_hp, replace = FALSE)
}

################################################################################

#Simulations
for (w in 1:length(parameters)) {
  nei = parameters[w]
  
  # Not disjoint population #
  Population = Population_ref
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_pop = length(v_pop_prob)
  
  # Initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # reach vector applying memory error
  
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
  
  Population = cbind(Population, reach = vect_reach)
  Population = cbind(Population, reach_memory = vect_reach_re)
  Population = cbind(Population, hp_total = vect_hp) 
  Population = cbind(Population, hp_survey = vect_hp_vis)
  
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population[net_sw[[i]][[1]],][,j+1])
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("kp_reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:(length(v_pop_prob))) {
    for (j in ind1){
      ind2 = dplyr::select(Population, starts_with("subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_memory_factor)) 
    }
    Population = cbind(Population, Subpoblacion_total = i_hp_vis)
    names(Population)[dim(Population)[2]] = str_c("kp_alters_",i)
  }
  
  #Vector with the number of people in each subpopulation
  
  v_pop_total = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total[k] = sum(Population[,k+1]) # N_k
  }
  
  
  
  # Disjoint population #
  Population_Poisson = Population_Poisson_ref
  
  Population_Poisson = cbind(Population_Poisson, reach = vect_reach)
  Population_Poisson = cbind(Population_Poisson, reach_memory = vect_reach_re)
  Population_Poisson = cbind(Population_Poisson, hp_total = vect_hp) 
  Population_Poisson = cbind(Population_Poisson, hp_survey = vect_hp_vis)
  
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population_Poisson[net_sw[[i]][[1]],][,j+1])
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population_Poisson = cbind(Population_Poisson,Subpoblacion_total = v_1)
    names(Population_Poisson)[dim(Population_Poisson)[2]] = str_c("kp_reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:(length(v_pop_prob))) {
    for (j in ind1){
      ind2 = dplyr::select(Population_Poisson, starts_with("subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_memory_factor)) 
    }
    Population_Poisson = cbind(Population_Poisson, Subpoblacion_total = i_hp_vis)
    names(Population_Poisson)[dim(Population_Poisson)[2]] = str_c("kp_alters_",i)
  }
  
  #Vector with the number of people in each subpopulation
  
  v_pop_total_Poisson = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total_Poisson[k] = sum(Population_Poisson[,k+1]) # N_k
  }
  
  
  # Visibility factor estimate
  survey_hp_vf = getSurvey_VF(n_survey_hp, Population, Mhp_vis, memory_factor)
  vf_estimate = VF_Estimate(survey_hp_vf)
  
  # Not disjoint population analysis #
  
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
    survey_hp = Population[Population$hidden_population == 1,][list_surveys_hp[[l]],]
    
    
    # Hidden population estimates
    Nh_real = sum(Population$hidden_population) 
    
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
    
    Nh_GNSUM   =  getNh_GNSUM(survey, survey_hp, v_pop_total, N)
    
    
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
  
  Nh_real_Poisson =  rep(NA,b) 
  
  Nh_basic_sum_Poisson = rep(NA,b) 
  #Nh_basicvis_sum_Poisson = rep(NA,b) 
  Nh_basic_mean_Poisson = rep(NA,b) 
  #Nh_basicvis_mean_Poisson = rep(NA,b)                                      
  
  Nh_PIMLE_Poisson = rep(NA,b) 
  #Nh_PIMLEvis_Poisson = rep(NA,b) 
  
  Nh_MLE_Poisson = rep(NA,b) 
  #Nh_MLEvis_Poisson = rep(NA,b) 
  
  Nh_MoS_Poisson = rep(NA,b) 
  #Nh_MoSvis_Poisson = rep(NA,b) 
  
  Nh_GNSUM_Poisson = rep(NA,b) 
  
  lista_sim_Poisson = list()
  
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population_Poisson[list_surveys[[l]],]
    survey_hp = Population_Poisson[Population_Poisson$hidden_population == 1,][list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_subpop = visibility_factor
    
    #Hidden population estimates
    Nh_real_Poisson = sum(Population_Poisson$hidden_population) 
    
    Nh_basic_sum_Poisson    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum_Poisson = getNh_basicvis_sum(survey,N,vf_subpop) 
    Nh_basic_mean_Poisson    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean_Poisson = getNh_basicvis_mean(survey,N,vf_subpop) 
    
    Nh_PIMLE_Poisson    = getNh_PIMLE(survey, v_pop_total_Poisson, N)
    #Nh_PIMLEvis_Poisson = getNh_PIMLEvis(survey, v_pop_total_Poisson, N, vf_subpop)
    
    Nh_MLE_Poisson     = getNh_MLE(survey, v_pop_total_Poisson)
    #Nh_MLEvis_Poisson  = getNh_MLEvis(survey, v_pop_total_Poisson, vf_subpop)
    
    Nh_MoS_Poisson     = getNh_MoS(survey, v_pop_total_Poisson, N)
    #Nh_MoSvis_Poisson  = getNh_MoSvis(survey, v_pop_total_Poisson, N, vf_subpop)
    
    Nh_GNSUM_Poisson   =  getNh_GNSUM(survey, survey_hp, v_pop_total, N)
    
    
    #Dataframe for saving the estimates
    sim_Poisson = data.frame(Nh_real = Nh_real_Poisson)
    names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_real_",l)
    
    sim_Poisson = cbind(sim_Poisson,Nh_basic_sum = Nh_basic_sum_Poisson)
    names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_basic_sum_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_basicvis_sum = Nh_basicvis_sum_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_basicvis_sum_",l)
    
    sim_Poisson = cbind(sim_Poisson,Nh_basic_mean = Nh_basic_mean_Poisson)
    names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_basic_mean_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_basicvis_mean = Nh_basicvis_mean_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_basicvis_mean_",l)
    
    sim_Poisson = cbind(sim_Poisson,Nh_PIMLE = Nh_PIMLE_Poisson)
    names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_PIMLE_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_PIMLEvis = Nh_PIMLEvis_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_PIMLEvis_",l)
    
    sim_Poisson = cbind(sim_Poisson,Nh_MLE = Nh_MLE_Poisson)
    names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MLE_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_MLEvis = Nh_MLEvis_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MLEvis_",l)
    
    sim_Poisson = cbind(sim_Poisson,Nh_MoS = Nh_MoS_Poisson)
    names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MoS_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_MoSvis = Nh_MoSvis_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MoSvis_",l)
    
    sim_Poisson = cbind(sim_Poisson, Nh_GNSUM = Nh_GNSUM_Poisson)
    names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_GNSUM_",l)
    
    lista_sim_Poisson[[l]] = sim_Poisson
  }
  simulacion_Poisson = bind_cols(lista_sim_Poisson)
  lista_simulacion_Poisson[[w]] = simulacion_Poisson
  print(w)
  
}

simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = parameters)

simulaciones_Poisson = bind_rows(lista_simulacion_Poisson)
simulaciones_Poisson = cbind(simulaciones_Poisson, data = parameters)



################################################################################
write.csv(simulaciones,                        # Data frame 
          file = "Simulation_networkneighbours_notdisjoint.csv", # Csv name
          row.names = TRUE )                   # Row names: TRUE or FALSE 
################################################################################

################################################################################
write.csv(simulaciones_Poisson,                        # Data frame 
          file = "Simulation_networkneighbours_Poisson.csv", # Csv name
          row.names = TRUE )                   # Row names: TRUE or FALSE 
################################################################################

timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=10000) (virtual machine)
#timer -> 6.755417 hours
###########################################################################