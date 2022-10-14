################################################################################################################
# Simulation based on the value of the memory factor of the reach variable, leaving the rest of parameters fixed
################################################################################################################

t = Sys.time()


#####################
## Simulation data ##
#####################
N = 10000                     # Population size
v_pop_prob = rep(1/10, 5)     # Probability of each subpopulation. As we are working with disjoint and no disjoint subpopulations
# sum(v_pop_prob) < 1. 
n_pop = length(v_pop_prob)    # Number of subpopulations
hp_prob = 0.1                 # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500                # Number of individuals we draw in the survey
n_survey_hp = 50              # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0         # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1         # Visibility factor (Binomial's probability)

################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
memory_factor = 0      #reach memory factor (parameter to change variance of the perturbations' normal)
################################################################################


seed = 207                    # Seed
set.seed(seed)

#Graph
dim = 1      # Graph dimension 
nei = 75     # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1    # Probability of randomize a connection. It is applied to all connections


###############################################################################################################################################################

## Populations ##

# Not disjoint population #

Graph_population_matrix = getData(N, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]       # Population´s graph
Population = Graph_population_matrix[[2]]   # Population
Mhp_vis = Graph_population_matrix[[3]]      # Population's visibility matrix

# Population number
v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(dplyr::select(Population, starts_with("subpop") & ends_with(as.character(k)) ) ) # N_k
}


# Disjoint population #

Population_Poisson =  genPopulation_Poisson(N,v_pop_prob, Population$hidden_population, Mhp_vis, sub_memory_factor, Population$reach, Population$reach_memory, Population$hp_total, Population$hp_survey)


v_pop_total_Poisson = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total_Poisson[k] = sum(dplyr::select(Population_Poisson, starts_with("subpop") & ends_with(as.character(k)) ) ) # N_k
}


################################################################################

## Auxiliary data for the simulation ##

# Study parameters
parameters = seq(from = 0, to = 0.5, length.out = 50)

#Dataframe to save the data
simulaciones = data.frame(data = parameters)
simulaciones_Poisson = data.frame(data = parameters)

# reach vector
vect_reach = Population$reach
vect_hp_vis = Population$hp_survey
vect_reach_re =  rep(NA, nrow(Population))

#Number of iterations for the simulation
b = 100 

lista_simulacion = list()
lista_simulacion_Poisson = list()

# Visibility factor estimate
survey_hp_vf = getSurvey_VF(50, Population, Mhp_vis, memory_factor)

# hp_total vector
vect_hp_vis = Population$hp_total

################################################################################

## Surveys ##

# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population), n_survey, replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population[Population$hidden_population == 1,]), n_survey_hp, replace = FALSE)
}

################################################################################

# Simulation

for (i in 1:length(parameters)) {
  
  ## Parameter implementation ##
  
  memory_factor = parameters[i]   
  
  for (j in 1:nrow(Population)) {
    vect_reach_re[j] =  max(1,round(rtruncnorm(1, a = vect_hp_vis[j] - 0.5, b = 2*vect_reach[j] - vect_hp_vis[j] + 0.5, mean = vect_reach[j], sd = memory_factor*vect_reach[j])))
    
  }
  
  Population$reach_memory = vect_reach_re
  Population_Poisson$reach_memory = vect_reach_re
  
  ##########################################  
  ##   Not disjoint population analysis   ##
  
  ## Variable reset ##
  
  Nh_real =  rep(NA,b) 
  
  Nh_basic_sum = rep(NA,b) 
  #Nh_basicvis_sum = rep(NA,b) 
  Nh_basic_mean = rep(NA,b) 
  #Nh_basicvis_mean = rep(NA,b)                                      
  
  #Nh_PIMLE = rep(NA,b) 
  #Nh_PIMLEvis = rep(NA,b) 
  
  #Nh_MLE = rep(NA,b) 
  #Nh_MLEvis = rep(NA,b) 
  
  #Nh_MoS = rep(NA,b) 
  #Nh_MoSvis = rep(NA,b) 
  
  #Nh_GNSUM = rep(NA,b) 
  
  lista_sim = list()
  
  # Population for the VF estimate
  Population_vf = getSurvey_VF(sum(Population$hidden_population), Population, Mhp_vis, memory_factor)
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$hidden_population == 1,][list_surveys_hp[[l]],]
    survey_hp_vf = Population_vf[list_surveys_hp[[l]],]
    
    # Visibility factor estimate
    vf_estimate = VF_Estimate(survey_hp_vf)
    
    # Hidden population estimates
    Nh_real = sum(Population$hidden_population) 
    
    Nh_basic_sum    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum = getNh_basicvis_sum(survey,N,vf_estimate) 
    Nh_basic_mean    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean = getNh_basicvis_mean(survey,N,vf_estimate) 
    
    #Nh_PIMLE    = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, vf_estimate)
    
    #Nh_MLE     = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, vf_estimate)
    
    #Nh_MoS     = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, vf_estimate)
    
    #Nh_GNSUM   =  getNh_GNSUM(survey, survey_hp, v_pop_total, N)
    
    
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
    
    #sim = cbind(sim,Nh_PIMLE = Nh_PIMLE)
    #names(sim)[dim(sim)[2]] = str_c("Nh_PIMLE_",l)
    
    #sim = cbind(sim,Nh_PIMLEvis = Nh_PIMLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_PIMLEvis_",l)
    
    #sim = cbind(sim,Nh_MLE = Nh_MLE)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MLE_",l)
    
    #sim = cbind(sim,Nh_MLEvis = Nh_MLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MLEvis_",l)
    
    #sim = cbind(sim,Nh_MoS = Nh_MoS)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MoS_",l)
    
    #sim = cbind(sim,Nh_MoSvis = Nh_MoSvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MoSvis_",l)
    
    #sim = cbind(sim,Nh_GNSUM = Nh_GNSUM)
    #names(sim)[dim(sim)[2]] = str_c("Nh_GNSUM_",l)
    
    lista_sim[[l]] = sim
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  
  
  ######################################
  ## Disjoint subpopulations analysis ##
  
  ## Variable reset ##
  
  Nh_real_Poisson =  rep(NA,b) 
  
  Nh_basic_sum_Poisson = rep(NA,b) 
  #Nh_basicvis_sum_Poisson = rep(NA,b) 
  Nh_basic_mean_Poisson = rep(NA,b) 
  #Nh_basicvis_mean_Poisson = rep(NA,b)                                      
  
  #Nh_PIMLE_Poisson = rep(NA,b) 
  #Nh_PIMLEvis_Poisson = rep(NA,b) 
  
  #Nh_MLE_Poisson = rep(NA,b) 
  #Nh_MLEvis_Poisson = rep(NA,b) 
  
  #Nh_MoS_Poisson = rep(NA,b) 
  #Nh_MoSvis_Poisson = rep(NA,b) 
  
  #Nh_GNSUM_Poisson = rep(NA,b) 
  
  lista_sim_Poisson = list()
  
  # Population for the VF estimate
  Population_Poisson_vf = getSurvey_VF(sum(Population_Poisson$hidden_population), Population_Poisson, Mhp_vis, memory_factor)
  
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population_Poisson[list_surveys[[l]],]
    survey_hp = Population_Poisson[Population_Poisson$hidden_population == 1,][list_surveys_hp[[l]],]
    survey_hp_vf = Population_Poisson_vf[list_surveys_hp[[l]],]
    
    
    #Visibility factor estimate
    vf_estimate = VF_Estimate(survey_hp_vf)
    
    #Hidden population estimates
    Nh_real_Poisson = sum(Population_Poisson$hidden_population) 
    
    Nh_basic_sum_Poisson    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum_Poisson = getNh_basicvis_sum(survey,N,vf_estimate) 
    Nh_basic_mean_Poisson    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean_Poisson = getNh_basicvis_mean(survey,N,vf_estimate) 
    
    #Nh_PIMLE_Poisson    = getNh_PIMLE(survey, v_pop_total_Poisson, N)
    #Nh_PIMLEvis_Poisson = getNh_PIMLEvis(survey, v_pop_total_Poisson, N, vf_estimate)
    
    #Nh_MLE_Poisson     = getNh_MLE(survey, v_pop_total_Poisson)
    #Nh_MLEvis_Poisson  = getNh_MLEvis(survey, v_pop_total_Poisson, vf_estimate)
    
    #Nh_MoS_Poisson     = getNh_MoS(survey, v_pop_total_Poisson, N)
    #Nh_MoSvis_Poisson  = getNh_MoSvis(survey, v_pop_total_Poisson, N, vf_estimate)
    
    #Nh_GNSUM_Poisson   =  getNh_GNSUM(survey, survey_hp, v_pop_total, N)
    
    
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
    
    #sim_Poisson = cbind(sim_Poisson,Nh_PIMLE = Nh_PIMLE_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_PIMLE_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_PIMLEvis = Nh_PIMLEvis_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_PIMLEvis_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_MLE = Nh_MLE_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MLE_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_MLEvis = Nh_MLEvis_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MLEvis_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_MoS = Nh_MoS_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MoS_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_MoSvis = Nh_MoSvis_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_MoSvis_",l)
    
    #sim_Poisson = cbind(sim_Poisson,Nh_GNSUM = Nh_GNSUM_Poisson)
    #names(sim_Poisson)[dim(sim_Poisson)[2]] = str_c("Nh_GNSUM_",l)
    
    lista_sim_Poisson[[l]] = sim_Poisson
  }
  simulacion_Poisson = bind_cols(lista_sim_Poisson)
  lista_simulacion_Poisson[[i]] = simulacion_Poisson
  
  print(i)
  
}

simulaciones = bind_rows(lista_simulacion)
simulaciones_Poisson = bind_rows(lista_simulacion_Poisson)

simulaciones["data"] = parameters
simulaciones_Poisson["data"] = parameters


################################################################################
write.csv(simulaciones,                                        # Data frame 
          file = "Simulations_memoryfactor_207.csv",   # Csv name
          row.names = TRUE )                      # Row names: TRUE or FALSE 
################################################################################


################################################################################
# write.csv(simulaciones_Poisson,                           # Data frame 
#           file = "Simulations_memoryfactor_Poisson.csv",  # Csv name
#           row.names = TRUE )                      # Row names: TRUE or FALSE 
################################################################################


timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC) (length(parameters) = 50)
#timer -> 11.14328 mins

# Computation time (N=10000) (office PC) (length(parameters) = 50)
#timer ->  

# Computation time (N=10000) (office PC) (length(parameters) = 50)
#timer ->  

###########################################################################
