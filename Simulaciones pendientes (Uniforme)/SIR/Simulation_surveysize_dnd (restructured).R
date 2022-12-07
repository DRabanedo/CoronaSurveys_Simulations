##########################################################################################
# Simulation based on the value of the survey's size, leaving the rest of parameters fixed
##########################################################################################

t = Sys.time()


N = 10000                      # Population size

v_pop_prob =  rep(1/10, 5)    # Probability of each subpopulation
n_pop = length(v_pop_prob)    # Number of subpopulations

hp_prob = 0.1                 # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500                # Number of individuals we draw in the survey
n_survey_hp = 50              # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0         # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0             # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1         # Visibility factor (Binomial's probability)
seed = 207                    # Seed
set.seed(seed)

#Graph
dim = 1    # Graph dimension 
nei = 18   # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections




################################################################################

# Network
net_model = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)

# Not disjoint population #

Graph_population_matrix = gen_Data_SIR(N, v_pop_prob, hp_prob, visibility_factor, memory_factor, sub_memory_factor, net = net_model)

net_sw     = Graph_population_matrix[[1]]      # PopulationÂ´s graph
Population = Graph_population_matrix[[2]]      # Population
Mhp_vis    = Graph_population_matrix[[3]]      # Population's visibility matrix

#Vector with the number of people in each subpopulation

v_pop_total = getV_pop(n_pop, Population)

################################################################################

# Disjoint population #

Population_disjoint = gen_Population_disjoint(N, net_model, v_pop_prob, Population$hidden_population, Mhp_vis, sub_memory_factor, Population$reach, Population$reach_memory, Population$hp_total, Population$hp_survey)

v_pop_total_disjoint =  getV_pop(n_pop, Population_disjoint)
################################################################################

## Auxiliar simulation data ##

# Number of simulations
b = 100 

# Study parameters
parameters    = round(seq(from = 1, to = N, length.out = 100))
parameters_hp = round(seq(from = 1, to = sum(Population$hidden_population), length.out = 100))

simulaciones = data.frame(data = parameters)
simulaciones_disjoint = data.frame(data = parameters)

################################################################################

#Simulation

for (l in 1:b) {
  
  ###########################
  # Not disjoint population #
  
  Nh_real =  rep(NA,length(parameters)) 
  
  Nh_basic_sum    = rep(NA,length(parameters)) 
  #Nh_basicvis_sum = rep(NA,length(parameters))  
  Nh_basic_mean    = rep(NA,length(parameters)) 
  #Nh_basicvis_mean = rep(NA,length(parameters)) 
  
  Nh_PIMLE = rep(NA,length(parameters)) 
  #Nh_PIMLEvis = rep(NA,length(parameters)) 
  
  Nh_MLE = rep(NA,length(parameters)) 
  #Nh_MLEvis = rep(NA,length(parameters)) 
  
  Nh_MoS = rep(NA,length(parameters)) 
  #Nh_MoSvis = rep(NA,length(parameters)) 
  
  Nh_GNSUM = rep(NA,length(parameters))  
  
  Nh_Direct = rep(NA,length(parameters))
  
  
  ########################
  # Disjoint populations #
  
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
  
  Nh_Direct_disjoint = rep(NA,b)
  
  
  for (i in 1:length(parameters)) {
    
    #Surveys variation
    survey_pop = gen_Survey(parameters[i], Population)
    survey = Population[survey_pop,]
    #Hidden's population survey
    survey_hp_pop  = gen_Survey(parameters_hp[i], Population[Population$hidden_population == 1,])
    survey_hp      = Population[Population$hidden_population == 1,][survey_hp_pop,]
    
    #Surveys variation
    survey_disjoint = Population_disjoint[survey_pop,]
    
    survey_hp_disjoint = Population_disjoint[Population_disjoint$hidden_population == 1,][survey_hp_pop,]
    
    
    #Visibility factor estimate
    # Population for the VF estimate
    Population_vf = gen_Survey_VF(sum(Population$hidden_population), Population, Mhp_vis, memory_factor)
    survey_hp_vf  = Population_vf[survey_hp_pop,]
    
    vf_estimate = VF_Estimate(survey_hp_vf)
    vf_estimate_disjoint = vf_estimate
    
    #Estimations
    Nh_real[i] = sum(Population$hidden_population) 
    
    Nh_basic_sum[i]    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum[i] = getNh_basicvis_sum(survey,N,vf_estimate)
    
    Nh_basic_mean[i]    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean[i] = getNh_basicvis_mean(survey,N,vf_estimate)
    
    Nh_PIMLE[i] = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis[i] = getNh_PIMLEvis(survey, v_pop_total, N, vf_estimate)
    
    Nh_MLE[i] = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis[i] = getNh_MLEvis(survey, v_pop_total, vf_estimate)
    
    Nh_MoS[i] = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis[i] = getNh_MoSvis(survey, v_pop_total, N, vf_estimate)
    
    Nh_GNSUM[i] =  getNh_GNSUM(survey, survey_hp, v_pop_total, N)
    
    Nh_Direct[i] = getNh_Direct(survey, N)
    
    
    #Estimations
    Nh_real_disjoint[i] = sum(Population_disjoint$hidden_population) 
    
    Nh_basic_sum_disjoint[i]    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum_disjoint[i] = getNh_basicvis_sum(survey,N,vf_estimate_disjoint)
    
    Nh_basic_mean_disjoint[i]    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean_disjoint[i] = getNh_basicvis_mean(survey,N,vf_estimate_disjoint)
    
    Nh_PIMLE_disjoint[i] = getNh_PIMLE(survey, v_pop_total_disjoint, N)
    #Nh_PIMLEvis_disjoint[i] = getNh_PIMLEvis(survey, v_pop_total_disjoint, N, vf_estimate_disjoint)
    
    Nh_MLE_disjoint[i] = getNh_MLE(survey, v_pop_total_disjoint)
    #Nh_MLEvis_disjoint[i] = getNh_MLEvis(survey, v_pop_total_disjoint, vf_estimate_disjoint)
    
    Nh_MoS_disjoint[i] = getNh_MoS(survey, v_pop_total_disjoint, N)
    #Nh_MoSvis_disjoint[i] = getNh_MoSvis(survey, v_pop_total_disjoint, N, vf_estimate_disjoint)
    
    Nh_GNSUM_disjoint[i] =  getNh_GNSUM(survey, survey_hp, v_pop_total_disjoint, N)
    
    Nh_Direct_disjoint[i] = getNh_Direct(survey, N)
  }
  
  #Dataframe construction
  
  simulaciones = cbind(simulaciones,Nh_real = Nh_real)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_real_",l)
  
  simulaciones = cbind(simulaciones,Nh_basic_sum = Nh_basic_sum)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basic_sum_",l)
  
  #simulaciones = cbind(simulaciones,Nh_basicvis_sum = Nh_basicvis_sum)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basicvis_sum_",l)
  
  simulaciones = cbind(simulaciones,Nh_basic_mean = Nh_basic_mean)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basic_mean_",l)
  
  #simulaciones = cbind(simulaciones,Nh_basicvis_mean = Nh_basicvis_mean)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basicvis_mean_",l)
  
  simulaciones = cbind(simulaciones,Nh_PIMLE = Nh_PIMLE)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_PIMLE_",l)
  
  #simulaciones = cbind(simulaciones,Nh_PIMLEvis = Nh_PIMLEvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_PIMLEvis_",l)
  
  simulaciones = cbind(simulaciones,Nh_MLE = Nh_MLE)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MLE_",l)
  
  #simulaciones = cbind(simulaciones,Nh_MLEvis = Nh_MLEvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MLEvis_",l)
  
  simulaciones = cbind(simulaciones,Nh_MoS = Nh_MoS)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MoS_",l)
  
  #simulaciones = cbind(simulaciones,Nh_MoSvis = Nh_MoSvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MoSvis_",l)
  
  simulaciones = cbind(simulaciones,Nh_GNSUM = Nh_GNSUM)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_GNSUM_",l)
  
  simulaciones = cbind(simulaciones,Nh_Direct = Nh_Direct)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_Direct_",l)
  
  
  #Dataframe construction
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_real = Nh_real)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_real_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basic_sum = Nh_basic_sum_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basic_sum_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basicvis_sum = Nh_basicvis_sum_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basicvis_sum_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basic_mean = Nh_basic_mean_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basic_mean_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basicvis_mean = Nh_basicvis_mean_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basicvis_mean__disjoint",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_PIMLE = Nh_PIMLE_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_PIMLE_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_PIMLEvis = Nh_PIMLEvis_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_PIMLEvis_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MLE = Nh_MLE_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MLE_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MLEvis = Nh_MLEvis_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MLEvis_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MoS = Nh_MoS_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MoS_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MoSvis = Nh_MoSvis_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MoSvis_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_GNSUM = Nh_GNSUM_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_GNSUM_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_Direct = Nh_Direct_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_Direct_",l)
  
  print(l)
  
}

################################################################################
write.csv(simulaciones,                       # Data frame
          file = "Simulation_surveysize_notdisjoint.csv",     # Csv name
          row.names = TRUE )                  # Rownames TRUE o FALSE
################################################################################

################################################################################
write.csv(simulaciones_disjoint,                       # Data frame
          file = "Simulation_surveysize_disjoint.csv",     # Csv name
          row.names = TRUE )                  # Rownames TRUE o FALSE
################################################################################


timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=10000) (my PC)
#timer ->  9.369394 mins

# Computation time (N=10000) (1-1000 // 400)
# timer -> 4.234875 hours
###########################################################################