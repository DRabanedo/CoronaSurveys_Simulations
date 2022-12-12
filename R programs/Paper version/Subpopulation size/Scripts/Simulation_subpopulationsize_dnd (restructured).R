######################################################################################
# Simulation based on the size of subpopulations, leaving the rest of parameters fixed
######################################################################################

t = Sys.time()


N = 10000                  # Population size

v_pop_prob =  rep(1/10, 5) # Probability of each subpopulation
n_pop = length(v_pop_prob) # Number of subpopulations

hp_prob = 0.1              # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500             # Number of individuals we draw in the survey
n_survey_hp = 50           # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0      # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0          # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1      # Visibility factor (Binomial's probability)
seed = 207                 # Seed
set.seed(seed)


#Graph
dim = 1    # Graph dimension 
nei = 18   # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


# Study parameters
parameters = list(rep(0.1,5), c(0.2, 0.1, 0.05, 0.05, 0.05, 0.05), c(0.4, rep(0.025, 4)), rep(0.05, 10), c(0.15, 0.1, 0.05, 0.2), c(0.15, 0.125, 0.1, 0.075, 0.05), rep(0.05, 5), c(0.5, 0.025, 0.025, 0.025, 0.025), c(rep(0.02,10), rep(0.04, 5), 0.08, 0.08, 0.16))

################################################################################

# Network
net_model = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)

## Populations ##
# Not disjoint population #

Graph_population_matrix = gen_Data_SIR(N, v_pop_prob, hp_prob, visibility_factor, memory_factor, sub_memory_factor, net = net_model)

net_sw = Graph_population_matrix[[1]]       # Population´s graph
Population = Graph_population_matrix[[2]]   # Population
Mhp_vis = Graph_population_matrix[[3]]      # Population's visibility matrix

#Vector with the number of people in each subpopulation

v_pop_total = getV_pop(n_pop, Population)

################################################################################

# Disjoint population #

Population_disjoint =  gen_Population_disjoint(N, net_model, v_pop_prob, Population$hidden_population, Mhp_vis, sub_memory_factor, Population$reach, Population$reach_memory, Population$hp_total, Population$hp_survey)

v_pop_total_disjoint =  getV_pop(n_pop, Population_disjoint)
################################################################################

## Auxiliar simulation data ##

# Number of simulations
b = 100 

# Variable creation
lista_simulacion = list()
lista_sim = list()

lista_simulacion_disjoint = list()
lista_sim_disjoint = list()


################################################################################
#Surveys

# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = gen_Survey(n_survey, Population)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = gen_Survey(n_survey_hp, Population[Population$hidden_population == 1,])
}
################################################################################

#Simulation

for (w in 1:length(parameters)) {
  
  # Parameter selection
  v_pop_prob = parameters[[w]]
  n_pop      = length(parameters[[w]])
  
  
  # Not disjoint population loop #
  subpop_df = gen_Subpopulation(N, v_pop_prob)
  
  population_buc  = data.frame(hidden_population = Population$hidden_population) #Hidden population
  population_buc  = cbind(population_buc, subpop_df) #Subpopulations
  population_buc  = cbind(population_buc, reach = Population$reach) #Reach variable
  population_buc  = cbind(population_buc, hp_total = Population$hp_total) # HP reach variable
  population_buc  = cbind(population_buc, hp_survey = Population$hp_survey) # HP reach recall error variable
  population_buc  = cbind(population_buc, reach_memory = Population$reach_memory) #Reach recall error variable
  population_buc  = gen_Subpopulation_memoryfactor(population_buc, Mhp_vis, sub_memory_factor, net_sw)
  population_buc  = gen_Subpopulation_alters_memoryfactor(population_buc, Mhp_vis, sub_memory_factor)
  
  Population = population_buc
  
  #Subpopulation number
  v_pop_total = getV_pop(n_pop, Population)
  
  
  # Disjoint population loop #
  subpop_disjoint_df = gen_Subpopulation(N, v_pop_prob)
  
  population_buc_disjoint  = data.frame(hidden_population = Population_disjoint$hidden_population) #Hidden population
  population_buc_disjoint  = cbind(population_buc_disjoint, subpop_disjoint_df) #Subpopulations
  population_buc_disjoint  = cbind(population_buc_disjoint, reach = Population_disjoint$reach) #Reach variable
  population_buc_disjoint  = cbind(population_buc_disjoint, hp_total = Population_disjoint$hp_total) # HP reach variable
  population_buc_disjoint  = cbind(population_buc_disjoint, hp_survey = Population_disjoint$hp_survey) # HP reach recall error variable
  population_buc_disjoint  = cbind(population_buc_disjoint, reach_memory = Population_disjoint$reach_memory) #Reach recall error variable
  population_buc_disjoint  = gen_Subpopulation_memoryfactor(population_buc_disjoint, Mhp_vis, sub_memory_factor, net_sw)
  population_buc_disjoint  = gen_Subpopulation_alters_memoryfactor(population_buc_disjoint, Mhp_vis, sub_memory_factor)
  
  Population_disjoint = population_buc_disjoint
  
  # Population number (disjoint)
  v_pop_total_disjoint = getV_pop(n_pop, Population_disjoint)
  
  
  #Variable reset
  Nh_real =  rep(NA,b) 
  
  #Nh_basic_sum     = rep(NA,b) 
  #Nh_basicvis_sum  = rep(NA,b) 
  #Nh_basic_mean    = rep(NA,b) 
  #Nh_basicvis_mean = rep(NA,b)                                      
  
  Nh_PIMLE     = rep(NA,b) 
  #Nh_PIMLEvis = rep(NA,b) 
  
  Nh_MLE     = rep(NA,b) 
  #Nh_MLEvis = rep(NA,b) 
  
  Nh_MoS     = rep(NA,b) 
  #Nh_MoSvis = rep(NA,b) 
  
  Nh_GNSUM   = rep(NA,b) 
  
  lista_sim = list() 
  
  # Population for the VF estimate
  Population_vf = gen_Survey_VF(sum(Population$hidden_population), Population, Mhp_vis, memory_factor)
  
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$hidden_population == 1,][list_surveys_hp[[l]],]
    survey_hp_vf = Population_vf[list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_estimate = VF_Estimate(survey_hp_vf)
    
    # Hidden population estimates
    Nh_real = sum(Population$hidden_population) 
    
    #Nh_basic_sum     = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum  = getNh_basicvis_sum(survey,N,vf_estimate) 
    #Nh_basic_mean    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean = getNh_basicvis_mean(survey,N,vf_estimate) 
    
    Nh_PIMLE     = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, vf_estimate)
    
    Nh_MLE      = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, vf_estimate)
    
    Nh_MoS      = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, vf_estimate)
    
    Nh_GNSUM    = getNh_GNSUM(survey, survey_hp, v_pop_total, N)    
    
    
    #Dataframe for saving the estimates
    sim = data.frame(Nh_real = Nh_real)
    names(sim)[dim(sim)[2]] = str_c("Nh_real_",l)
    
    #sim = cbind(sim,Nh_basic_sum = Nh_basic_sum)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basic_sum_",l)
    
    #sim = cbind(sim,Nh_basicvis_sum = Nh_basicvis_sum)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basicvis_sum_",l)
    
    #sim = cbind(sim,Nh_basic_mean = Nh_basic_mean)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basic_mean_",l)
    
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
  
  #Nh_basic_sum_disjoint     = rep(NA,b) 
  #Nh_basicvis_sum_disjoint  = rep(NA,b) 
  #Nh_basic_mean_disjoint    = rep(NA,b) 
  #Nh_basicvis_mean_disjoint = rep(NA,b)                                      
  
  Nh_PIMLE_disjoint     = rep(NA,b) 
  #Nh_PIMLEvis_disjoint = rep(NA,b) 
  
  Nh_MLE_disjoint     = rep(NA,b) 
  #Nh_MLEvis_disjoint = rep(NA,b) 
  
  Nh_MoS_disjoint     = rep(NA,b) 
  #Nh_MoSvis_disjoint = rep(NA,b) 
  
  Nh_GNSUM_disjoint   = rep(NA,b) 
  
  lista_sim_disjoint  = list()
  
  # Population for the visibility factor (vf) estimate
  Population_disjoint_vf = cbind(Population_disjoint[Population_disjoint$hidden_population == 1,], reach_hp = Population_vf$reach_hp)
  Population_disjoint_vf = cbind(Population_disjoint_vf, reach_hp_memory = Population_vf$reach_hp_memory)
  
  #Iterations
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey       = Population_disjoint[list_surveys[[l]],]
    survey_hp    = Population_disjoint[Population_disjoint$hidden_population == 1,][list_surveys_hp[[l]],]
    survey_hp_vf = Population_vf[list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_estimate = VF_Estimate(survey_hp_vf)
    
    #Hidden population estimates
    Nh_real_disjoint = sum(Population_disjoint$hidden_population) 
    
    #Nh_basic_sum_disjoint     = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum_disjoint  = getNh_basicvis_sum(survey,N,vf_estimate) 
    #Nh_basic_mean_disjoint    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean_disjoint = getNh_basicvis_mean(survey,N,vf_estimate) 
    
    Nh_PIMLE_disjoint     = getNh_PIMLE(survey, v_pop_total_disjoint, N)
    #Nh_PIMLEvis_disjoint = getNh_PIMLEvis(survey, v_pop_total_disjoint, N, vf_estimate)
    
    Nh_MLE_disjoint      = getNh_MLE(survey, v_pop_total_disjoint)
    #Nh_MLEvis_disjoint  = getNh_MLEvis(survey, v_pop_total_disjoint, vf_estimate)
    
    Nh_MoS_disjoint      = getNh_MoS(survey, v_pop_total_disjoint, N)
    #Nh_MoSvis_disjoint  = getNh_MoSvis(survey, v_pop_total_disjoint, N, vf_estimate)
    
    Nh_GNSUM_disjoint    = getNh_GNSUM(survey, survey_hp, v_pop_total_disjoint, N)
    
    
    #Dataframe for saving the estimates
    sim_disjoint = data.frame(Nh_real = Nh_real_disjoint)
    names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_real_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_basic_sum = Nh_basic_sum_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_basic_sum_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_basicvis_sum = Nh_basicvis_sum_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_basicvis_sum_",l)
    
    #sim_disjoint = cbind(sim_disjoint,Nh_basic_mean = Nh_basic_mean_disjoint)
    #names(sim_disjoint)[dim(sim_disjoint)[2]] = str_c("Nh_basic_mean_",l)
    
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
simulaciones = cbind(simulaciones, data = 1:length(parameters))

simulaciones_disjoint = bind_rows(lista_simulacion_disjoint)
simulaciones_disjoint = cbind(simulaciones_disjoint, data = 1:length(parameters))

################################################################################
file_name = str_c("Simulations_subpopulationsize_notdisjoint_", seed,".csv")
write.csv(simulaciones,                 # Data frame
          file = file_name,             # CSV name
          row.names = FALSE )           # Row names: TRUE or FALSE
################################################################################

################################################################################
file_name_disjoint = str_c("Simulations_subpopulationsize_disjoint_", seed,".csv")
write.csv(simulaciones,                    # Data frame
          file = file_name_disjoint,       # CSV name
          row.names = FALSE )              # Row names: TRUE or FALSE
################################################################################

timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS #############################
# Computation time (N=10000) (office PC)
#timer -> 
############################################################################