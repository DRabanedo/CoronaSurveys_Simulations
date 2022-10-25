##########################################################################################
# Simulation based on the value of the survey's size, leaving the rest of parameters fixed
##########################################################################################

t = Sys.time()

N = 10000                 # Population size
v_pop = c(1:5)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)   # Number of subpopulations
v_pop_prob = rep(1/10, 5) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

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
parameters = round(seq(from = 1, to = 1000, length.out = 50))

################################################################################

# Not disjoint population #

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

#Vector with the number of people in each subpopulation

v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}

################################################################################

# Disjoint population #

Population_disjoint =  genPopulation_disjoint(N,v_pop_prob, Population$Hidden_Population)

Population_disjoint = cbind(Population_disjoint, Reach = Population$Reach)
Population_disjoint = cbind(Population_disjoint, Reach_memory = Population$Reach_memory)
Population_disjoint = cbind(Population_disjoint, HP_total = Population$HP_total)
Population_disjoint = cbind(Population_disjoint, HP_total_apvis = Population$HP_total_apvis)


for(j in 1:length(v_pop_prob)){
  v_1 = rep(NA,N)
  for(i in 1:N) {
    vis_pob = sum(Population_disjoint[net_sw[[i]][[1]],][,j+1])
    v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
  }
  
  Population_disjoint = cbind(Population_disjoint,Subpoblacion_total = v_1)
  names(Population_disjoint)[dim(Population_disjoint)[2]] = str_c("KP_total_apvis_",j)
}

v_pop_total_disjoint = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total_disjoint[k] = sum(Population_disjoint[,k+1]) # N_k
}
################################################################################

## Auxiliar simulation data ##

# Number of simulations
b = 100 

simulaciones = data.frame(data = parameters)
simulaciones_disjoint = data.frame(data = parameters)

################################################################################
#Surveys

# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:length(parameters)) {
  list_surveys[[h]] = sample(nrow(Population), parameters[h] , replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population[Population$Hidden_Population == 1,]), n_survey_hp, replace = FALSE)
}

################################################################################

#Simulation

for (l in 1:b) {

  ###########################
  # Not disjoint population #
  
  #Hidden's population survey
  survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
  
  #Visibility factor estimate
  vf_subpop = visibility_factor
  
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
  
  
  for (i in 1:length(parameters)) {
    
    #Surveys variation
    survey = Population[list_surveys[[i]],]
    
    #Estimations
    Nh_real[i] = sum(Population$Hidden_Population) 
    
    Nh_basic_sum[i]    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum[i] = getNh_basicvis_sum(survey,N,vf_subpop)
    
    Nh_basic_mean[i]    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean[i] = getNh_basicvis_mean(survey,N,vf_subpop)
    
    Nh_PIMLE[i] = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis[i] = getNh_PIMLEvis(survey, v_pop_total, N, vf_subpop)
    
    Nh_MLE[i] = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis[i] = getNh_MLEvis(survey, v_pop_total, vf_subpop)
    
    Nh_MoS[i] = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis[i] = getNh_MoSvis(survey, v_pop_total, N, vf_subpop)
    
    Nh_GNSUM[i] =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N, sub_memory_factor)
    
    Nh_Direct[i] = getNh_Direct(survey, N)
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
  
  
  ########################
  # Disjoint populations #
  
  survey_hp = Population_disjoint[Population_disjoint$Hidden_Population == 1,][list_surveys_hp[[l]],]
  
  
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
    survey = Population_disjoint[list_surveys[[i]],]
    
    #Estimations
    Nh_real_disjoint[i] = sum(Population_disjoint$Hidden_Population) 
    
    Nh_basic_sum_disjoint[i]    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum_disjoint[i] = getNh_basicvis_sum(survey,N,vf_subpop)
    
    Nh_basic_mean_disjoint[i]    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean_disjoint[i] = getNh_basicvis_mean(survey,N,vf_subpop)
    
    Nh_PIMLE_disjoint[i] = getNh_PIMLE(survey, v_pop_total_disjoint, N)
    #Nh_PIMLEvis_disjoint[i] = getNh_PIMLEvis(survey, v_pop_total_disjoint, N, vf_subpop)
    
    Nh_MLE_disjoint[i] = getNh_MLE(survey, v_pop_total_disjoint)
    #Nh_MLEvis_disjoint[i] = getNh_MLEvis(survey, v_pop_total_disjoint, vf_subpop)
    
    Nh_MoS_disjoint[i] = getNh_MoS(survey, v_pop_total_disjoint, N)
    #Nh_MoSvis_disjoint[i] = getNh_MoSvis(survey, v_pop_total_disjoint, N, vf_subpop)
    
    Nh_GNSUM_disjoint[i] =  getNh_GNSUM(Population_disjoint, survey, survey_hp, Mhp_vis, v_pop_total_disjoint, N, sub_memory_factor)
    
    Nh_Direct_disjoint[i] = getNh_Direct(survey, N)
    
  }
  
  #Dataframe construction
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_real_disjoint = Nh_real_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_real_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basic_sum_disjoint = Nh_basic_sum_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basic_sum_disjoint_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basicvis_sum_disjoint = Nh_basicvis_sum_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basicvis_sum_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basic_mean_disjoint = Nh_basic_mean_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basic_mean_disjoint",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_basicvis_mean_disjoint = Nh_basicvis_mean_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_basicvis_mean__disjoint",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_PIMLE_disjoint = Nh_PIMLE_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_PIMLE_disjoint_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_PIMLEvis_disjoint = Nh_PIMLEvis_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_PIMLEvis_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MLE_disjoint = Nh_MLE_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MLE_disjoint_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MLEvis_disjoint = Nh_MLEvis_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MLEvis_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MoS_disjoint = Nh_MoS_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MoS_disjoint_",l)
  
  #simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_MoSvis_disjoint = Nh_MoSvis_disjoint)
  #names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_MoSvis_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_GNSUM_disjoint = Nh_GNSUM_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_GNSUM_disjoint_",l)
  
  simulaciones_disjoint = cbind(simulaciones_disjoint,Nh_Direct_disjoint = Nh_Direct_disjoint)
  names(simulaciones_disjoint)[dim(simulaciones_disjoint)[2]] = str_c("Nh_Direct_disjoint_",l)

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

# Computation time (N=1000) (my PC)
#timer ->  1.511148 mins   

# Computation time (N=10000) (my PC)
#timer ->  9.369394 mins

# Computation time (N=10000) (1-1000 // 400)
# timer -> 4.234875 hours
###########################################################################