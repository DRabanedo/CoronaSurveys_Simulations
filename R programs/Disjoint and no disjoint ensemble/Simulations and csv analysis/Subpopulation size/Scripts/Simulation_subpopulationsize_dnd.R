######################################################################################
# Simulation based on the size of subpopulations, leaving the rest of parameters fixed
######################################################################################

t = Sys.time()

N = 1000                  # Population size
v_pop = c(1:5)             # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)      # Number of subpopulations
v_pop_prob =  rep(1/10, 5) # Probability of each subpopulation
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
nei = 75   # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
           # before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


# Study parameters
parameters = list(rep(1/10, 9), c(rep(1/20, 8),11/20), c(rep(1/20, 7),6/20,6/20), c(rep(1/20, 4), rep(4/20, 3)), c(rep(1/20, 3), rep(4/20, 4)),  
                  c(rep(1/40, 7),4/20, 4/20, 8/20), c(rep(1/20, 4), 2/20, 3/20, 4/20, 5/20), c(1/20, 2/20, 2/20, 3/20, 4/20, 5/20))


################################################################################

## Population ##

# Not disjoint population #

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix



# Disjoint population #

n_columnas = ncol(Population)

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


################################################################################

## Auxiliar simulation data ##

# Number of simulations
b = 5 

# Variable creation
lista_simulacion = list()
lista_sim = list()

lista_simulacion_disjoint = list()
lista_sim_disjoint = list()

k = length(v_pop)


################################################################################
#Surveys

# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population), n_survey, replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population[Population$Hidden_Population == 1,]), n_survey_hp, replace = FALSE)
}

################################################################################

#Simulation

for (w in 1:length(parameters)) {
  # Parameter selection
  v_pop_prob = parameters[[w]]
  m_pob = length(parameters[[w]])
  
  ## Simulation modifications ##
  n_columnas = ncol(Population)
  v_pop = c(1:m_pob)
  n_pop = length(v_pop)
  
  # Not disjoint populations # 
    
  Population_buc = data.frame(Hidden_Population = Population$Hidden_Population)
  for (i in 1:length(v_pop_prob)) {
    Population_buc = cbind(Population_buc, Subpopulation = sample(c(0,1), N, replace = TRUE, p = c(1-v_pop_prob[i],v_pop_prob[i])))
    names(Population_buc)[dim(Population_buc)[2]] = str_c("Subpopulation_",i)
  }
  
  Population_buc = cbind(Population_buc, Reach = Population$Reach)
  Population_buc = cbind(Population_buc, Reach_memory = Population$Reach_memory)
  Population_buc = cbind(Population_buc, HP_total = Population$HP_total)
  Population_buc = cbind(Population_buc, HP_total_apvis = Population$HP_total_apvis)
  
  Population = Population_buc
    
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population[net_sw[[i]][[1]],][,j+1])
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
    
  k = length(v_pop)
    
  v_pop_total = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total[k] = sum(Population[,k+1]) # N_k
  }
  
  
  
  # Disjoint populations # 
  
  Population_disjoint_buc = data.frame(Hidden_Population = Population_disjoint$Hidden_Population)
  for (i in 1:length(v_pop_prob)) {
    Population_disjoint_buc = cbind(Population_disjoint_buc, Subpopulation = sample(c(0,1), N, replace = TRUE, p = c(1-v_pop_prob[i],v_pop_prob[i])))
    names(Population_disjoint_buc)[dim(Population_disjoint_buc)[2]] = str_c("Subpopulation_",i)
  }
  
  Population_disjoint_buc = cbind(Population_disjoint_buc, Reach = Population_disjoint$Reach)
  Population_disjoint_buc = cbind(Population_disjoint_buc, Reach_memory = Population_disjoint$Reach_memory)
  Population_disjoint_buc = cbind(Population_disjoint_buc, HP_total = Population_disjoint$HP_total)
  Population_disjoint_buc = cbind(Population_disjoint_buc, HP_total_apvis = Population_disjoint$HP_total_apvis)
  
  Population_disjoint = Population_disjoint_buc
  
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population_disjoint[net_sw[[i]][[1]],][,j+1])
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population_disjoint = cbind(Population_disjoint,Subpoblacion_total = v_1)
    names(Population_disjoint)[dim(Population_disjoint)[2]] = str_c("KP_total_apvis_",j)
  }
  
  k = length(v_pop)
  
  v_pop_total_disjoint = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total_disjoint[k] = sum(Population_disjoint[,k+1]) # N_k
  }
  
  
  #Variable reset
  
  Nh_real =  rep(NA,b) 
  
  #Nh_basic_sum = rep(NA,b) 
  #Nh_basicvis_sum = rep(NA,b) 
  #Nh_basic_mean = rep(NA,b) 
  #Nh_basicvis_mean = rep(NA,b)                                      
  
  Nh_PIMLE = rep(NA,b) 
  #Nh_PIMLEvis = rep(NA,b) 
  
  Nh_MLE = rep(NA,b) 
  #Nh_MLEvis = rep(NA,b) 
  
  Nh_MoS = rep(NA,b) 
  #Nh_MoSvis = rep(NA,b) 
  
  Nh_GNSUM = rep(NA,b) 
  
  lista_sim = list() 
  
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_subpop = vf_subpop_es(survey_hp, Population, Mhp_vis, sub_memory_factor)
    
    # Hidden population estimates
    Nh_real = sum(Population$Hidden_Population) 
    
    #Nh_basic_sum    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum = getNh_basicvis_sum(survey,N,vf_subpop) 
    #Nh_basic_mean    = getNh_basic_mean(survey,N) 
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
  
  #Nh_basic_sum_disjoint = rep(NA,b) 
  #Nh_basicvis_sum_disjoint = rep(NA,b) 
  #Nh_basic_mean_disjoint = rep(NA,b) 
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
    vf_subpop = vf_subpop_es(survey_hp, Population_disjoint, Mhp_vis, sub_memory_factor)
    
    #Hidden population estimates
    Nh_real_disjoint = sum(Population_disjoint$Hidden_Population) 
    
    #Nh_basic_sum_disjoint    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum_disjoint = getNh_basicvis_sum(survey,N,vf_subpop) 
    #Nh_basic_mean_disjoint    = getNh_basic_mean(survey,N) 
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
write.csv(simulaciones,                                # Data frame
          file = "Simulations_subpopulationsize_notdisjoint",      # CSV name
          row.names = TRUE )                           # Row names: TRUE or FALSE
################################################################################

################################################################################
write.csv(simulaciones,                                # Data frame
          file = "Simulations_subpopulationsize_disjoint",      # CSV name
          row.names = TRUE )                           # Row names: TRUE or FALSE
################################################################################

timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS #############################

# Computation time (N=1000) (my PC)
#timer -> 2.826157 mins    

# Computation time (N=10000) (office PC)
#timer -> 13.3138 mins

############################################################################