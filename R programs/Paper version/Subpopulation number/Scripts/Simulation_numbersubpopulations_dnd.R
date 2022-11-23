########################################################################################
# Simulation based on the number of subpopulations, leaving the rest of parameters fixed
########################################################################################

t = Sys.time()

N = 10000                 # Population size
v_pop_prob = rep(1/10,5)  # Probability of each subpopulation. sum(v_pop_prob) < 1
n_pop = length(v_pop_prob)

hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 


sub_memory_factor = 0     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1     # Visibility factor (Binomial's probability)

seed = 207                # Seed
set.seed(seed)

#Graph
dim = 1    # Graph dimension 
nei = 18   # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
           # before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


################################################################################
## Populations ##

# Not disjoint population #

Graph_population_matrix = getData(N, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

# Population number
v_pop_total = getV_pop(n_pop, Population)

# Disjoint population #

Population_disjoint =  genPopulation_Disjoint(N, net_sw, v_pop_prob, Population$hidden_population, Mhp_vis, sub_memory_factor, Population$reach, Population$reach_memory, Population$hp_total, Population$hp_survey)

# Population number (disjoint)
v_pop_total_disjoint = getV_pop(n_pop, Population_disjoint)


################################################################################

# Auxiliary simulation data

# Study parameters
parameters = round(seq(from = 1, to = 20, length.out = 10))

#Dataframe to save the data
simulaciones          = data.frame(data = parameters)
simulaciones_disjoint = data.frame(data = parameters)

#Number of simulations
b = 100

#Variable creation
lista_simulacion = list()
lista_sim = list()

lista_simulacion_disjoint = list()
lista_sim_disjoint = list()

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
for (w in 1:length(parameters)) {
  n_pop = parameters[w]
  v_pop_prob = c(rep(0.75/n_pop, n_pop))
  
  population_buc = data.frame(hidden_population = Population$hidden_population)
  
  subpop_vect = round(N * v_pop_prob)
  
  for (k in 1:length(subpop_vect)) {
    # Index belonging to the subpopulation 
    subpop_ind = sample(1:N, subpop_vect[k], replace = FALSE)
    
    # Index transformed into a 0 & 1 vector
    subpop = rep(NA, N)
    for (i in 1:N){
      if (as.logical(sum(i %in% subpop_ind))){
        subpop[i] = 1
      }
      else{
        subpop[i] = 0
      }
      
    }
    
    #Dataframe append
    population_buc = cbind(population_buc, Subpopulation = subpop)
    names(population_buc)[dim(population_buc)[2]] = str_c("subpopulation_",k)
  }
  
  population_buc = cbind(population_buc, reach = Population$reach)
  population_buc = cbind(population_buc, reach_memory = Population$reach_memory)
  population_buc = cbind(population_buc, hp_total = Population$hp_total) 
  population_buc = cbind(population_buc, hp_survey = Population$hp_survey)
  
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(population_buc[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(str_c("_", as.character(j) )))) 
      vis_yij = sum(population_buc[net_sw[[i]][[1]],]["hidden_population"][as.logical(dplyr::select(population_buc[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(str_c("_", as.character(j) )))[,1]),]) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    population_buc = cbind(population_buc,Subpoblacion_total = v_1)
    names(population_buc)[dim(population_buc)[2]] = str_c("kp_reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:length(v_pop_prob)) {
    for (j in ind1){
      ind2 = dplyr::select(population_buc, starts_with("subpop") & ends_with(str_c("_", as.character(i) )))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_memory_factor)) 
    }
    population_buc = cbind(population_buc, Subpoblacion_total = i_hp_vis)
    names(population_buc)[dim(population_buc)[2]] = str_c("kp_alters_",i)
  }
  
  Population = population_buc
  
  # Population number
  v_pop_total = getV_pop(n_pop, Population)

  
  
  # Disjoint population #
  
  population_disjoint_buc = data.frame(hidden_population = Population_disjoint$hidden_population)
  
  subpop_vect = round(N * v_pop_prob)
  rownames(population_disjoint_buc) <- c(1:N)
  
  # Variables for the loop
  sampling_vect = 1:n
  gen_subpop = rep(0, n)
  n_vect = 1:n
  
  for (k in 1:length(subpop_vect)) {
    # Index belonging to the subpopulation k
    subpop_ind = sample(sampling_vect, subpop_vect[k], replace = FALSE)
    
    # Index transformed into a 0 & 1 vector to represent the populations
    subpop = rep(NA, n)
    
    for (i in 1:n){
      if (as.logical(sum(i %in% subpop_ind))){
        subpop[i] = 1
        
        # for k in 1:n appends 1 if a population is assigned
        gen_subpop[i] = 1
      }
      else{
        subpop[i] = 0
      }
      
    }
    
    # Creates a vector with the people who does not have a population 
    sampling_vect = n_vect[gen_subpop == 0]
    
    #Dataframe append population k
    population_disjoint_buc = cbind(population_disjoint_buc, Subpopulation = subpop)
    names(population_disjoint_buc)[dim(population_disjoint_buc)[2]] = str_c("subpopulation_",k)
  }

  
  
  population_disjoint_buc = cbind(population_disjoint_buc, reach = Population$reach)
  population_disjoint_buc = cbind(population_disjoint_buc, reach_memory = Population$reach_memory)
  population_disjoint_buc = cbind(population_disjoint_buc, hp_total = Population$hp_total) 
  population_disjoint_buc = cbind(population_disjoint_buc, hp_survey = Population$hp_survey)
  
  for(j in 1:length(v_pop_prob)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(population_disjoint_buc[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(str_c("_", as.character(j) )))) 
      vis_yij = sum(population_disjoint_buc[net_sw[[i]][[1]],]["hidden_population"][as.logical(dplyr::select(population_disjoint_buc[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(str_c("_", as.character(j) )))[,1]),]) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    population_disjoint_buc = cbind(population_disjoint_buc,Subpoblacion_total = v_1)
    names(population_disjoint_buc)[dim(population_disjoint_buc)[2]] = str_c("kp_reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:length(v_pop_prob)) {
    for (j in ind1){
      ind2 = dplyr::select(population_disjoint_buc, starts_with("subpop") & ends_with(str_c("_", as.character(i) )))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_memory_factor)) 
    }
    population_disjoint_buc = cbind(population_disjoint_buc, Subpoblacion_total = i_hp_vis)
    names(population_disjoint_buc)[dim(population_disjoint_buc)[2]] = str_c("kp_alters_",i)
  }
  
  Population_disjoint = population_disjoint_buc
  
  # Population number (disjoint)
  v_pop_total_disjoint = getV_pop(n_pop, Population_disjoint)
  
  
  ###########################
  ## Not disjoint analysis ##
  
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
  
  # Population for the VF estimate
  Population_vf = getSurvey_VF(sum(Population$hidden_population), Population, Mhp_vis, memory_factor)
  
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$hidden_population == 1,][list_surveys_hp[[l]],]
    survey_hp_vf = Population_vf[list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_estimate = VF_Estimate(survey_hp_vf)
    
    #Hidden population estimates
    Nh_real = sum(Population$hidden_population) 
    
    #Nh_basic_sum    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum = getNh_basicvis_sum(survey,N,vf_estimate) 
    #Nh_basic_mean    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean = getNh_basicvis_mean(survey,N,vf_estimate)
    
    Nh_PIMLE    = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, vf_estimate)
    
    Nh_MLE     = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, vf_estimate)
    
    Nh_MoS     = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, vf_estimate)
    
    Nh_GNSUM   =  getNh_GNSUM(survey, survey_hp, v_pop_total, N)    
    
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
  
  # Population for the visibility factor (vf) estimate
  Population_disjoint_vf = cbind(Population_disjoint[Population_disjoint$hidden_population == 1,], reach_hp = Population_vf$reach_hp)
  Population_disjoint_vf = cbind(Population_disjoint_vf, reach_hp_memory = Population_vf$reach_hp_memory)
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population_disjoint[list_surveys[[l]],]
    survey_hp = Population_disjoint[Population_disjoint$hidden_population == 1,][list_surveys_hp[[l]],]
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
simulaciones = cbind(simulaciones, data = parameters)

simulaciones_disjoint = bind_rows(lista_simulacion_disjoint)
simulaciones_disjoint = cbind(simulaciones_disjoint, data = parameters)




################################################################################
file_name = str_c("Simulation_subpopulationnumber_notdisjoint", seed,".csv")
write.csv(simulaciones,                      # Data frame
          file = file_name,                  # CSV name
          row.names = FALSE )                 # row names: TRUE or FALSE 

################################################################################


################################################################################
file_name_disjoint = str_c("Simulation_subpopulationnumber_disjoint", seed,".csv")
write.csv(simulaciones_disjoint,              # Data frame
          file = file_name_disjoint,          # CSV name
          row.names = FALSE )                  # row names: TRUE or FALSE 

################################################################################

timer = Sys.time() - t
timer

####################### Network analysis #######################################
###### Links to the hidden population distribution & Degree distribution #######
plot_name = str_c("Network_numbersubpopulations_", seed_number, ".png")

png(filename = plot_name,
    width = 1000, height = 1000)
net_analysis(net_sw, Population, p, 2*nei)
dev.off()

#################### COMPUTATION TIME ANALYSIS ###########################
# Computation time (N=10000) (virtual machine)
#timer ->  1.112204 hours
###########################################################################
