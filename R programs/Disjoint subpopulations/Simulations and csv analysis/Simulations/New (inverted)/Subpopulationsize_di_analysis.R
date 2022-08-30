######################################################################################
# Simulation based on the size of subpopulations, leaving the rest of parameters fixed
######################################################################################

N = 10000                 # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob =  rep(1/length(v_pop), length(v_pop)) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1     # Visibility factor (Binomial's probability)
seed = 207                # Seed
set.seed(seed)


#Graph
dim = 1    # Graph dimension 
nei = 75   # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


# Study parameters
parameters = list(rep(1/10, 10), c(rep(1/20, 9),11/20), c(rep(1/20, 8),6/20,6/20), c(rep(4/20, 4), rep(1/20, 4)), c(rep(1/20, 4), rep(4/20, 4)),  
                  c(rep(1/40, 8),4/20, 4/20, 8/20), c(2/20, rep(1/20, 4), 2/20, 3/20, 4/20, 5/20),
                  c(3/20, 1/20, 2/20, 2/20, 3/20, 4/20, 5/20))


#Population and Survey

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

#Auxiliar data for the simulation
b = 50 #Number of simulations
lista_simulacion = list()
lista_sim = list()
k = length(v_pop)

lista_simulacion_summary_MLE = list()
lista_simulacion_summary_MoS = list()

lista_simulacion_summary_MLE_Nh = list()
lista_simulacion_summary_MoS_Nh = list()

#Surveys
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

for (i in 1:length(parameters)) {
  # Parameter selection
  v_pop_prob = parameters[[i]]
  m_pob = length(parameters[[i]])-1
  
  ## Simulation modifications ##
  n_columnas = ncol(Population)
  v_pop = c(0:m_pob)
  n_pop = length(v_pop)-1 
  
  Population$Population = sample(v_pop, N, replace = TRUE, p = v_pop_prob)
  Population = Population[,1:(ncol(Population)-k)]
  
  for(j in 0:n_pop){
    v_1 = rep(NA,N)
    for(v in 1:N) {
      vis_pob = sum(Population[net_sw[[v]][[1]],]$Population == j)
      v_1[v] = round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob))
    }
    
    Population = cbind(Population,SubPopulation_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
  
  k = length(v_pop)
  
  v_pop_total = rep(NA, n_pop)
  for (j in 1:n_pop) {
    v_pop_total[j] = sum(Population$Population == j)
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
  
  v_summary_MoS = rep(0,5)
  v_summary_MLE = rep(0,5)
  
  v_summary_MoS_Nh = rep(0,5)
  v_summary_MLE_Nh = rep(0,5)
  
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_subpop = vf_subpop_es(survey_hp, Population, Mhp_vis,sub_memory_factor)
    
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
    
    # Degrees
    
    # Real reach
    d_real = survey$Reach
    
    # Reach applying the memory fector
    d_basic = survey$Reach_memory
    
    # Estimated Reach by MLE and PIMLE
    d_iest = c()
    for (j in 1:nrow(survey)) {
      d_iest[j] = N * sum(survey[j,tail(names(survey),length(v_pop_total))])/sum(v_pop_total)
    }
    d_MLE = d_iest 
    
    # Estimated Reach by MoS
    d_i_est = rep(NA, nrow(survey))
    for (j in 1:nrow(survey)) {
      d_i_est[j] = (sum((survey[j,tail(names(survey),length(v_pop_total))])/v_pop_total))/length(v_pop_total) * N
    }
    d_MoS=d_i_est
    
    
    #Estimated Nh for MoS
    Nh_MoS_analysis = survey$HP_total_apvis/d_i_est
    
    #Estimated Nh for MLE
    Nh_MLE_analysis = survey$HP_total_apvis/d_iest
    
    Nh_real_analysis = sum(Population$Hidden_Population)/N
    
    v_summary_MoS = v_summary_MoS + c(min(abs(d_real - d_MoS)),max(abs(d_real - d_MoS)),mean(abs(d_real - d_MoS)),median(abs(d_real - d_MoS)), sd(abs(d_real - d_MoS)))
    v_summary_MLE = v_summary_MLE + c(min(abs(d_real - d_MLE)),max(abs(d_real - d_MLE)),mean(abs(d_real - d_MLE)),median(abs(d_real - d_MLE)), sd(abs(d_real - d_MLE)))
    
    v_summary_MoS_Nh = v_summary_MoS_Nh + c(min(N*abs(Nh_MoS_analysis-Nh_real_analysis)),max(N*abs(Nh_MoS_analysis-Nh_real_analysis)),mean(N*abs(Nh_MoS_analysis-Nh_real_analysis)),median(N*abs(Nh_MoS_analysis-Nh_real_analysis)), sd(N*abs(Nh_MoS_analysis-Nh_real_analysis)))
    v_summary_MLE_Nh = v_summary_MLE_Nh + c(min(N*abs(Nh_MLE_analysis-Nh_real_analysis)),max(N*abs(Nh_MLE_analysis-Nh_real_analysis)),mean(N*abs(Nh_MLE_analysis-Nh_real_analysis)),median(N*abs(Nh_MLE_analysis-Nh_real_analysis)), sd(N*abs(Nh_MLE_analysis-Nh_real_analysis)))
    
  }
  
  v_summary_MLE = v_summary_MLE/b
  v_summary_MoS = v_summary_MoS/b
  
  v_summary_MoS_Nh = v_summary_MoS_Nh/b
  v_summary_MLE_Nh = v_summary_MLE_Nh/b
  
  sim_summary_MoS = data.frame("Min" = v_summary_MoS[1], "Max" = v_summary_MoS[2], "Mean" = v_summary_MoS[3], "Median" = v_summary_MoS[4], "sd" = v_summary_MoS[5])
  lista_simulacion_summary_MoS[[i]] = sim_summary_MoS
  
  sim_summary_MLE = data.frame("Min" = v_summary_MLE[1], "Max" = v_summary_MLE[2], "Mean" = v_summary_MLE[3], "Median" = v_summary_MLE[4], "sd" = v_summary_MLE[5])
  lista_simulacion_summary_MLE[[i]] = sim_summary_MLE
  
  sim_summary_MoS_Nh = data.frame("Min" = v_summary_MoS_Nh[1], "Max" = v_summary_MoS_Nh[2], "Mean" = v_summary_MoS_Nh[3], "Median" = v_summary_MoS_Nh[4], "sd" = v_summary_MoS_Nh[5])
  lista_simulacion_summary_MoS_Nh[[w]] = sim_summary_MoS_Nh
  
  sim_summary_MLE_Nh = data.frame("Min" = v_summary_MLE_Nh[1], "Max" = v_summary_MLE_Nh[2], "Mean" = v_summary_MLE_Nh[3], "Median" = v_summary_MLE_Nh[4], "sd" = v_summary_MLE_Nh[5])
  lista_simulacion_summary_MLE_Nh[[w]] = sim_summary_MLE_Nh
  
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  
  print(i)
  
}

simulaciones_summary_MoS = bind_rows(lista_simulacion_summary_MoS)
simulaciones_summary_MoS = cbind(simulaciones_summary_MoS, data = 1:length(parameters))

simulaciones_summary_MLE = bind_rows(lista_simulacion_summary_MLE)
simulaciones_summary_MLE = cbind(simulaciones_summary_MLE, data = 1:length(parameters))

simulaciones_summary_Nh_MoS = bind_rows(lista_simulacion_summary_MoS_Nh)
simulaciones_summary_Nh_MoS = cbind(simulaciones_summary_Nh_MoS, data = 1:length(parameters))

simulaciones_summary_Nh_MLE = bind_rows(lista_simulacion_summary_MLE_Nh)
simulaciones_summary_Nh_MLE = cbind(simulaciones_summary_Nh_MLE, data = 1:length(parameters))

simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = 1:length(parameters))

####################################################################################
write.csv(simulaciones,                                # Data frame                #
          file = "Simulations_subpopulationsize",      # CSV name                  #
          row.names = TRUE )                           # Row names: TRUE or FALSE  #
####################################################################################

####################################################################################
write.csv(simulaciones_summary_MLE,                                # Data frame    #
          file = "Simulations_subpopulationsize_summary_MLE",      # CSV name      #
          row.names = TRUE )                           # Row names: TRUE or FALSE  #
####################################################################################

####################################################################################
write.csv(simulaciones_summary_MoS,                                # Data frame    #
          file = "Simulations_subpopulationsize_summary_MoS",      # CSV name      #
          row.names = TRUE )                           # Row names: TRUE or FALSE  #
####################################################################################

####################################################################################
write.csv(simulaciones_summary_Nh_MLE,                                # Data frame #
          file = "Simulations_subpopulationsize_summary_MLE_Nh",      # CSV name   #
          row.names = TRUE )                           # Row names: TRUE or FALSE  #
####################################################################################

####################################################################################
write.csv(simulaciones_summary_Nh_MoS,                                # Data frame #
          file = "Simulations_subpopulationsize_summary_MoS_Nh",      # CSV name   #
          row.names = TRUE )                           # Row names: TRUE or FALSE  #
####################################################################################
