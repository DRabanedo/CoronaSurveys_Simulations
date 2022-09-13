########################################################################################
# Simulation based on the size of subpopulations, leaving the rest of parameters fixed #
########################################################################################

t = Sys.time()

N = 10000                 # Population size
v_pop = c(1:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)   # Number of subpopulations
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
nei = 75   # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


# Study parameters
parameters = list(rep(1/10, 9), c(rep(1/20, 8),11/20), c(rep(1/20, 7),6/20,6/20), c(rep(1/20, 4), rep(4/20, 3)), c(rep(1/20, 3), rep(4/20, 4)),  
                  c(rep(1/40, 7),4/20, 4/20, 8/20), c(rep(1/20, 4), 2/20, 3/20, 4/20, 5/20),
                  c(1/20, 2/20, 2/20, 3/20, 4/20, 5/20))

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

for (w in 1:length(parameters)) {
  # Parameter selection
  v_pop_prob = parameters[[w]]
  m_pob = length(parameters[[w]])
  
  ## Simulation modifications ##
  n_columnas = ncol(Population)
  v_pop = c(1:m_pob)
  n_pop = length(v_pop)
  
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
    vf_subpop = vf_subpop_es(survey_hp,Population, Mhp_vis,sub_memory_factor)
    
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
  
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[w]] = simulacion
  
  v_summary_MLE = v_summary_MLE/b
  v_summary_MoS = v_summary_MoS/b
  
  v_summary_MoS_Nh = v_summary_MoS_Nh/b
  v_summary_MLE_Nh = v_summary_MLE_Nh/b
  
  sim_summary_MoS = data.frame("Min" = v_summary_MoS[1], "Max" = v_summary_MoS[2], "Mean" = v_summary_MoS[3], "Median" = v_summary_MoS[4], "sd" = v_summary_MoS[5])
  lista_simulacion_summary_MoS[[w]] = sim_summary_MoS
  
  sim_summary_MLE = data.frame("Min" = v_summary_MLE[1], "Max" = v_summary_MLE[2], "Mean" = v_summary_MLE[3], "Median" = v_summary_MLE[4], "sd" = v_summary_MLE[5])
  lista_simulacion_summary_MLE[[w]] = sim_summary_MLE
  
  sim_summary_MoS_Nh = data.frame("Min" = v_summary_MoS_Nh[1], "Max" = v_summary_MoS_Nh[2], "Mean" = v_summary_MoS_Nh[3], "Median" = v_summary_MoS_Nh[4], "sd" = v_summary_MoS_Nh[5])
  lista_simulacion_summary_MoS_Nh[[w]] = sim_summary_MoS_Nh
  
  sim_summary_MLE_Nh = data.frame("Min" = v_summary_MLE_Nh[1], "Max" = v_summary_MLE_Nh[2], "Mean" = v_summary_MLE_Nh[3], "Median" = v_summary_MLE_Nh[4], "sd" = v_summary_MLE_Nh[5])
  lista_simulacion_summary_MLE_Nh[[w]] = sim_summary_MLE_Nh
  
  print(w)
}

simulaciones_summary_di_MoS = bind_rows(lista_simulacion_summary_MoS)
simulaciones_summary_di_MoS = cbind(simulaciones_summary_di_MoS, data = 1:length(parameters))

simulaciones_summary_di_MLE = bind_rows(lista_simulacion_summary_MLE)
simulaciones_summary_di_MLE = cbind(simulaciones_summary_di_MLE, data = 1:length(parameters))

simulaciones_summary_Nh_MoS = bind_rows(lista_simulacion_summary_MoS_Nh)
simulaciones_summary_Nh_MoS = cbind(simulaciones_summary_Nh_MoS, data = 1:length(parameters))

simulaciones_summary_Nh_MLE = bind_rows(lista_simulacion_summary_MLE_Nh)
simulaciones_summary_Nh_MLE = cbind(simulaciones_summary_Nh_MLE, data = 1:length(parameters))


simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = 1:length(parameters))

################################################################################
write.csv(simulaciones,                                # Data frame
          file = "Simulations_subpopulationsize_notdisjoint",      # CSV name
          row.names = TRUE )                           # Row names: TRUE or FALSE
################################################################################


####################################################################################
write.csv(simulaciones_summary_di_MLE,                             # Data frame    #
          file = "Simulations_subpopulationsize_summary_MLE",      # CSV name      #
          row.names = TRUE )                           # Row names: TRUE or FALSE  #
####################################################################################

####################################################################################
write.csv(simulaciones_summary_di_MoS,                             # Data frame    #
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



timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS #############################

# Computation time (N=1000) (my PC)
#timer -> 2.826157 mins    

# Computation time (N=10000) (office PC)
#timer -> 13.3138 mins

############################################################################