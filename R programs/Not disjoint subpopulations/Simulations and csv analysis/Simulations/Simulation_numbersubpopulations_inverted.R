########################################################################################
# Simulation based on the number of subpopulations, leaving the rest of parameters fixed
########################################################################################

t = Sys.time()

N = 10000                 # Population size
v_pop = c(1:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)   # Number of subpopulations
v_pop_prob = rep(1/length(v_pop), length(v_pop)) # Probability of each subpopulation
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
nei = 75   # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections



# Study parameters
parameters = round(seq(from = 1, to = 30, length.out = 30))


#Population

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]          # Population´s graph
Population = Graph_population_matrix[[2]]      # Population of refeence fixed in each loop
Mhp_vis = Graph_population_matrix[[3]]         # Population's visibility matrix

Population_ref = Population
v_pop_ref = v_pop

# Surveys representing the different iterations. 
# The surveys are fixed so the variance and bias can be calculated.

b = 50 #Number of simulations
lista_simulacion = list()
lista_sim = list()

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population), n_survey, replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population[Population$Hidden_Population == 1,]), n_survey_hp, replace = FALSE)
}

################################################################################


for (w in 1:length(parameters)) {
  k = length(v_pop)
  m_pob = parameters[w]

  ## Simulation modifications ##
  n_columnas = ncol(Population)
  v_pop = c(1:m_pob)
  n_pop = length(v_pop)
  
  
  v_pop_prob = rep(1/length(v_pop), length(v_pop))
  

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
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Hidden population estimates
    Nh_real = sum(Population$Hidden_Population) 
    
    #Nh_basic_sum    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum = getNh_basicvis_sum(survey,N,visibility_factor) 
    #Nh_basic_mean    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean = getNh_basicvis_mean(survey,N,visibility_factor)
    
    Nh_PIMLE    = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
    
    Nh_MLE     = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, visibility_factor)
    
    Nh_MoS     = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
    
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
  print(w)
  
}

simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = parameters)

################################################################################
write.csv(simulaciones,                                # Data frame
          file = "Simulation_subpopulationnumber",     # CSV name
          row.names = TRUE )                           # row names: TRUE or FALSE 

################################################################################

timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (virtual machine)
#timer ->   5.898855 mins  

# Computation time (N=10000) (virtual machine)
#timer ->  1.112204 hours

###########################################################################
