########################################################################################
# Simulation based on the number of subpopulations, leaving the rest of parameters fixed
########################################################################################


t = Sys.time()

N = 1000                  # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
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
parameters = round(seq(from = 2, to = 20, length.out = 10))


#Population

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]          # Population´s graph
Population = Graph_population_matrix[[2]]  # Population of refeence fixed in each loop
Mhp_vis = Graph_population_matrix[[3]]         # Population's visibility matrix

Population_ref = Population
v_pop_ref = v_pop

# Surveys representing the different iterations. 
# The surveys are fixed so the variance and bias can be calculated.

b = 25 #Number of simulations

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population), n_survey, replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population[Population$Hidden_Population == 1,]), n_survey_hp, replace = FALSE)
}


# Dataframe
simulaciones = data.frame(data = parameters)
  
################################################################################

for (l in 1:25) {
  
  Nh_real =  rep(NA,length(parameters)) 
  
  Nh_basic = rep(NA,length(parameters)) 
  Nh_basicvis = rep(NA,length(parameters)) 
  
  Nh_PIMLE = rep(NA,length(parameters)) 
  Nh_PIMLEvis = rep(NA,length(parameters)) 
  
  Nh_MLE = rep(NA,length(parameters)) 
  Nh_MLEvis = rep(NA,length(parameters)) 
  
  Nh_MoS = rep(NA,length(parameters)) 
  Nh_MoSvis = rep(NA,length(parameters)) 
  
  Nh_GNSUM = rep(NA,length(parameters)) 

  
  Population = Population_ref
  v_pop = v_pop_ref
  #Aux variable to select the right cols of the dataframe
  k = length(v_pop)
  
  
  for (i in 1:length(parameters)) {
    m_pob = parameters[i]
    n_columnas = ncol(Population)
    v_pop = c(0:m_pob)
    n_pop = length(v_pop)-1 
    v_pop_prob = rep(1/length(v_pop), length(v_pop))
    print(ncol(Population))
    Population$Population = sample(v_pop, N, replace = TRUE, p = v_pop_prob)
    Population = Population[,1:(ncol(Population)-k)]
    print(ncol(Population))
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
    print(ncol(Population))
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Estimate
    #Nh_real[i] = sum(Population$Hidden_Population) 
    
    #Nh_basic[i] = getNh_basic(survey,N) 
    #Nh_basicvis[i] = getNh_basicvis(survey,N,visibility_factor) 
    
    Nh_PIMLE[i] = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis[i] = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
    
    Nh_MLE[i] = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis[i] = getNh_MLEvis(survey, v_pop_total, visibility_factor)
    
    Nh_MoS[i] = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis[i] = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
    
    Nh_GNSUM[i] =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
  }
  
  
  simulaciones = cbind(simulaciones,Nh_real = Nh_real)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_real_",l)
  
  #simulaciones = cbind(simulaciones,Nh_basic = Nh_basic)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basic",l)
  
  #simulaciones = cbind(simulaciones,Nh_basicvis = Nh_basicvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basicvis_",l)
  
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
  
  
  
  print(l)
  
}


################################################################################
simulaciones
write.csv(simulaciones,                                # Data frame
          file = "Simulaciones_numerosubPopulationes", # CSV name
          row.names = TRUE )                           # row names: TRUE or FALSE 

################################################################################

timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC)
#timer ->  41.22529 mins  

# Computation time (N=10000) (my PC)
#timer ->  

###########################################################################

