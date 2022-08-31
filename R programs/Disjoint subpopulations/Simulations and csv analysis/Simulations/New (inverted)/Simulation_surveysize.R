##########################################################################################
# Simulation based on the value of the survey's size, leaving the rest of parameters fixed
##########################################################################################

t = Sys.time()

N = 10000                 # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = rep(1/length(v_pop), length(v_pop)) #Probability of each subpopulation
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
parameters = round(seq(from = 1, to = 1000, length.out = 400))


#Population and Survey generation
Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix


#Vector with the number of people in each subpopulation

v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population$Population == k) # N_k
  
}

################################################################################

simulaciones = data.frame(data = parameters)

b = 50

for (l in 1:b) {

  #Hidden's population survey
  survey_hp = getSurvey(n_survey_hp, Population[Population$Hidden_Population==1,])
  
  #Visibility factor estimate
  vf_subpop = vf_subpop_es(survey_hp, Population, Mhp_vis, sub_memory_factor)
  
  Nh_real =  rep(NA,length(parameters)) 
  
  Nh_basic_sum    = rep(NA,length(parameters)) 
  Nh_basicvis_sum = rep(NA,length(parameters))  
  Nh_basic_mean    = rep(NA,length(parameters)) 
  Nh_basicvis_mean = rep(NA,length(parameters)) 
  
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
    n_survey = parameters[i]
    survey = getSurvey(n_survey,Population)
    
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
    
    Nh_GNSUM[i] =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
    
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
  
  
  print(l)
  
}

################################################################################
write.csv(simulaciones,                       # Data frame
          file = "Simulation_surveysize",     # Csv name
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

