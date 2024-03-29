################################################################################################################
# Simulation based on the value of the memory factor of the Reach variable, leaving the rest of parameters fixed
################################################################################################################

################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
memory_factor = 0      #Reach memory factor (parameter to change variance of the perturbations' normal)
################################################################################

t = Sys.time()
N = 1000                 # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1     # Visibility factor (Binomial's probability)
seed = 207                # Seed
set.seed(seed)

#Graph
dim = 1    # Graph dimension 
nei = 75   # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections



#Population and Survey
Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

survey = getSurvey(n_survey,Population)
survey_hp = getSurvey(n_survey_hp,Population[Population$Hidden_Population==1,])


#Vector with the number of people in each subpopulation
v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population$Population == k) # N_k
  
}

# Study parameters
parameters = seq(from = 0, to = 0.50, length.out = 41)

#Dataframe to save the data
simulaciones = data.frame(data = parameters)


################################################################################

#Number of simulations
b = 25


for (l in 1:b) {
  Nh_real =  rep(NA,length(parameters)) 
  
  Nh_basic = rep(NA,length(parameters)) 
  #Nh_basicvis = rep(NA,length(parameters)) 
  
  #Nh_PIMLE = rep(NA,length(parameters)) 
  #Nh_PIMLEvis = rep(NA,length(parameters)) 
  
  #Nh_MLE = rep(NA,length(parameters)) 
  #Nh_MLEvis = rep(NA,length(parameters)) 
  
  #Nh_MoS = rep(NA,length(parameters)) 
  #Nh_MoSvis = rep(NA,length(parameters)) 
  
  #Nh_GNSUM = rep(NA,length(parameters))    
  
  
  #General and hidden population's surveys
  survey = getSurvey(n_survey,Population)
  survey_hp = getSurvey(n_survey_hp, Population[Population$Hidden_Population==1,])
  

  # Auxiliar data for the simulation
  vect_reach = survey$Reach
  vect_reach_re =  rep(NA, nrow(survey))
  
  
  for (i in 1:length(parameters)) {
    
    memory_factor = parameters[i]   
    
    for (j in 1:nrow(survey)) {
      vect_reach_re[j] = round(max(rnorm(1,mean = vect_reach[j], sd = memory_factor*vect_reach[j]),1))
    }
    survey$Reach_memory = vect_reach_re
    
    Nh_real[i] = sum(Population$Hidden_Population) 
    
    Nh_basic[i] = getNh_basic(survey,N) 
    #Nh_basicvis[i] = getNh_basicvis(survey,N,visibility_factor) 
    
    #Nh_PIMLE[i] = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis[i] = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
    
    #Nh_MLE[i] = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis[i] = getNh_MLEvis(survey, v_pop_total, visibility_factor)
    
    #Nh_MoS[i] = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis[i] = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
    
    #Nh_GNSUM[i] =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
  }
  
  
  simulaciones = cbind(simulaciones,Nh_real = Nh_real)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_real_",l)
  
  simulaciones = cbind(simulaciones,Nh_basic = Nh_basic)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basic",l)
  
  #simulaciones = cbind(simulaciones,Nh_basicvis = Nh_basicvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basicvis_",l)
  
  #simulaciones = cbind(simulaciones,Nh_PIMLE = Nh_PIMLE)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_PIMLE_",l)
  
  #simulaciones = cbind(simulaciones,Nh_PIMLEvis = Nh_PIMLEvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_PIMLEvis_",l)
  
  #simulaciones = cbind(simulaciones,Nh_MLE = Nh_MLE)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MLE_",l)
  
  #simulaciones = cbind(simulaciones,Nh_MLEvis = Nh_MLEvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MLEvis_",l)
  
  #simulaciones = cbind(simulaciones,Nh_MoS = Nh_MoS)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MoS_",l)
  
  #simulaciones = cbind(simulaciones,Nh_MoSvis = Nh_MoSvis)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MoSvis_",l)
  
  #simulaciones = cbind(simulaciones,Nh_GNSUM = Nh_GNSUM)
  #names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_GNSUM_",l)
  
}

################################################################################


simulaciones
write.csv(simulaciones,                        # Data frame 
          file = "Simulaciones_memoryfactor", # Csv name
          row.names = TRUE )                   # Row names: TRUE or FALSE 


timer = Sys.time() - t
timer


#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC)
#timer -> 11.06498 secs not saving all the unnecessary estimators 
#timer -> 16.75759 mins saving all the unnecessary estimators     

# Computation time (N=10000) (office PC)
#timer ->  8.939757 mins not saving all the unnecessary estimators 
#timer ->  24.34588 mins saving all the unnecessary estimators

#Problem: MoS and PIMLE have computation time of 0.2 per iteration
# 0.2 * 25 * 41 = 200 sec -> 3,33 min
# 3,33 * 4 = 13,3 min

###########################################################################



