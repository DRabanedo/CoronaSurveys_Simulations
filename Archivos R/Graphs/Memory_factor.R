###########################################################################################################
# Graph based on the value of the memory factor of the Reach variable, leaving the rest of parameters fixed
###########################################################################################################
t = Sys.time()
################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
memory_factor = 0      #Reach memory factor (parameter to change variance of the perturbations' normal)
################################################################################


N = 1000                  # Population size
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


# Variable reset
Nh_real =  rep(NA,length(parameters)) 

Nh_basic_sum = rep(NA,length(parameters)) 
Nh_basic_mean = rep(NA,length(parameters)) 
Nh_basicvis_sum = rep(NA,length(parameters)) 
Nh_basicvis_mean = rep(NA,length(parameters)) 

#Nh_PIMLE = rep(NA,length(parameters)) 
#Nh_PIMLEvis = rep(NA,length(parameters)) 

#Nh_MLE = rep(NA,length(parameters)) 
#Nh_MLEvis = rep(NA,length(parameters)) 

#Nh_MoS = rep(NA,length(parameters)) 
#Nh_MoSvis = rep(NA,length(parameters)) 

#Nh_GNSUM = rep(NA,length(parameters))  


#Population and Survey
Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

survey = getSurvey(n_survey,Population)    # Survey
survey_hp = getSurvey(n_survey_hp,Population[Population$Hidden_Population==1,]) # Hiddden population survey


#Vector with the number of people in each subpopulation

v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population$Population == k) # N_k
  
}


# Study parameters
parameters = seq(from = 0, to = 0.50, length.out = 41)

################################################################################
# Estimation based on the different parameters

for (i in 1:length(parameters)) {
  # Parameter choice
  memory_factor = parameters[i]   
  
  for (j in 1:nrow(survey)) {
    vect_reach_re[j] = round(max(rnorm(1,mean = vect_reach[j], sd = memory_factor*vect_reach[j]),1))
  }
  survey$Reach_memory = vect_reach_re
  
  # Hidden population estimates
  Nh_real[i] = sum(Population$Hidden_Population) 
  
  Nh_basic_sum[i] = getNh_basic_sum(survey,N) 
  Nh_basicvis_sum[i] = getNh_basicvis_sum(survey,N,visibility_factor) 
  Nh_basic_mean[i] = getNh_basic_mean(survey,N) 
  Nh_basicvis_mean[i] = getNh_basicvis_mean(survey,N,visibility_factor)
  
  #Nh_PIMLE[i] = getNh_PIMLE(survey, v_pop_total, N)
  #Nh_PIMLEvis[i] = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
  
  #Nh_MLE[i] = getNh_MLE(survey, v_pop_total)
  #Nh_MLEvis[i] = getNh_MLEvis(survey, v_pop_total, visibility_factor)
  
  #Nh_MoS[i] = getNh_MoS(survey, v_pop_total, N)
  #Nh_MoSvis[i] = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
  
  #Nh_GNSUM[i] =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
}

################################################################################

################################################################################

# Graph 
x_1 = parameters
ggplot() + 
  geom_line(aes(x = x_1, y =  Nh_basic_sum, col = "Basic_sum")) + 
  geom_line(aes(x = x_1, y =  Nh_basicvis_sum, col = "Basic_vis_sum")) + 
  geom_line(aes(x = x_1, y =  Nh_basic_mean, col = "Basic_mean")) + 
  geom_line(aes(x = x_1, y =  Nh_basicvis_mean, col = "Basic_vis_mean")) + 
  #geom_line(aes(x = x_1, y =  Nh_PIMLEvis, col = "PIMLE_vis")) + 
  #geom_line(aes(x = x_1, y =  Nh_PIMLE, col = "PIMLE")) + 
  #geom_line(aes(x = x_1, y =  Nh_MLE, col = "MLE")) + 
  #geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "MLE_vis")) + 
  #geom_line(aes(x = x_1, y =  Nh_MoS, col = "MoS")) + 
  #geom_line(aes(x = x_1, y =  Nh_MoSvis, col = "MoS_vis")) + 
  #geom_line(aes(x = x_1, y =  Nh_GNSUM, col = "GNSUM")) +
  
  geom_line(aes(x = x_1, y =  Nh_real, col = "Real value")) +
  scale_color_discrete("Estimators") + 
  labs(title = "Prediction variability according to the memory factor",
       x = "Memory factor",
       y = "Hidden population estimate")
################################################################################

timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS ###########################
# Computation time (N=1000)  (my PC)
#timer ->  10.95956 secs
###########################################################################


