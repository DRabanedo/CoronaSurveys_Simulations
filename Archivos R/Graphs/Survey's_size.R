#############################################################################
# Graph based on the size of the survey, leaving the rest of parameters fixed
#############################################################################


N = 1000                 # Population size
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
parameters = round(seq(from = 1, to = 100, length.out = 20))


#Population and Survey

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

survey_hp = getSurvey(n_survey_hp,Population[Population$Hidden_Population==1,])

#Vector with the number of people in each subpopulation

v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population$Population == k) # N_k
  
}


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

Nh_Direct = rep(NA,length(parameters))  


for (i in 1:length(parameters)) {
  
  n_survey = parameters[i]
  survey = getSurvey(n_survey,Population)
  
  
  
  Nh_real[i] = sum(Population$Hidden_Population) 
  
  Nh_basic[i] = getNh_basic(survey,N) 
  Nh_basicvis[i] = getNh_basicvis(survey,N,visibility_factor) 
  
  Nh_PIMLE[i] = getNh_PIMLE(survey, v_pop_total, N)
  Nh_PIMLEvis[i] = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
  
  Nh_MLE[i] = getNh_MLE(survey, v_pop_total)
  Nh_MLEvis[i] = getNh_MLEvis(survey, v_pop_total, visibility_factor)
  
  Nh_MoS[i] = getNh_MoS(survey, v_pop_total, N)
  Nh_MoSvis[i] = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
  
  Nh_GNSUM[i] =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
  
  Nh_Direct[i] = getNh_Direct(survey,N)
}


################################################################################



x_1 = parameters
ggplot() + 
  geom_line(aes(x = x_1, y =  Nh_basic, col = "Basic")) + 
  #geom_line(aes(x = x_1, y =  Nh_basicvis, col = "Basic_vis")) + 
  #geom_line(aes(x = x_1, y =  Nh_PIMLEvis, col = "PIMLE_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLE, col = "PIMLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLE, col = "MLE")) + 
  #geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "MLE_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_MoS, col = "MoS")) + 
  #geom_line(aes(x = x_1, y =  Nh_MoSvis, col = "MoS_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_GNSUM, col = "GNSUM")) +
  geom_line(aes(x = parameters, y =  Nh_Direct, col = "Direct")) +
  geom_line(aes(x = x_1, y =  Nh_real, col = "Real value")) +
  scale_color_discrete("Estimators") + 
  labs(title = "Prediction variability according to the survey's size",
       x = "Survey's size",
       y = "Hidden population estimate")




