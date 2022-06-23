###################################################################################
# Graph based on the number of subpopulations, leaving the rest of parameters fixed
###################################################################################


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
parameters = round(seq(from = 2, to = 20, length.out = 10))


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




#Population and Survey

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

survey = getSurvey(n_survey,Population)
survey_hp = getSurvey(n_survey_hp,Population[Population$Hidden_Population==1,])

# Auxiliar data for the simulation
k = length(v_pop)

for (i in 1:length(parameters)) {
  m_pop = parameters[i]
  n_colum = ncol(Population)
  v_pop = c(0:m_pop)
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
  
  #Surveys
  survey = getSurvey(n_survey,Population)
  survey_hp = getSurvey(n_survey_hp, Population[Population$Hidden_Population==1,])
  
  
  #Hidden population estimate
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
  print(ncol(Population))
}


x_1 = parameters
ggplot() + 
  geom_line(aes(x = x_1, y =  Nh_basic, col = "Basic")) + 
  geom_line(aes(x = x_1, y =  Nh_basicvis, col = "Basic_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLEvis, col = "PIMLE_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLE, col = "PIMLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLE, col = "MLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "MLE_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_MoS, col = "MoS")) + 
  geom_line(aes(x = x_1, y =  Nh_MoSvis, col = "MoS_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_GNSUM, col = "GNSUM")) +
  
  geom_line(aes(x = x_1, y =  Nh_real, col = "Real value")) +
  scale_color_discrete("Estimators") + 
  labs(title = "Prediction variability according to the subpopulation's number",
       x = "Subpopulation's number",
       y = "Hidden population estimate")


