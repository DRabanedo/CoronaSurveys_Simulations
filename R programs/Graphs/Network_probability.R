##########################################################################
# Graph based on the graph structure, leaving the rest of parameters fixed
##########################################################################

t = Sys.time()

N = 1000                 # Population size
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
parameters = seq(from = 0.05, to = 1, length.out = 21)


# Variable reset
Nh_real =  rep(NA,length(parameters)) 

Nh_basic_sum = rep(NA,length(parameters)) 
Nh_basic_mean = rep(NA,length(parameters)) 
#Nh_basicvis_sum = rep(NA,length(parameters)) 
#Nh_basicvis_mean = rep(NA,length(parameters)) 

Nh_PIMLE = rep(NA,length(parameters)) 
#Nh_PIMLEvis = rep(NA,length(parameters)) 

Nh_MLE = rep(NA,length(parameters)) 
#Nh_MLEvis = rep(NA,length(parameters)) 

Nh_MoS = rep(NA,length(parameters)) 
#Nh_MoSvis = rep(NA,length(parameters)) 

Nh_GNSUM = rep(NA,length(parameters))  

#Aux data for the simulation
Population_ref = genPopulation(N, v_pop, v_pop_prob,hp_prob)
survey_hp_ref = sample(nrow(Population[Population$Hidden_Population==1,]), n_survey_hp, replace = FALSE)
survey_ref = sample(nrow(Population), n_survey, replace = FALSE)

################################################################################
# Estimation based on the different parameters


for (i in 1:length(parameters)) {
  # Parameter Choice
  p = parameters[i]
  
  #Population
  Population = Population_ref
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_pop = length(v_pop)-1
  
  # Initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # Reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP, p = visibility_factor)
  
  
  for (h in 1:N) {
    # net_sw[[h]], list with one element, the list of the adjacent vertices to h
    
    vect_hp[h] = sum(Mhp[h,])
    vect_reach[h] = length(net_sw[[h]][[1]])
    vect_hp_vis[h] = sum(Mhp_vis[h,])
    
    vect_reach_re[h] = round(rnorm(1, mean = vect_reach[h], sd = memory_factor*vect_reach[h]))
  }
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_memory = vect_reach_re)
  Population = cbind(Population, HP_total = vect_hp) 
  Population = cbind(Population, HP_total_apvis = vect_hp_vis)
  
  for(j in 0:n_pop){
    v_1 = rep(NA,N)
    for(f in 1:N) {
      vis_pob = sum(Population[net_sw[[f]][[1]],]$Population == j)
      # Visibility of population j by f, applying a normal in order to represent the real visibility
      v_1[f] = round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob))
    }
    
    Population = cbind(Population,SubPopulation_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
  
  #Surveys
  survey = getSurvey(n_survey,Population)    # Population survey
  survey_hp = getSurvey(n_survey_hp,Population[Population$Hidden_Population==1,]) # Hiddden population survey
  
  #Vector with the number of people in each subpopulation
  
  v_pop_total = rep(NA, n_pop)
  for (k in 1:n_pop) {
    v_pop_total[k] = sum(Population$Population == k) # N_k
    
  }
  
  # Hidden population estimates
  Nh_real[i] = sum(Population$Hidden_Population) 
  
  Nh_basic_sum[i] = getNh_basic_sum(survey,N) 
  #Nh_basicvis_sum[i] = getNh_basicvis_sum(survey,N,visibility_factor) 
  Nh_basic_mean[i] = getNh_basic_mean(survey,N) 
  #Nh_basicvis_mean[i] = getNh_basicvis_mean(survey,N,visibility_factor) 
  
  Nh_PIMLE[i] = getNh_PIMLE(survey, v_pop_total, N)
  #Nh_PIMLEvis[i] = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
  
  Nh_MLE[i] = getNh_MLE(survey, v_pop_total)
  #Nh_MLEvis[i] = getNh_MLEvis(survey, v_pop_total, visibility_factor)
  
  Nh_MoS[i] = getNh_MoS(survey, v_pop_total, N)
  #Nh_MoSvis[i] = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
  
  Nh_GNSUM[i] =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
}
################################################################################


################################################################################
# Graph 
x_1 = parameters
ggplot() + 
  geom_line(aes(x = x_1, y =  Nh_basic_sum, col = "Basic_sum")) + 
  #geom_line(aes(x = x_1, y =  Nh_basicvis_sum, col = "Basic_vis_sum")) + 
  geom_line(aes(x = x_1, y =  Nh_basic_mean, col = "Basic_mean")) + 
  #geom_line(aes(x = x_1, y =  Nh_basicvis_mean, col = "Basic_vis_mean")) + 
  #geom_line(aes(x = x_1, y =  Nh_PIMLEvis, col = "PIMLE_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLE, col = "PIMLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLE, col = "MLE")) + 
  #geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "MLE_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_MoS, col = "MoS")) + 
  #geom_line(aes(x = x_1, y =  Nh_MoSvis, col = "MoS_vis")) + 
  geom_line(aes(x = x_1, y =  Nh_GNSUM, col = "GNSUM")) +
  geom_line(aes(x = x_1, y =  Nh_real, col = "Real value")) +
  scale_color_discrete("Estimators") + 
  labs(title = "Prediction variability according to network structure",
       x = "Aleatorization of connections",
       y = "Hidden population estimate")
################################################################################

Sys.time() - t

#################### COMPUTATION TIME ANALYSIS ###########################
# Computation time (N=1000)  (my PC)
#timer ->  3.394592 mins 
###########################################################################


