##############################################################################################
# Simulation based on the value of the visibility factor, leaving the rest of parameters fixed
##############################################################################################


################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
visibility_factor = 1  
################################################################################


N = 1000                  # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0     # Visibility factor (Binomial's probability)
seed = 207                # Seed
set.seed(seed)

#Graph
dim = 1    # Graph dimension 
nei = 75   # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


# Study parameters
parameters = seq(from = 0.1, to = 1, length.out = 19)


#Dataframe to save the data
simulaciones = data.frame()


#Population and Survey
Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population

########################### WARNING #######################################
#We use Mhp not Mhp_vis as a model in order to change it in each iteration
Mhp = Graph_population_matrix[[3]]         # Population's visibility matrix
###########################################################################

survey_hp = getSurvey(n_survey_hp, Population[Population$Hidden_Population==1,])



#Vector with the number of people in each subpopulation
v_pop_total = rep(NA, n_pop)
for (j in 1:n_pop) {
  v_pop_total[j] = sum(Population$Population == j)
}


################################################################################

# AUXILIARY DATA FOR THE SIMULATION

vect_hp_vis = rep(NA,nrow(Population))
# The data is fixed to use it as a reference in every loop
vect_hp = Population$HP_total
lista_simulacion = list()

b = 25 #Number of iterations for the simulation

# Surveys representing the different iterations. 
# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population), n_survey, replace = FALSE)
}

# Simulation

for (i in 1:length(parameters)) {
  
  visibility_factor = parameters[i] #Study parameter
  Mhp_vis =  apply(Mhp,c(1,2), berHP, p = visibility_factor) #New visibility applied
  
  # Hidden population known by k after applying the factor
  k = 1
  for (j in as.numeric(rownames(Population))) {
    vect_hp_vis[k] = sum(Mhp_vis[j,])
    k = k +1
  }
  
  Population$HP_total_apvis = vect_hp_vis
  
  Nh_real =  rep(NA,b) 
  
  Nh_basic = rep(NA,b) 
  Nh_basicvis = rep(NA,b)                                      
  
  Nh_PIMLE = rep(NA,b) 
  Nh_PIMLEvis = rep(NA,b) 
  
  Nh_MLE = rep(NA,b) 
  Nh_MLEvis = rep(NA,b) 
  
  Nh_MoS = rep(NA,b) 
  Nh_MoSvis = rep(NA,b) 

  Nh_GNSUM = rep(NA,b) 
  
  lista_sim = list()
  
  for (l in 1:b) {
      
      #We choose the same survey for each l in order to calculate the bias and variance
      survey  = Population[list_surveys[[l]],]
    
      Nh_real = sum(Population$Hidden_Population) 
    
      Nh_basic    = getNh_basic(survey,N) 
      Nh_basicvis = getNh_basicvis(survey,N,visibility_factor) 
    
      Nh_PIMLE    = getNh_PIMLE(survey, v_pop_total, N)
      Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
    
      Nh_MLE     = getNh_MLE(survey, v_pop_total)
      Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, visibility_factor)
    
      Nh_MoS     = getNh_MoS(survey, v_pop_total, N)
      Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
    
      Nh_GNSUM   =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
      
      sim = data_frame(Nh_real = Nh_real)
      names(sim)[dim(sim)[2]] = str_c("Nh_real_",l)
      
      sim = cbind(sim,Nh_basic = Nh_basic)
      names(sim)[dim(sim)[2]] = str_c("Nh_basic",l)
      
      sim = cbind(sim,Nh_basicvis = Nh_basicvis)
      names(sim)[dim(sim)[2]] = str_c("Nh_basicvis_",l)
      
      sim = cbind(sim,Nh_PIMLE = Nh_PIMLE)
      names(sim)[dim(sim)[2]] = str_c("Nh_PIMLE_",l)
      
      sim = cbind(sim,Nh_PIMLEvis = Nh_PIMLEvis)
      names(sim)[dim(sim)[2]] = str_c("Nh_PIMLEvis_",l)
      
      sim = cbind(sim,Nh_MLE = Nh_MLE)
      names(sim)[dim(sim)[2]] = str_c("Nh_MLE_",l)
      
      sim = cbind(sim,Nh_MLEvis = Nh_MLEvis)
      names(sim)[dim(sim)[2]] = str_c("Nh_MLEvis_",l)
      
      sim = cbind(sim,Nh_MoS = Nh_MoS)
      names(sim)[dim(sim)[2]] = str_c("Nh_MoS_",l)
      
      sim = cbind(sim,Nh_MoSvis = Nh_MoSvis)
      names(sim)[dim(sim)[2]] = str_c("Nh_MoSvis_",l)
      
      sim = cbind(sim,Nh_GNSUM = Nh_GNSUM)
      names(sim)[dim(sim)[2]] = str_c("Nh_GNSUM_",l)
      
      lista_sim[[l]] = sim
      print(l)
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  

}

simulaciones = bind_rows(lista_simulacion)


################################################################################
simulaciones
write.csv(simulaciones,                          # Data frame 
          file = "Simulations_visibilityfactor", # Csv name
          row.names = TRUE )                     # Row names: TRUE o FALSE 
################################################################################

#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC)
#timer -> 16.54465 mins  

# Computation time (N=10000) (office PC)
#timer -> 

###########################################################################


