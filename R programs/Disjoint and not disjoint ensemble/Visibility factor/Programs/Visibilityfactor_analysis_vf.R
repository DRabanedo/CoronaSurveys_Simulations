##############################################################################################
# Simulation based on the value of the visibility factor, leaving the rest of parameters fixed
##############################################################################################

t = Sys.time()
################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
visibility_factor = 1  
################################################################################


N = 10000                     # Population size
v_pop = c(1:5)                # Subpopulations vector 
n_pop = length(v_pop)         # Number of subpopulations
v_pop_prob = rep(1/10, 5)     # Probability of each subpopulation. As we are working with disjoint and no disjoint subpopulations
                              # sum(v_pop_prob) < 1. 
hp_prob = 0.1                 # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500                # Number of individuals we draw in the survey
n_survey_hp = 50              # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0         # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0             
seed = 207                    # Seed
set.seed(seed)

#Graph
dim = 1    # Graph dimension 
nei = 75   # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


# Study parameters
parameters = seq(from = 0.1, to = 0.95, length.out = 50)


#Population and Survey
Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population

########################### WARNING #######################################
#We use Mhp not Mhp_vis as a model in order to change it in each iteration
Mhp = Graph_population_matrix[[3]]         # Population's visibility matrix
###########################################################################



#Vector with the number of people in each subpopulation
v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}


################################################################################

# AUXILIARY DATA FOR THE SIMULATION

vect_hp_vis = rep(NA,nrow(Population))
# The data is fixed to use it as a reference in every loop
vect_hp = Population$HP_total
lista_simulacion = list()

b = 100 #Number of iterations for the simulation

# Surveys representing the different iterations. 
# The surveys are fixed so the variance and bias can be calculated.

list_surveys = list()
for (h in 1:b) {
  list_surveys[[h]] = sample(nrow(Population), n_survey, replace = FALSE)
}

list_surveys_hp = list()
for (h in 1:b) {
  list_surveys_hp[[h]] = sample(nrow(Population[Population$Hidden_Population == 1,]), n_survey_hp, replace = FALSE)
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
  
  lista_sim = list()
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    
    # Binomial 
    vf_bin_est = vf_bin(survey_hp,Population, Mhp_vis, memory_factor)
    
    # Double truncation
    vf_dt_bin_est = vf_dt_bin(survey_hp,Population, Mhp_vis, memory_factor)
    vf_dt_fbin_est = vf_dt_fbin(survey_hp,Population, Mhp_vis, memory_factor) 
    vf_dt_dt_est = vf_dt_dt(survey_hp,Population, Mhp_vis, memory_factor) 
    vf_dt_ut_est = vf_dt_ut(survey_hp,Population, Mhp_vis, memory_factor) 
    
    
    
    # Unilateral truncation 
    vf_ut_bin_est = vf_ut_bin(survey_hp,Population, Mhp_vis, memory_factor) 
    vf_ut_fbin_est = vf_ut_fbin(survey_hp,Population, Mhp_vis, memory_factor)
    vf_ut_dt_est = vf_ut_dt(survey_hp,Population, Mhp_vis, memory_factor)
    vf_ut_ut_est = vf_ut_ut(survey_hp,Population, Mhp_vis, memory_factor)
    
    
    
    #Dataframe for saving the estimate
    
    sim = data.frame(vf_dt_bin = vf_dt_bin_est)
    names(sim)[dim(sim)[2]] = str_c("vf_dt_bin_",l)
    
    sim = cbind(sim,vf_bin = vf_bin_est)
    names(sim)[dim(sim)[2]] = str_c("vf_bin_",l)
    
    sim = cbind(sim,vf_dt_fbin = vf_dt_fbin_est)
    names(sim)[dim(sim)[2]] = str_c("vf_dt_fbin_",l)
    
    sim = cbind(sim,vf_dt_dt = vf_dt_dt_est)
    names(sim)[dim(sim)[2]] = str_c("vf_dt_dt_",l)
    
    sim = cbind(sim,vf_dt_ut = vf_dt_ut_est)
    names(sim)[dim(sim)[2]] = str_c("vf_dt_ut_",l)
    
    sim = cbind(sim,vf_ut_bin = vf_ut_bin_est)
    names(sim)[dim(sim)[2]] = str_c("vf_ut_bin_",l)
    
    sim = cbind(sim,vf_ut_fbin = vf_ut_fbin_est)
    names(sim)[dim(sim)[2]] = str_c("vf_ut_fbin_",l)
    
    sim = cbind(sim,vf_ut_dt = vf_ut_dt_est)
    names(sim)[dim(sim)[2]] = str_c("vf_ut_dt_",l)
    
    sim = cbind(sim,vf_ut_ut = vf_ut_ut_est)
    names(sim)[dim(sim)[2]] = str_c("vf_ut_ut_",l)
    
    lista_sim[[l]] = sim
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  print(i)
  
}

simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = parameters)

timer = Sys.time() - t
timer
#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC)
#timer -> 16.54465 mins  

################################################################################
write.csv(simulaciones,                          # Data frame 
          file = "Simulations_visibilityfactorestimate_vf_207.csv", # Csv name
          row.names = TRUE )                     # Row names: TRUE o FALSE 
################################################################################


################################################################################

# Computation time (N=1000) (office PC)
#timer -> 12.12937 mins

# Computation time (N=10000) (office PC)
#timer ->  2.061217 hours

################################################################################

