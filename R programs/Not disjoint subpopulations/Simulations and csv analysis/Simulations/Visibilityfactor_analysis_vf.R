##############################################################################################
# Simulation based on the value of the visibility factor, leaving the rest of parameters fixed
##############################################################################################

t = Sys.time()
################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
visibility_factor = 1  
################################################################################


N = 10000                 # Population size
v_pop = c(1:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)     # Number of subpopulations
v_pop_prob = rep(1/length(v_pop), length(v_pop)) # Probability of each subpopulation
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
parameters = seq(from = 0.1, to = 1, length.out = 89)


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
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}


################################################################################

# AUXILIARY DATA FOR THE SIMULATION

vect_hp_vis = rep(NA,nrow(Population))
# The data is fixed to use it as a reference in every loop
vect_hp = Population$HP_total
lista_simulacion = list()

b = 50 #Number of iterations for the simulation

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
  
  # Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],] #Survey
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],] #Hidden population survey
    
    #Visibility factor estimate
    vf_subpop = vf_subpop_es(survey_hp, Population, Mhp_vis, sub_memory_factor)
    vf_subpop_out =vf_subpop_es_out(survey_hp, Population, Mhp_vis, sub_memory_factor)
    
    vf_reach = vf_reach_es(survey_hp, Population, Mhp_vis, memory_factor)
    vf_reach_out = vf_reach_es_out(survey_hp, Population, Mhp_vis, memory_factor)
    
    
    
    #Dataframe for saving the estimates
    sim = data.frame(vf_subpop = vf_subpop)
    names(sim)[dim(sim)[2]] = str_c("vf_subpop",l)
    
    sim = cbind(vf_subpop_out = vf_subpop_out)
    names(sim)[dim(sim)[2]] = str_c("vf_subpop_out",l)
    
    sim = cbind(vf_reach = vf_reach)
    names(sim)[dim(sim)[2]] = str_c("vf_reach_",l)
    
    sim = cbind(vf_reach_out = vf_reach_out)
    names(sim)[dim(sim)[2]] = str_c("vf_reach_out_",l)
    
    lista_sim[[l]] = sim
    print(l)
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  
  
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
          file = "Simulations_visibilityfactorestimate_vf_notdisjoint", # Csv name
          row.names = TRUE )                     # Row names: TRUE o FALSE 
################################################################################



# Computation time (N=1000) (office PC)
#timer -> 12.12937 mins

# Computation time (N=10000) (office PC)
#timer ->  2.061217 hours

#Computation time (N = 10000) (office PC) (89)
#timer -> 6.721095 hours
###########################################################################

