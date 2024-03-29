#########################################
# Visibility factor estimate simulation #
#########################################

t = Sys.time()
################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
sub_memory_factor = 0      #Subpopulation memory factor (parameter to change variance of the perturbations' normal)
################################################################################

N = 10000                 # Population size
v_pop = c(1:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)   # Number of subpopulations
v_pop_prob = rep(1/length(v_pop), length(v_pop)) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
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



#Vector with the number of people in each subpopulation

v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}

# Study parameters
parameters = seq(from = 0, to = 0.5, length.out = 89)

################################################################################

# AUXILIARY DATA FOR THE SIMULATION

vis_pob_reset = Population[,(ncol(Population)-(n_pop-1)):ncol(Population)]
b = 50 #Number of iterations for the simulation
lista_simulacion = list()

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


#Simulation  

for (i in 1:length(parameters)) {
  
  sub_memory_factor = parameters[i]   
  
  for(j in 1:n_pop){
    v_1 = rep(NA,nrow(Population))
    vis_pob = Population[,ncol(Population)-(n_pop-j)]
    for(k in 1:nrow(Population)) {
      v_1[k] = max(round(rnorm(1, mean = vis_pob[k], sd = sub_memory_factor*vis_pob[k])),0)
    }
    Population[,ncol(Population)-(n_pop-j)] = v_1
  }
  
  lista_sim = list() 
  
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],] #Survey
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],] #Hidden population survey
    
    #Visibility factor estimate
    vf_subpop = vf_subpop_es(survey_hp, Population, Mhp_vis, sub_memory_factor)
    vf_subpop_out =vf_subpop_es_out(survey_hp, Population, Mhp_vis, sub_memory_factor)
    
    #Dataframe for saving the estimates
    sim = data.frame(vf_subpop = vf_subpop)
    names(sim)[dim(sim)[2]] = str_c("vf_subpop",l)
    
    sim = cbind(sim,vf_subpop_out = vf_subpop_out)
    names(sim)[dim(sim)[2]] = str_c("vf_subpop_out",l)
    
    lista_sim[[l]] = sim
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  Population[,(ncol(Population)-(n_pop-1)):ncol(Population)] = vis_pob_reset
  print(i)
  
}

simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = parameters)



################################################################################
write.csv(simulaciones,                                  # Data frame
          file = "Simulation_visibilityfactorestimate_subpop_notdisjoint", # Csv's name
          row.names = TRUE )                             # Row names: TRUE o FALSE 
################################################################################



timer = Sys.time() - t
timer
#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (virtual machine)
#timer ->  44.03315 secs mins   

# Computation time (N=10000) (virtual machine) IN PROCESS
#timer ->  9.609865 mins

###########################################################################