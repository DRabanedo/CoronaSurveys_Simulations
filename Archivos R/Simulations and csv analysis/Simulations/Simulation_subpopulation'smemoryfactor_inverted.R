###########################################################################################################
# Graph based on the value of the memory factor of the subpopulations, leaving the rest of parameters fixed
###########################################################################################################

t = Sys.time()
################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
sub_memory_factor = 0      #Subpopulation memory factor (parameter to change variance of the perturbations' normal)
################################################################################

N = 1000                  # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probability of each subpopulation
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
  v_pop_total[k] = sum(Population$Population == k) # N_k
  
}


# Study parameters
parameters = seq(from = 0, to = 1, length.out = 21)

#Dataframe to save the data
simulaciones = data.frame(datos =  parameters)

################################################################################

# AUXILIARY DATA FOR THE SIMULATION

vis_pob_reset = Population[,(ncol(Population)-(n_pop)):ncol(Population)]
b = 25 #Number of iterations for the simulation
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


  
  
  
for (i in 1:length(parameters)) {
    
  sub_memory_factor = parameters[i]   
    
  for(j in 0:n_pop){
    v_1 = rep(NA,nrow(Population))
    vis_pob = Population[,ncol(Population)-(n_pop-j)]
    for(k in 1:nrow(Population)) {
      v_1[k] = max(round(rnorm(1, mean = vis_pob[k], sd = sub_memory_factor*vis_pob[k])),0)
    }
    Population[,ncol(Population)-(n_pop-j)] = v_1
  }
  
  Nh_real =  rep(NA,b) 
  
  #Nh_basic = rep(NA,b) 
  #Nh_basicvis = rep(NA,b)                                      
  
  Nh_PIMLE = rep(NA,b) 
  #Nh_PIMLEvis = rep(NA,b) 
  
  Nh_MLE = rep(NA,b) 
  #Nh_MLEvis = rep(NA,b) 
  
  Nh_MoS = rep(NA,b) 
  #Nh_MoSvis = rep(NA,b) 
  
  #Nh_GNSUM = rep(NA,b) 
  
  lista_sim = list() 
  
  for (l in 1:b) {
    #We choose the same survey for each l in order to calculate the bias and variance
    survey  = Population[list_surveys[[l]],]
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[list_surveys_hp[[l]],]
    
    Nh_real = sum(Population$Hidden_Population) 
    
    #Nh_basic    = getNh_basic(survey,N) 
    #Nh_basicvis = getNh_basicvis(survey,N,visibility_factor) 
    
    Nh_PIMLE    = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, visibility_factor)
    
    Nh_MLE     = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, visibility_factor)
    
    Nh_MoS     = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, visibility_factor)
    
    #Nh_GNSUM   =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
    
    sim = data_frame(Nh_real = Nh_real)
    names(sim)[dim(sim)[2]] = str_c("Nh_real_",l)
    
    #sim = cbind(sim,Nh_basic = Nh_basic)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basic",l)
    
    #sim = cbind(sim,Nh_basicvis = Nh_basicvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basicvis_",l)
    
    sim = cbind(sim,Nh_PIMLE = Nh_PIMLE)
    names(sim)[dim(sim)[2]] = str_c("Nh_PIMLE_",l)
    
    #sim = cbind(sim,Nh_PIMLEvis = Nh_PIMLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_PIMLEvis_",l)
    
    sim = cbind(sim,Nh_MLE = Nh_MLE)
    names(sim)[dim(sim)[2]] = str_c("Nh_MLE_",l)
    
    #sim = cbind(sim,Nh_MLEvis = Nh_MLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MLEvis_",l)
    
    sim = cbind(sim,Nh_MoS = Nh_MoS)
    names(sim)[dim(sim)[2]] = str_c("Nh_MoS_",l)
    
    #sim = cbind(sim,Nh_MoSvis = Nh_MoSvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MoSvis_",l)
    
    #sim = cbind(sim,Nh_GNSUM = Nh_GNSUM)
    #names(sim)[dim(sim)[2]] = str_c("Nh_GNSUM_",l)
    
    lista_sim[[l]] = sim
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  Population[,(ncol(Population)-(n_pop)):ncol(Population)] = vis_pob_reset
  print(i)
  
}
  
simulaciones = bind_rows(lista_simulacion)
  



################################################################################
simulaciones
write.csv(simulaciones,                        # Data frame
          file = "Simulaciones_subpopulation'smemoryfactor", # Csv's name
          row.names = TRUE )                   # Row names: TRUE o FALSE 
################################################################################

timer = Sys.time() - t
timer
#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC)
#timer -> 1.044038 mins    

# Computation time (N=10000) (office PC)
#timer ->  

###########################################################################