################################################################################################################
# Simulation based on the value of the memory factor of the Reach variable, leaving the rest of parameters fixed
################################################################################################################

t = Sys.time()


################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
memory_factor = 0      #Reach memory factor (parameter to change variance of the perturbations' normal)
################################################################################


N = 500                 # Population size
v_pop = c(1:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)   # Number of subpopulations
v_pop_prob = rep(1/length(v_pop), length(v_pop)) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 100            # Number of individuals we draw in the survey
n_survey_hp = 10          # Number of individuals we draw in the hidden population survey 

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

survey_hp = getSurvey(n_survey_hp,Population[Population$Hidden_Population==1,])


#Vector with the number of people in each subpopulation
v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}

# Study parameters
parameters = seq(from = 0, to = 1, length.out = 161)

#Dataframe to save the data
#simulaciones = data.frame(data = parameters)
#degrees = data.frame(data = parameters)



################################################################################
# AUXILIARY DATA FOR THE SIMULATION

vect_reach = Population$Reach
vect_reach_re =  rep(NA, nrow(Population))
b = 50 #Number of iterations for the simulation
lista_simulacion = list()
lista_degrees = list()

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
  
  memory_factor = parameters[i]   
  
  for (j in 1:nrow(Population)) {
    vect_reach_re[j] = round(max(rnorm(1,mean = vect_reach[j], sd = memory_factor*vect_reach[j]),1))
  }
  Population$Reach_memory = vect_reach_re
  
  #Variable reset
  Nh_real =  rep(NA,b) 
  
  Nh_basic_sum = rep(NA,b) 
  #Nh_basicvis_sum = rep(NA,b) 
  Nh_basic_mean = rep(NA,b) 
  #Nh_basicvis_mean = rep(NA,b)                                      
  
  #Nh_PIMLE = rep(NA,b) 
  #Nh_PIMLEvis = rep(NA,b) 
  
  #Nh_MLE = rep(NA,b) 
  #Nh_MLEvis = rep(NA,b) 
  
  #Nh_MoS = rep(NA,b) 
  #Nh_MoSvis = rep(NA,b) 
  
  #Nh_GNSUM = rep(NA,b) 
  
  lista_sim = list()
  lista_deg =list()
  
  
  #Iterations
  for (l in 1:b) {
    
    #We choose the same survey for each l in order to calculate the bias and variance
    #Surveys
    survey = Population[list_surveys[[l]],]
    survey_hp = Population[Population$Hidden_Population == 1,][list_surveys_hp[[l]],]
    
    #Visibility factor estimate
    vf_subpop = vf_subpop_es(survey_hp, Population, Mhp_vis, sub_memory_factor)
    
    #Hidden population estimates
    Nh_real = sum(Population$Hidden_Population) 
    
    Nh_basic_sum    = getNh_basic_sum(survey,N) 
    #Nh_basicvis_sum = getNh_basicvis_sum(survey,N,vf_subpop) 
    Nh_basic_mean    = getNh_basic_mean(survey,N) 
    #Nh_basicvis_mean = getNh_basicvis_mean(survey,N,vf_subpop) 
    
    #Nh_PIMLE    = getNh_PIMLE(survey, v_pop_total, N)
    #Nh_PIMLEvis = getNh_PIMLEvis(survey, v_pop_total, N, vf_subpop)
    
    #Nh_MLE     = getNh_MLE(survey, v_pop_total)
    #Nh_MLEvis  = getNh_MLEvis(survey, v_pop_total, vf_subpop)
    
    #Nh_MoS     = getNh_MoS(survey, v_pop_total, N)
    #Nh_MoSvis  = getNh_MoSvis(survey, v_pop_total, N, vf_subpop)
    
    #Nh_GNSUM   =  getNh_GNSUM(Population, survey, survey_hp, Mhp_vis, v_pop_total, N)
    
    
    #Dataframe for saving the estimates
    sim = data.frame(Nh_real = Nh_real)
    names(sim)[dim(sim)[2]] = str_c("Nh_real_",l)
    
    sim = cbind(sim,Nh_basic_sum = Nh_basic_sum)
    names(sim)[dim(sim)[2]] = str_c("Nh_basic_sum_",l)
    
    #sim = cbind(sim,Nh_basicvis_sum = Nh_basicvis_sum)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basicvis_sum_",l)
    
    sim = cbind(sim,Nh_basic_mean = Nh_basic_mean)
    names(sim)[dim(sim)[2]] = str_c("Nh_basic_mean_",l)
    
    #sim = cbind(sim,Nh_basicvis_mean = Nh_basicvis_mean)
    #names(sim)[dim(sim)[2]] = str_c("Nh_basicvis_mean_",l)
    
    #sim = cbind(sim,Nh_PIMLE = Nh_PIMLE)
    #names(sim)[dim(sim)[2]] = str_c("Nh_PIMLE_",l)
    
    #sim = cbind(sim,Nh_PIMLEvis = Nh_PIMLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_PIMLEvis_",l)
    
    #sim = cbind(sim,Nh_MLE = Nh_MLE)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MLE_",l)
    
    #sim = cbind(sim,Nh_MLEvis = Nh_MLEvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MLEvis_",l)
    
    #sim = cbind(sim,Nh_MoS = Nh_MoS)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MoS_",l)
    
    #sim = cbind(sim,Nh_MoSvis = Nh_MoSvis)
    #names(sim)[dim(sim)[2]] = str_c("Nh_MoSvis_",l)
    
    #sim = cbind(sim,Nh_GNSUM = Nh_GNSUM)
    #names(sim)[dim(sim)[2]] = str_c("Nh_GNSUM_",l)
    
    lista_sim[[l]] = sim
    
    # Degrees
    
    # Real reach
    d_real = survey$Reach
    
    # Reach applying the memory fector
    d_basic = survey$Reach_memory
    
    # Estimated Reach by MLE and PIMLE
    d_iest = c()
    for (j in 1:nrow(survey)) {
      d_iest[j] = N * sum(survey[j,tail(names(survey),length(v_pop_total))])/sum(v_pop_total)
    }
    d_MLE = d_iest 
    
    # Estimated Reach by MoS
    d_i_est = rep(NA, nrow(survey))
    for (j in 1:nrow(survey)) {
      d_i_est[j] = (sum((survey[j,tail(names(survey),length(v_pop_total))])/v_pop_total))/length(v_pop_total) * N
    }
    d_MoS=d_i_est
    
    deg = matrix(c(d_real,d_basic,d_MLE,d_MoS),nrow = 1)
    deg = as.data.frame(deg)
    for (ind in 1:length(d_real)) {
      colnames(deg)[ind]=str_c("d_real_",ind,"_",l)
      colnames(deg)[ind+length(d_real)]=str_c("d_basic_",ind,"_",l)
      colnames(deg)[ind+2*length(d_real)]=str_c("d_MLE_",ind,"_",l)
      colnames(deg)[ind+3*length(d_real)]=str_c("d_basic_",ind,"_",l)
    }
    lista_deg[[l]]=deg
  }
  simulacion = bind_cols(lista_sim)
  lista_simulacion[[i]] = simulacion
  print(i)
  
  degr = bind_cols(lista_deg)
  lista_degrees[[i]] =degr
}

simulaciones = bind_rows(lista_simulacion)
simulaciones = cbind(simulaciones, data = parameters)


degrees=bind_rows(lista_degrees)
degrees =cbind(degrees,data=parameters)

################################################################################
write.csv(simulaciones,                           # Data frame 
          file = "Simulations_memoryfactor_notdisjoint",   # Csv name
          row.names = TRUE )                      # Row names: TRUE or FALSE 
################################################################################




timer = Sys.time() - t
timer

#################### COMPUTATION TIME ANALYSIS ###########################

# Computation time (N=1000) (my PC)
#timer -> 29.46286 secs not saving all the unnecessary estimators 

# Computation time (N=10000) (office PC) (length(parameters) = 41)
#timer -> 6.755314 mins not saving all the unnecessary estimators 

# Computation time (N=10000) (office PC) (length(parameters) = 161)
#timer -> 8.205191 mins not saving all the unnecessary estimators


#Problem: MoS and PIMLE have computation time of 0.2 per iteration
# 0.2 * 25 * 41 = 200 sec -> 3,33 min
# 3,33 * 4 = 13,3 min

###########################################################################




