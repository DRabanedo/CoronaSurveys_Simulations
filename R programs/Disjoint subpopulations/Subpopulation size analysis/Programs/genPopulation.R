getPopulation = function(N, dis_populations,prob_vector,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor, sub_memory_factor){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  #N: population size
  #dis_populations: vector with the populations
  #prob_vector: vector with the population's probabilities
  #PropHiddenPop: Hidden Population proportion
  #dim: Integer constant, the dimension of the starting lattice.
  #nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p: Real constant between zero and one, the rewiring probability.
  #visibility_factor: the visibility factor
  #memory factor: numeric value, Reach divided by the standard deviation of the normal we use to correct the Reach
  # sub_memory_factor: the subpopulations visibility divided by the standard deviation of the normal distributions we use to correct the subpopulations visibility
  #                    it is applied to each subpopulation
  # sub_memory_factor: the subpopulations' visibility divided by the standard deviation of the normals we use to correct the subpopulations' visibility
  
  Population = genPopulation(N, dis_populations, prob_vector,PropHiddenPop)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(dis_populations)-1
  # initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # Reach vector applying memory error
  
  for (i in 1:N) {
    # net_sw[[i]], list with one element, the list of the adjacent vertices to i
    
    vect_hp[i] = sum(Population[net_sw[[i]][[1]],]$Hidden_Population)
    vect_reach[i] = length(net_sw[[i]][[1]])
    vect_hp_vis[i] = rbinom(1,vect_hp[i],visibility_factor)
    
    vect_reach_re[i] = max(1,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
  }
  
  
  Population = cbind(Population, Reach = vect_reach)
  Population = cbind(Population, Reach_memory = vect_reach_re)
  Population = cbind(Population, HP_total = vect_hp) 
  Population = cbind(Population, HP_total_apvis = vect_hp_vis)
  
  for(j in 0:n_populations){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Population[net_sw[[i]][[1]],]$Population == j)
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("KP_total_apvis_",j)
  }
  
  
  return(Population)
}


################################################################################################################
# Simulation based on the value of the memory factor of the Reach variable, leaving the rest of parameters fixed
################################################################################################################

t = Sys.time()


################ WARNING #########################################################
# IT IS IMPORTANT TO LEAVE THIS FACTOR AS 0, SINCE THE FIRST POPULATION THAT WE ARE 
# GOING TO BUILD WILL BE A BASIS FOR THE REST
memory_factor = 0      #Reach memory factor (parameter to change variance of the perturbations' normal)
################################################################################


N = 1000                 # Population size
v_pop = c(0:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1   # Number of subpopulations
v_pop_prob = rep(1/length(v_pop), length(v_pop)) #Probability of each subpopulation
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



#Population and Survey
Graph_population_matrix = getPopulation(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)
