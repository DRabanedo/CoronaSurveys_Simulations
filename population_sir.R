# Population using pandemic model HP distribution #

# SIR Hidden Population generation #

gen_SIRpop = function(N, net ,beta, gamma, chosen_nodes, infected_people){
  # SIR method for determination of the hidden population
  
  # beta: infection rate
  # gamma: removal rate
  # infected_people: number of infected people
  # chosen_nodes: nodes from which we start
  
  # Loop variables
  infected_nodes = c()
  infected_recovered = c()
  
  infected_nodes = c(infected_nodes,chosen_nodes)
  
  while (length(infected_nodes) < infected_people) {
    
    #Infection stage
    for (node in infected_nodes) {
      for (neighbor in net[[node]][[1]]){
        if (runif(1) < beta & !(neighbor %in% infected_nodes) & !(neighbor %in% infected_recovered )) {
          infected_nodes = c(infected_nodes,neighbor)
          #print(infected_nodes)
        }
      }
    }
    
    # Removal stage
    infected_survivors = c()
    
    for (node in infected_nodes) {
      if (runif(1) < gamma){
        infected_recovered = c(infected_recovered,node)
      }
      else{
        infected_survivors = c(infected_survivors,node)
      }
    }
    infected_nodes = infected_survivors
  }
  
  final_infected_nodes = infected_nodes[1:infected_people]
  
  hp_vector = rep(NA, N)
  
  for (i in 1:N){
    if (as.logical(sum(i %in% final_infected_nodes))){
      hp_vector[i] = 1
    }
    else{
      hp_vector[i] = 0
    }
  }
  return(hp_vector)
}


# Population function #

getData = function(N, prob_vect, prob_hp, dim, nei, p, vis_factor, mem_factor, sub_mem_factor, beta, gamma, chosen_nodes, infected_people){
  # list, contains the network, the population data and the matrix for the GNSUM
  
  # N:  Population size
  # prob_vect:  Vector with the population's probabilities
  # prob_hp: Hidden Population proportion
  # dim: Integer constant, the dimension of the starting lattice.
  # nei: Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  # p: Real constant between zero and one, the rewiring probability.
  # vis_factor: Visibility factor
  # mem_factor: Memory factor
  # sub_mem_factor: Subpopulation memory factor
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  # SIR generation of the hidden population 
  hidden_population = gen_SIRpop(N, net_sw ,beta, gamma, chosen_nodes, infected_people)
  
  # First we generate the population, then we remplace the apropiate column
  Population = genPopulation(N, prob_vect,prob_hp)
  Population$hidden_population = hidden_population

  n_populations = length(prob_vect)
  # initializes the vectors
  vect_hp = rep(NA,N)       # number of hidden population individuals known for each person
  vect_hp_vis = rep(NA,N)   # vect_hp applying visibility
  vect_reach = rep(NA,N)    # the degrees of each individual
  vect_reach_re = rep(NA,N) # reach vector applying memory error
  
  # Matrix representing the directed graph that connects individuals with the people of the Hidden Population they know 
  Mhp = matrixHP(net_sw,Population)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = vis_factor)
  
  # Mhp_hp calculation for making two different Bernouillis
  Mhp_vis_comp = 2*Mhp_vis - 1
  Mhp_hp = 1*(t(Mhp_vis) ==  Mhp_vis_comp)
  Mhp_hp_vis = apply(Mhp_hp, c(1,2), berHP, p = (0.5 +  0.5*vis_factor))
  
  # Matrix for knowing with bernoilli should be applied
  Mhp_sir = Mhp_vis + 2*Mhp_hp + 2*Mhp_hp_vis
  print(Mhp_sir)
  # Final visibility matrix
  Mhp_vis = apply(Mhp_sir, c(1,2), to_matrix_SIR)
  
  for (i in 1:N) {
    # net_sw[[i]], list with one element, the list of the adjacent vertices to i
    
    vect_hp[i] = sum(Mhp[i,])
    vect_reach[i] = length(net_sw[[i]][[1]])
    vect_hp_vis[i] = round(rtruncnorm(1, a = max(-0.5,  2 * sum(Mhp_vis[i,]) - vect_reach[i] + 0.5 ) , b = min(2 * sum(Mhp_vis[i,]) + 0.5, vect_reach[i]-0.5), mean = sum(Mhp_vis[i,]), sd = mem_factor*sum(Mhp_vis[i,])))
    
    vect_reach_re[i] = max(1,round(rtruncnorm(1, a = vect_hp_vis[i] - 0.5, b = 2*vect_reach[i] - vect_hp_vis[i] + 0.5, mean = vect_reach[i], sd = mem_factor*vect_reach[i])))
  }
  
  Population = cbind(Population, reach = vect_reach)
  Population = cbind(Population, reach_memory = vect_reach_re)
  Population = cbind(Population, hp_total = vect_hp) 
  Population = cbind(Population, hp_survey = vect_hp_vis)
  
  for(j in 1:length(prob_vect)){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(as.character(j)))) 
      vis_yij = sum(Population[net_sw[[i]][[1]],]["hidden_population"][as.logical(dplyr::select(Population[net_sw[[i]][[1]],],starts_with("subpop") & ends_with(as.character(j)))[,1]),]) 
      # Visibility of population j by i, applying a normal in order to represent the real visibility
      
      v_1[i] = max(0,round(rtruncnorm(1, a = vis_yij - 0.5 , b = 2*vis_pob - vis_yij + 0.5,  mean = vis_pob, sd = sub_mem_factor*vis_pob)))
    }
    
    Population = cbind(Population,Subpoblacion_total = v_1)
    names(Population)[dim(Population)[2]] = str_c("kp_reach_",j)
  }
  
  ind1 = 1:N
  i_hp_vis = rep(NA,N)
  for (i in 1:length(prob_vect)) {
    for (j in ind1){
      ind2 = dplyr::select(Population, starts_with("subpop") & ends_with(as.character(i)))[,1] != 0
      i_hp_vis[j] = round(rtruncnorm(1, a = -0.5, b =  2*sum(Mhp_vis[ind2,j]) + 0.5, mean = sum(Mhp_vis[ind2,j]), sd = sum(Mhp_vis[ind2,j])*sub_mem_factor)) 
    }
    Population = cbind(Population, Subpoblacion_total = i_hp_vis)
    names(Population)[dim(Population)[2]] = str_c("kp_alters_",i)
  }
  
  
  return(list(net_sw, Population, Mhp_vis))
}

################################################################################

N = 1000                      # Population size
v_pop_prob = rep(1/10, 5)     # Probability of each subpopulation. As we are working with disjoint and no disjoint subpopulations
# sum(v_pop_prob) < 1. 
n_pop = length(v_pop_prob)    # Number of subpopulations
hp_prob = 0.1                 # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 50                 # Number of individuals we draw in the survey
n_survey_hp = 50              # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0         # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1         # Visibility factor (Binomial's probability)
memory_factor = 0             # Reach memory factor (parameter to change variance of the perturbations' normal)

seed = 207                    # Seed
set.seed(seed)

#Graph
dim = 1      # Graph dimension 
nei = 75     # Number of neighbors that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
             # before applying the randomization.
p   = 0.1    # Probability of randomize a connection. It is applied to all connections


# Variables #
beta = 0.1                              # Set the infection rate
gamma = 0                               # Set the removal rate
infected_people = 100                   # Set number of infected people
chosen_nodes = seq(1, 1000, 100)        # Set of nodes from which we start


Graph_population_matrix = getData(N, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor, beta, gamma, chosen_nodes, infected_people)

net_sw = Graph_population_matrix[[1]]       # Population´s graph
Population = Graph_population_matrix[[2]]   # Population
Mhp_vis = Graph_population_matrix[[3]]      # Population's visibility matrix
