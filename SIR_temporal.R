#####################
# Libraries used #  #
library(igraph)     #
library(stringr)    #
library(ggplot2)    #
library(sampler)    #
library(dplyr)      #
#####################

############################## Functions #######################################

gen_SIRpop = function(N, net ,beta, gamma, chosen_nodes, infected_recovered){
  # SIR method for determination of the hidden population
  
  # beta: infection rate
  # gamma: removal rate
  # infected_people: number of infected people
  # chosen_nodes: nodes from which we start
  
  # Loop variables
  infected_nodes = c()
  infected_nodes = c(infected_nodes,chosen_nodes)
  
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

  return(list(infected_nodes, infected_recovered))
}


get_Time_matrix = function(N, nei, t, dim = 1, p = 0.1, beta = 0.1, gamma = 0.03, chosen_nodes = 1){
  ## Function parameters ##
  
  # t       Number of time iterations
  # nei     Number of neighbors (mean value)
  # dim     Dimension of the network
  # p       Probability of randomising a connection
  # beta    Beta parameter of SIR model
  # gamma   Gamma parameter of SIR model
  
  # Network creation #
  network     = sample_smallworld(dim, size = N, nei, p, loops = FALSE, multiple = FALSE)
  data_matrix = matrix(rep(NA, t*N*2), nrow = N, ncol = t*2)
  
  # Loop variables
  infected_recovered = c()
  
  for (k in 1:t){
    # Nodes infected by the SIR model #
    SIR_data             = gen_SIRpop(N, network, beta, gamma, chosen_nodes, infected_recovered)
    final_infected_nodes = SIR_data[[1]]
    infected_recovered   = SIR_data[[2]]
    
    # Infected vector for the next step
    chosen_nodes = final_infected_nodes
    
    # Hidden population links 
    links_hp = rep(NA, N)
    for (j in 1:N){
      count = 0
      for (l in network[[j]][[1]]){
        if (as.logical(sum(l %in% final_infected_nodes))){
          count = count + 1
        }
      }
      links_hp[j] = count
    }
    
    # Output data creation #
    hp_vector = rep(NA, N)
    for (i in 1:N){
      if (as.logical(sum(i %in% final_infected_nodes))){
        hp_vector[i] = 1
      }
      else{
        hp_vector[i] = 0
      }
    }
    
    # Final data matrix #
    data_matrix[,(2*k - 1)] = hp_vector
    data_matrix[,(2*k)]     = links_hp

  }
  
  return(data_matrix)
}

################################################################################
# Simulationes en tiempo #

# Variables para la estructura de la red #
N   = 10000         # Tamaño de la población a estudiar
nei = 36            # Media de vecinos que tiene cada nodo (introducir un número par)

# Variables del modelo SIR
R_0     = 3         # R_0 del modelo SIR, en este caso elegimos 3
beta    = 0.03      # Probabilidad de infectar a un vecino en cada unidad de tiempo
gamma   = beta/R_0  # Probabilidad de que una persona enferma se recupere en cada unidad de tiempo
n_nodes = 1         # Número de nodos infectados al comenzar el proceso (número de focos)

t   = 30   # Número de iteraciones temporales

# Generación de los datos #
matrix_data = get_Time_matrix(N, round(nei/2), t, beta = beta, gamma = gamma, dim = 1, p = 0.1, chosen_nodes = round(seq(1,N,N/n_nodes)))

################################################################################

# Gráfica de la proporción de la población oculta en función del tiempo
hp_prop  = rep(NA, t)
for (i in 1:t){
  hp_prop[i] = sum(matrix_data[,2*i - 1])/N
}

ggplot() + 
  geom_line(aes(x = 1:t, y =  hp_prop)) +
  scale_color_discrete("Legend") + 
  labs(title = "Infected people distribution",
       x = "Time",
       y = "Infected people proportion")

################################################################################







