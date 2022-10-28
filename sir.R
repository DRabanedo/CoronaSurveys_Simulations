dim=1
nei=15
p=0.1
N=100
net = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
infected_nodes = c()
infected_recovered = c()
# Set the infection rate
beta = 0.1
# Set the removal rate
gamma = 0.01
# Set how many timesteps you want to pass through
n_timesteps = 20
# Start from the node you have chosen
chosen_node = 1
infected_nodes = c(infected_nodes,chosen_node)
for (i in 1:n_timesteps) {
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
