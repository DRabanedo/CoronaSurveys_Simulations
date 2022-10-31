#######################
# Libraries used #    #
library(igraph)       #
library(tidyverse)    #
library(stringr)      #
library(ggplot2)      #
library(sampler)      #
library(dplyr)        #
library(truncnorm)    #
#######################

# Function #

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

# Network construction #
dim = 1
nei = 75
p = 0.1
N = 10000
net = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)


# Variables #

beta = 0.1               # Set the infection rate
gamma = 0                # Set the removal rate
infected_people = 1000     # Set number of infected people
chosen_nodes = seq(1, 10000, 1000)         # Set of nodes from which we start

gen_SIRpop(N, net ,beta, gamma, chosen_nodes, infected_people)

