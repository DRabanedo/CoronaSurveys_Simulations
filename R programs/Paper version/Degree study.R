# Network degree analysis #

##### Packages #####
library(igraph)    #
library(stringr)   #
library(ggplot2)   #
####################

############################# igraph library ###################################

# Info

# sample_smallworld(dim, size, nei, p, loops = FALSE, multiple = FALSE)
# dim   #Integer constant, the dimension of the starting lattice.
# size  #Integer constant, the size of the lattice along each dimension.
# nei   #Integer constant, the neighborhood within which the vertices of 
#the lattice will be connected.
# p     #Real constant between zero and one, the rewiring probability.

# loops	    #Logical scalar, whether loops edges are allowed in the generated graph.
# multiple	#Logical scalar, whether multiple edges are allowed int the generated graph.

#First a lattice is created with the given dim, size and nei arguments. 
#Then the edges of the lattice are rewired uniformly randomly with probability p.

#Note that this function might create graphs with loops and/or multiple edges. 
#You can use simplify to get rid of these.

################################################################################

## Variables ##

# Seed #
seed = 207
set.seed(seed)

# Network #
dim = 1
size = 10000
nei = 75
loops = FALSE
multiple = FALSE



# Loop #
iterations = 1
p_iter = 0.1

degree_df = data.frame(row.names = c(1:(iterations*size)))
degree_vect = c()
degree_loop = rep(NA, size)


# Degree calculation #
for (h in 1:length(p_iter)){
  # Variable reset #
  degree_vect = c()
  degree_loop = rep(NA, size)
  
  for (i in 1:iterations){
    network = sample_smallworld(dim, size, nei, p = p_iter[h], loops, multiple)
    for (j in 1:size){
      degree_loop[j] = length(network[[j]][[1]])
    }
    
    degree_vect = append(degree_vect, degree_loop)
  }
  degree_df = cbind(degree_df, Probability = degree_vect)
  names(degree_df)[dim(degree_df)[2]] = str_c("Probability_", p_iter[h])
  
}
getwd()

# Graph representation #

degree_summary =  degree_df$Probability_0.1
degree_var = round(var(degree_summary), digits = 2)
degree_median = median(degree_summary)
degree_max = max(degree_summary)
degree_min = min(degree_summary)

sub_title = str_c("Probability = ", p_iter, ", median = ", degree_median, ", var = ", degree_var,", min = ", degree_min, ", max = ", degree_max, ". Small World model")
plot_name = str_c("Smallword_", p_iter,".png")

png(filename = plot_name,
    width = 1000, height = 600)

degree_graph = ggplot(degree_df) +
  geom_histogram( aes(x = Probability_0.1, y = ..count../sum(..count..)), binwidth = 1, color = "black", fill = "grey", alpha = 0.4) +
  #geom_histogram( aes(x = Probability_0.5, y=..count../sum(..count..)), binwidth = 2, color = "red", fill = "red", alpha = 0.2) +
  #geom_histogram( aes(x = Probability_0.75, y=..count../sum(..count..)), binwidth = 2, color = "blue", fill = "blue", alpha = 0.2) +
  #geom_histogram( aes(x = Probability_1, y=..count../sum(..count..)), binwidth = 2, color = "yellow", fill = "yellow", alpha = 0.2) +
  labs(title = "Network degree simulation",
       subtitle = sub_title,
       x = "Number of neighbors",
       y = "Proportion")

degree_graph

dev.off()
