library(igraph)
sample_smallworld(dim, size, nei, p, loops = FALSE, multiple = FALSE)
dim = 1   #Integer constant, the dimension of the starting lattice.
size= 50  #Integer constant, the size of the lattice along each dimension.
nei = 1   #Integer constant, the neighborhood within which the vertices of 
          #the lattice will be connected.
p = 0.1   #Real constant between zero and one, the rewiring probability.

loops	    #Logical scalar, whether loops edges are allowed in the generated graph.
multiple	#Logical scalar, whether multiple edges are allowed int the generated graph.

#First a lattice is created with the given dim, size and nei arguments. 
#Then the edges of the lattice are rewired uniformly randomly with probability p.

#Note that this function might create graphs with loops and/or multiple edges. 
#You can use simplify to get rid of these.

net_sw = sample_smallworld(1, 20, 4, 0.1, loops = FALSE, multiple = FALSE)
plot(net_sw, xlab = "Small world model")
net_sw

net_sw2 <- watts.strogatz.game(1, 500, 1, 0.35, loops = FALSE, multiple = FALSE)
plot(net_sw2, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Small world model")

# Verification
total_edges = 0
for (i in 1:20) {
  total_edges = total_edges + length(net_sw[[i]][[1]])
  print(length(net_sw[[i]][[1]]))
}
total_edges/20
