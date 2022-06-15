library(igraph)
erdos.renyi.game(n, p.or.m, type = c("gnp", "gnm"), directed = FALSE, loops = FALSE)

n = 50      #Número de nodos en el grafo
p.or.m = p  #Either the probability for drawing an edge between two arbitrary 
#vertices (G(n,p) graph), or the number of edges in the graph 
#(for G(n,m) graphs).
type = gnp  #The type of the random graph to create, either gnp (G(n,p) graph) 
#or gnm (G(n,m) graph).
directed    #Si el grafo es dirigido, TRUE o FALSE.
loops       #Logical, whether to add loop edges, defaults to FALSE

# G(n,p) graphs, the graph has ‘n’ vertices and for each edge the probability 
#that it is present in the graph is ‘p’.

#In G(n,m) graphs, the graph has ‘n’ vertices and ‘m’ edges, and the ‘m’ edges 
#are chosen uniformly randomly from the set of all possible edges. This set
#includes loop edges as well if the loops parameter is TRUE.

net_gnm <- erdos.renyi.game(500, 350, type = "gnm")
plot(net_gnm, xlab = "Random Network: G(N,L) model")
#Quitando lo que sobra para ver bien el grafo
plot(net_gnm, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Random Network: G(N,m) model")

net_gnp <- erdos.renyi.game(50, 0.35, type = "gnp")
plot(net_gnp, xlab = "Random Network: G(N,p) model")
#Quitando lo que sobra para ver bien el grafo
plot(net_gnp, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Random Network: G(N,p) model")


