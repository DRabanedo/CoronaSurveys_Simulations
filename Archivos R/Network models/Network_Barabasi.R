#Dinamic mode

#In this model, an edge is most likely to attach to nodes with higher degrees. 
#The network begins with an initial network of m nodes.
#In the Barabasi-Albert model, new nodes are added to the network one at a time.
#Each new node is connected to m existing nodes with a probability that is 
#proportional to the number of links that the existing nodes already have.
net_bardin <- barabasi.game(500, power = 1.2, m = NULL, out.dist = NULL, out.seq = NULL,
                   out.pref = FALSE, zero.appeal = 1, directed = FALSE,
                   algorithm ="psumtree", start.graph = NULL)
plot(net_bardin, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Scale-free network model")

#Static mode

#In this case static.power.law.game() function can be used with known power law 
#exponent, # of edges, # of nodes.
net_barstat <- static.power.law.game(500, 500, exponent.out= 2.2, exponent.in = -1, loops = FALSE, multiple = FALSE, finite.size.correction = TRUE) 
plot(net_barstat, vertex.label= NA, edge.arrow.size=0.02,vertex.size = 0.5, xlab = "Scale-free network model (static)")





