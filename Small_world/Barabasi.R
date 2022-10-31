# Network degree analysis #

##### Packages #####
library(igraph)    #
library(stringr)   #
library(ggplot2)   #
####################

############################# igraph library ###################################

# Info

# sample_smallworld(dim, size, nei, p, loops = FALSE, multiple = FALSE)

# n           Number of vertices.
# power	      The power of the preferential attachment, the default is one, ie. linear preferential attachment.
# m	          Numeric constant, the number of edges to add in each time step This argument is only used if both 
              # out.dist and out.seq are omitted or NULL.
# out.dist  	Numeric vector, the distribution of the number of edges to add in each time step. This argument is
              # only used if the out.seq argument is omitted or NULL.
# out.seq	    Numeric vector giving the number of edges to add in each time step. Its first element is ignored as 
              # no edges are added in the first time step. 
# out.pref	  Logical, if true the total degree is used for calculating the citation probability, otherwise the in-degree is used.
# zero.appeal The ‘attractiveness’ of the vertices with no adjacent edges. See details below.
# directed	  Whether to create a directed graph.
# algorithm	  The algorithm to use for the graph generation. psumtree uses a partial prefix-sum tree to generate the graph, this 
              # algorithm can handle any power and zero.appeal values and never generates multiple edges. psumtree-multiple also 
              # uses a partial prefix-sum tree, but the generation of multiple edges is allowed. Before the 0.6 version igraph used 
              # this algorithm if power was not one, or zero.appeal was not one. bag is the algorithm that was previously (before version 0.6) 
              # used if power was one and zero.appeal was one as well. It works by putting the ids of the vertices into a bag 
              # (multiset, really), exactly as many times as their (in-)degree, plus once more. Then the required number of cited vertices are 
              # drawn from the bag, with replacement. This method might generate multiple edges. It only works if power and zero.appeal are equal one.

# start.graph	NULL or an igraph graph. If a graph, then the supplied graph is used as a starting graph for the preferential attachment
              # algorithm. The graph should have at least one vertex. If a graph is supplied here and the out.seq argument is not NULL, 
              # then it should contain the out degrees of the new vertices only, not the ones in the start.graph.




# sample_pa generates a directed graph by default, set directed to FALSE to generate an undirected graph. Note that even if an undirected graph is generated
# denotes the number of adjacent edges not initiated by the vertex itself and not the total (in- + out-) degree of the vertex, unless the out.pref argument is set to TRUE.

################################################################################

## Variables ##

# Seed #
seed = 207
set.seed(seed)
getwd()

# Network #
size = 10000
size_sw = 0.1
nei = 20
n_links = 9

links =  c(rep(30, size * (1-size_sw)/n_links), rep(40, size * (1-size_sw)/n_links), rep(50, size * (1-size_sw)/n_links),
          rep(60, size * (1-size_sw)/n_links), rep(70, size * (1-size_sw)/n_links), rep(80, size * (1-size_sw)/n_links), 
          rep(90, size * (1-size_sw)/n_links),  rep(65, 1/4* size * (1-size_sw)/n_links), rep(70, 1/4 * size * (1-size_sw)/n_links), 
          rep(75, 1/4* size * (1-size_sw)/n_links), rep(80, 1/4* size * (1-size_sw)/n_links), rep(85, 1/4* size * (1-size_sw)/n_links),
          rep(90, 1/4* size * (1-size_sw)/n_links), rep(95, 1/4* size * (1-size_sw)/n_links), rep(100, 1/4* size * (1-size_sw)/n_links))
  
m = 80
power_sa = seq(0.4,1.5,0.1)

# Loop #
degree_df = data.frame(row.names = c(1:size))
degree_vect = c()
degree_loop = rep(NA, size)

# Using sample_pa #
par(mfrow=c(3,4))

for (i in 1:length(power_sa)){
  network = sample_pa(n = size, 
                      out.seq = links,
                      #m = m,
                      power = power_sa[i],
                      directed = FALSE,
                      algorithm = 'psumtree',
                      start.graph = sample_smallworld(dim = 1, size = size*size_sw, nei = nei, p = 0, loops = FALSE, multiple = FALSE)
                      )
  for (j in 1:size){
    degree_loop[j] = length(network[[j]][[1]])
  }
    
  degree_vect = append(degree_vect, degree_loop)
  title = str_c("Power = ", power_sa[i])
  hist(degree_vect, main = title, ylab = "Count", xlab = "Neighbors", breaks = 50)
}

################################################################################

## Variables ##

# Seed #
getwd()

# Network #
size = 10000
size_sw = 0.1
nei = 10

links =  round( 1.0658 * c(rep(40, size * (1-size_sw)/n_links), rep(50, size * (1-size_sw)/n_links),
           rep(60, size * (1-size_sw)/n_links), rep(70, size * (1-size_sw)/n_links), rep(80, size * (1-size_sw)/n_links), 
           rep(90, size * (1-size_sw)/n_links), rep(100, size * (1-size_sw)/n_links),rep(50, 1/16* size * (1-size_sw)/n_links), rep(55, 1/16 * size * (1-size_sw)/n_links), 
           rep(60, 1/16* size * (1-size_sw)/n_links), rep(65, 1/16* size * (1-size_sw)/n_links), rep(70, 1/16 * size * (1-size_sw)/n_links), 
           rep(75, 1/16* size * (1-size_sw)/n_links), rep(80, 1/8* size * (1-size_sw)/n_links), rep(85, 1/8* size * (1-size_sw)/n_links),
           rep(90, 1/8* size * (1-size_sw)/n_links), rep(95, 1/4* size * (1-size_sw)/n_links), rep(100, 1/4* size * (1-size_sw)/n_links), 
           rep(110, 1/4* size * (1-size_sw)/n_links), rep(120, 1/4* size * (1-size_sw)/n_links), rep(130, 1/4* size * (1-size_sw)/n_links), 60, 60, 68))

length(links)

sum(links) + nei*2*size_sw*size
power_sa = 1

# Loop #
degree_df = data.frame(row.names = c(1:size))
degree_vect = c()
degree_loop = rep(NA, size)

# Using sample_pa #
par(mfrow=c(1,1))
for (i in 1:length(power_sa)){
  network = sample_pa(n = size, 
                      out.seq = links,
                      #m = m,
                      power = power_sa[i],
                      directed = FALSE,
                      algorithm = 'psumtree',
                      start.graph = sample_smallworld(dim = 1, size = size*size_sw, nei = nei, p = 0, loops = FALSE, multiple = FALSE)
  )
  for (j in 1:size){
    degree_loop[j] = length(network[[j]][[1]])
  }
  
  degree_vect = append(degree_vect, degree_loop)
  title = str_c("Power = ", power_sa[i])
  hist(degree_vect, main = title, ylab = "Count", xlab = "Neighbors", breaks = 50)
}

degree_summary = degree_vect
degree_var     = round(var(degree_summary), digits = 2)
degree_median  = median(degree_summary)
degree_max     = max(degree_summary)
degree_min     = min(degree_summary)
degree_mean    = mean(degree_summary)

degree_distribution(network)

sub_title = str_c("Power = ", power_sa,", mean = ", degree_mean, ", median = ", degree_median, ", var = ", degree_var,", min = ", degree_min, ", max = ", degree_max, ". Barabasi model")
plot_name = str_c("Smallword_", p_iter,".png")
degree_df = data.frame(Probability_0.1 = degree_vect)

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
