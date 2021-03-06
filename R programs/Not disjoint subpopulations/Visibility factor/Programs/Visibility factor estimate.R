##############################
# Visibility factor estimate #
##############################


N = 1000                 # Population size
v_pop = c(1:10)           # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)   # Number of subpopulations
v_pop_prob =  rep(1/length(v_pop), length(v_pop)) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0         # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 1     # Visibility factor (Binomial's probability)
seed = 207                # Seed
set.seed(seed)


#Graph
dim = 1    # Graph dimension 
nei = 75   # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1  # Probability of randomize a connection. It is applied to all connections


#Population and Survey

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]      # Population´s graph
Population = Graph_population_matrix[[2]]  # Population
Mhp_vis = Graph_population_matrix[[3]]     # Population's visibility matrix

survey_hp = Population[Population$Hidden_Population == 1,][sample(nrow(Population[Population$Hidden_Population == 1,]), n_survey_hp, replace = FALSE),]


v_pop_total = rep(NA, n_pop)
for (k in 1:n_pop) {
  v_pop_total[k] = sum(Population[,k+1]) # N_k
}


##### Visibility factor estimate (using subpopulations) #####

# People from the subpopulations who know that the people from the survey_hp belong to the hidden population

vf_subpop_es = function(survey_hp,Population, Mhp_vis){

    ind_survey = as.numeric(rownames(survey_hp))
  sum_pop_hp = 0
  for (i in 1:length(v_pop_total)) {
    ind_subpop = Population[,i+1] != 0
    sum_pop_hp = sum(Mhp_vis[ind_subpop,ind_survey]) + sum_pop_hp
  }
  
  #People from the subpopulations known by the hidden population in survey_hp  
  sum_pop = sum(select(Population, starts_with("KP_"))[ind_survey,])
  
  # Final estimate
  vf_subpop = sum_pop_hp/sum_pop
  
  return(vf_subpop)
}

vf_subpop = vf_subpop_es(survey_hp,Population, Mhp_vis)




##### Visibility factor estimate (using Reach) #####

# People from the general population who know that the people from the survey_hp belong to the hidden population
vf_reach_es = function(survey_hp,Population, Mhp_vis) {
  
  ind_survey = as.numeric(rownames(survey_hp))
  sum_reach_hp = sum(Mhp_vis[,ind_survey])
  
  #People from the subpopulations known by the hidden population in survey_hp  
  sum_reach = sum(Population$Reach_memory[ind_survey])
  
  # Final estimate
  vf_reach = sum_reach_hp/sum_reach
  return(vf_reach)
}

vf_reach = vf_reach_es(survey_hp,Population, Mhp_vis)



vf_subpop
vf_reach
visibility_factor
