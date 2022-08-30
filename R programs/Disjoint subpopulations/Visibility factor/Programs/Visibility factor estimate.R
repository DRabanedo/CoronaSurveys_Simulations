
##################################################################################
# In this simulation we represent a way to simulate the process of obtaining an  #
# estimate of the visibility factor. We use the formulas used in coming-out-rate # 
# and game of contacts calculus. We also introduce some conditions that made the #
# answers obtained on the survey consistent.                                     #
##################################################################################

#####################################################
# Visibility factor prediction using subpopulations #
#####################################################

vf_subpop_es = function(survey_hp,Population, Mhp_vis,sub_memory_factor){
  # Total sums variables creation
  sum_pop = 0
  sum_pop_hp = 0
  
  # Number of subpopulations calculus 
  n_pop = length(names(dplyr::select(Population, starts_with("KP_"))))-1
  
  # Index of the people who has been surveyed
  ind_survey = as.numeric(rownames(survey_hp))
  
  # Loop for studying each answered value for each subpopulation on the survey 
  for (k in 1:n_pop) {
    
    #People who belong to the subpopulation k
    ind_subpop = as.numeric(rownames(Population[Population$Population == k,]))
    
    # Vector representing how many people of subpopulation k knows that the surveyed
    # person belongs to the hidden population
    vect_pop_hp = colSums(Mhp_vis[ind_subpop,ind_survey])
    
    # Application of a memory factor to that answer
    mem_vect_pop_hp = rep(NA,length(vect_pop_hp))
    for (i in 1:length(vect_pop_hp)) {
      mem_vect_pop_hp[i] = max(0,round(rnorm(1,mean = vect_pop_hp[i],sub_memory_factor*vect_pop_hp[i])))
    }
    
    # Vector who represent how many people each surveyed person knows from subpopulation k 
    # (with memory factor applied)
    mem_vect_pop = dplyr::select(Population, starts_with("KP_") & ends_with(as.character(k)))[,1][ind_survey]
    mem_vect_pop
    for (j in 1:length(mem_vect_pop)) {
      
      # As the people a person knows is always bigger that the people that a person knows
      # AND knows its belonging to the hidden population, we make a loop that stops
      # when the assumption commented is verified.
      while (mem_vect_pop_hp[j] > mem_vect_pop[j]) {
        
        mem_vect_pop_hp[j] = max(0,round(rnorm(1,mean = vect_pop_hp[j],sub_memory_factor*vect_pop_hp[j])))
        
        # It makes it converge easier (reduces computation time)
        if (mem_vect_pop_hp[j]<vect_pop_hp[j])
          vect_pop_hp[j] = mem_vect_pop_hp[j]
      }
      
    }
    
    # Sum of the results on subpopulation k to the corresponding general variables
    sum_pop_hp = sum_pop_hp + sum(mem_vect_pop_hp)
    sum_pop = sum_pop + sum(mem_vect_pop)
    
  }
  
  
  # Visibility factor estimate
  vf_subpop = sum_pop_hp/sum_pop
  
  return(vf_subpop)
}

################################################################################

#Outliers detection
#####################################################
# Visibility factor prediction using subpopulations #
#####################################################

vf_subpop_es_out = function(survey_hp,Population, Mhp_vis,sub_memory_factor){
  # Total sums variables creation
  sum_pop = 0
  sum_pop_hp = 0
  
  # Number of subpopulations calculus 
  n_pop = length(names(dplyr::select(Population, starts_with("KP_"))))-1
  
  # Index of the people who has been surveyed
  ind_survey = as.numeric(rownames(survey_hp))
  
  # Loop for studying each answered value for each subpopulation on the survey 
  for (k in 1:n_pop) {
    
    #People who belong to the subpopulation k
    ind_subpop = as.numeric(rownames(Population[Population$Population == k,]))
    
    # Vector representing how many people of subpopulation k knows that the surveyed
    # person belongs to the hidden population
    vect_pop_hp = colSums(Mhp_vis[ind_subpop,ind_survey])
    
    # Application of a memory factor to that answer
    mem_vect_pop_hp = rep(NA,length(vect_pop_hp))
    for (i in 1:length(vect_pop_hp)) {
      mem_vect_pop_hp[i] = max(0,round(rnorm(1,mean = vect_pop_hp[i],sub_memory_factor*vect_pop_hp[i])))
    }
    
    # Vector who represent how many people each surveyed person knows from subpopulation k 
    # (with memory factor applied)
    mem_vect_pop = dplyr::select(Population, starts_with("KP_") & ends_with(as.character(k)))[,1][ind_survey]
    mem_vect_pop
    for (j in 1:length(mem_vect_pop)) {
      
      # As the people a person knows is always bigger that the people that a person knows
      # AND knows its belonging to the hidden population, we make a loop that stops
      # when the assumption commented is verified.
      if (mem_vect_pop_hp[j] > mem_vect_pop[j]){
        mem_vect_pop_hp[j] =  0
        mem_vect_pop[j] = 0
      }
      
    }
    
    # Sum of the results on subpopulation k to the corresponding general variables
    sum_pop_hp = sum_pop_hp + sum(mem_vect_pop_hp)
    sum_pop = sum_pop + sum(mem_vect_pop)
    
  }
  
  
  # Visibility factor estimate
  vf_subpop = sum_pop_hp/sum_pop
  
  return(vf_subpop)
}

################################################################################

#############################################
# Visibility factor predictions using reach #
#############################################

vf_reach_es = function(survey_hp,Population, Mhp_vis, memory_factor) {
  # People from the general population who know that the people from the survey_hp belong to the hidden population
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach_hp[i] = max(0,round(rnorm(1,mean = vect_reach_hp[i],memory_factor*vect_reach_hp[i])))
  }
  
  
  
  #People from the subpopulations known by the hidden population in survey_hp  
  mem_vect_reach = Population$Reach_memory[ind_survey]
  
  for (i in 1:length(mem_vect_reach)) {
    while (mem_vect_reach[i] < mem_vect_reach_hp[i]){
      mem_vect_reach_hp[i] = max(0,round(rnorm(1,mean = vect_reach_hp[i],memory_factor*vect_reach_hp[i])))
      
      if (mem_vect_reach_hp[i]<vect_reach_hp[i])
        vect_reach_hp[i] = mem_vect_reach_hp[i]
    }
    
  }
  
  # Final estimate
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}


################################################################################

#Outliers detection
#############################################
# Visibility factor predictions using reach #
#############################################

vf_reach_es_out = function(survey_hp,Population, Mhp_vis, memory_factor) {
  # People from the general population who know that the people from the survey_hp belong to the hidden population
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach_hp[i] = max(0,round(rnorm(1,mean = vect_reach_hp[i],memory_factor*vect_reach_hp[i])))
  }
  
  
  
  #People from the subpopulations known by the hidden population in survey_hp  
  mem_vect_reach = Population$Reach_memory[ind_survey]
  
  for (j in 1:length(mem_vect_reach)) {
    if (mem_vect_reach_hp[j] > mem_vect_reach[j]){
      mem_vect_reach_hp[j] =  0
      mem_vect_reach[j] = 0
    }
    
  }
  
  # Final estimate
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}

################################################################################

# Population generation #

N = 1000                 # Population size
v_pop = c(0:10)         # Subpopulations vector. They are disjoint and 0 corresponds to not classifying the individual in any of them
n_pop = length(v_pop)-1    # Number of subpopulations
v_pop_prob =  rep(1/length(v_pop), length(v_pop)) #Probability of each subpopulation
hp_prob = 0.1             # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 300            # Number of individuals we draw in the survey
n_survey_hp = 50          # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0.1     # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
memory_factor = 0.1         # Reach memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 0.9     # Visibility factor (Binomial's probability)

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

#Surveys
survey = getSurvey(300, Population)
survey_hp = getSurvey(50,Population[Population$Hidden_Population == 1,])


################################################################################


# Predictions

vf_subpop_es(survey_hp, Population, Mhp_vis, sub_memory_factor)
vf_subpop_es_out(survey_hp, Population, Mhp_vis, sub_memory_factor)

vf_reach_es(survey_hp, Population, Mhp_vis, memory_factor)
vf_reach_es_out(survey_hp, Population, Mhp_vis, memory_factor)

visibility_factor
