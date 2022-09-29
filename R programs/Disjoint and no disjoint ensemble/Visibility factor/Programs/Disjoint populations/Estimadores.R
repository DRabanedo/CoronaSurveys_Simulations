# Visibility factor estimate #

####################
# Binomial aproach #
####################

vf_bin = function(survey_hp,Population, Mhp_vis,memory_factor){
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double trunc + binomial calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rnorm(1 , mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = rbinom(1, size = mem_vect_reach[i], p = vect_reach_hp[i]/vect_reach[i])
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}




###########################
# Double truncation first #
###########################

################################################################################

# Double truncation + binomial #

vf_dt_bin = function(survey_hp,Population, Mhp_vis,memory_factor){
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
 
  # Double trunc + binomial calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + 2 * vect_reach[i] - vect_reach_hp[i], mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = rbinom(1, size = mem_vect_reach[i], p = vect_reach_hp[i]/vect_reach[i])
  }

  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}



################################################################################

# Double truncation + fake binomial #

vf_dt_fbin = function(survey_hp,Population, Mhp_vis,memory_factor){
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double trunc + binomial calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + 2 * vect_reach[i] - vect_reach_hp[i], mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = rbinom(1, size = mem_vect_reach[i], p = vect_reach_hp[i]/mem_vect_reach[i])
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}


################################################################################

# Double truncation + double truncation #

vf_dt_dt = function(survey_hp, Population, Mhp_vis,memory_factor){
  
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double truncation + double truncation calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i]    = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + 2 * vect_reach[i] - vect_reach_hp[i], mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = max(0,round(rtruncnorm(1, a = vect_reach_hp[i] - (mem_vect_reach[i]-vect_reach_hp[i]), b = mem_vect_reach[i], mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}

################################################################################

# Double truncation + unilateral truncation # 

vf_dt_ut = function(survey_hp,Population, Mhp_vis,memory_factor){
  
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i]    = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + 2 * vect_reach[i] - vect_reach_hp[i], mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = max(0,round(rtruncnorm(1, a =  0 , b = mem_vect_reach[i], mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}

################################################################################

# Double truncation + unilateral truncation corrected # 

vf_dt_utc = function(survey_hp,Population, Mhp_vis,memory_factor){
  
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i]    = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + 2 * vect_reach[i] - vect_reach_hp[i], mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = max(0,round(rnorm(1, mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))
    if (mem_vect_reach_hp[i] > mem_vect_reach[i]){
      while ((mem_vect_reach_hp[i] > mem_vect_reach[i]) | (mem_vect_reach_hp[i]  <= vect_reach_hp[i])) {
        mem_vect_reach_hp[i] = max(0,rnorm(1, mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor))
      }
    }
  }
  
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}

################################################################################



###############################
# Unilateral truncation first #
###############################

################################################################################

# Unilateral truncation + binomial #

vf_ut_bin = function(survey_hp,Population, Mhp_vis,memory_factor){
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double trunc + binomial calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + Inf, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = rbinom(1, size = mem_vect_reach[i], p = vect_reach_hp[i]/vect_reach[i])
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}



################################################################################

# Unilateral truncation + binomial #

vf_ut_fbin = function(survey_hp,Population, Mhp_vis,memory_factor){
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double trunc + binomial calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + Inf, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = rbinom(1, size = mem_vect_reach[i], p = vect_reach_hp[i]/mem_vect_reach[i])
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}
################################################################################


# Unilateral truncation + double truncation #

vf_ut_dt = function(survey_hp,Population, Mhp_vis,memory_factor){
  
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double truncation + double truncation calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b = (-1) + Inf, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = max(0,round(rtruncnorm(1, a = vect_reach_hp[i] - (mem_vect_reach[i]-vect_reach_hp[i]), b =  mem_vect_reach[i], mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}


################################################################################


# Unilateral truncation + unilateral truncation # 

vf_ut_ut = function(survey_hp,Population, Mhp_vis,memory_factor){
  
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double truncation + unilateral truncation calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rtruncnorm(1, a = 1 + vect_reach_hp[i] , b =Inf, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    mem_vect_reach_hp[i] = max(0,round(rtruncnorm(1, a = 0 , b = mem_vect_reach[i], mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}

################################################################################


################################################
# Unilateral truncation  with correction first #
################################################

# Unilateral truncation corrected + unilateral truncation corrected #

vf_utc_utc = function(survey_hp,Population, Mhp_vis,memory_factor){
  
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    if (mem_vect_reach[i]  <= vect_reach_hp[i]){
      while ((mem_vect_reach[i]  <= vect_reach_hp[i]) | (mem_vect_reach[i] > vect_reach[i])) {
        mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
      }
    }
    
    mem_vect_reach_hp[i] = max(0,round(rnorm(1, mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))
    if (mem_vect_reach_hp[i] > mem_vect_reach[i]){
      while ((mem_vect_reach_hp[i] > mem_vect_reach[i]) | (mem_vect_reach_hp[i]  <= vect_reach_hp[i])) {
        mem_vect_reach_hp[i] = max(0,rnorm(1, mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor))
      }
    }
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}


################################################################################

# Unilateral truncation corrected + double truncation #

vf_utc_dt = function(survey_hp,Population, Mhp_vis,memory_factor){
  
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double truncation + double truncation calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    if (mem_vect_reach[i]   <= vect_reach_hp[i]){
      while ((mem_vect_reach[i]   <= vect_reach_hp[i]) | (mem_vect_reach[i] > vect_reach[i])) {
        mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
      }
    }
    mem_vect_reach_hp[i] = max(0,round(rtruncnorm(1, a = vect_reach_hp[i] - (mem_vect_reach[i]-vect_reach_hp[i]), b = mem_vect_reach[i], mean = vect_reach_hp[i] , sd = vect_reach_hp[i]*memory_factor)))
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}


################################################################################

# Unilateral truncation corrected + fake binomial #

vf_utc_fbin = function(survey_hp,Population, Mhp_vis,memory_factor){
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double trunc + binomial calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    if (mem_vect_reach[i]  <= vect_reach_hp[i]){
      while ((mem_vect_reach[i]  <= vect_reach_hp[i]) | (mem_vect_reach[i] > vect_reach[i])) {
        mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
      }
    }
    mem_vect_reach_hp[i] = rbinom(1, size = mem_vect_reach[i], p = vect_reach_hp[i]/mem_vect_reach[i])
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}


################################################################################

# Unilateral truncation corrected + binomial #

vf_utc_bin = function(survey_hp,Population, Mhp_vis,memory_factor){
  # Variables #
  ind_survey = as.numeric(rownames(survey_hp))
  vect_reach_hp = colSums(Mhp_vis[,ind_survey])
  vect_reach = Population$Reach[ind_survey]
  mem_vect_reach_hp = rep(NA,length(vect_reach_hp))
  mem_vect_reach = rep(NA,length(vect_reach_hp))
  
  # Double trunc + binomial calculate
  for (i in 1:length(vect_reach_hp)) {
    mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
    if (mem_vect_reach[i]  <= vect_reach_hp[i]){
      while ((mem_vect_reach[i]  <= vect_reach_hp[i]) | (mem_vect_reach[i] > vect_reach[i])) {
        mem_vect_reach[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
      }
    }
    mem_vect_reach_hp[i] = rbinom(1, size = mem_vect_reach[i], p = vect_reach_hp[i]/vect_reach[i])
  }
  
  vf_reach = sum(mem_vect_reach_hp)/sum(mem_vect_reach)
  return(vf_reach)
}


################################################################################


################################################################################
#####################
## Simulation data ##
#####################

N = 1000                      # Population size
v_pop = c(1:5)                # Subpopulations vector 
n_pop = length(v_pop)         # Number of subpopulations
v_pop_prob = rep(1/10, 5)     # Probability of each subpopulation. As we are working with disjoint and no disjoint subpopulations
                              # sum(v_pop_prob)  <= 1. 
hp_prob = 0.1                 # Probability for an individual to be in the hidden population (People who have COVID-19)
n_survey = 500                # Number of individuals we draw in the survey
n_survey_hp = 50              # Number of individuals we draw in the hidden population survey 

sub_memory_factor = 0         # Subpopulation memory factor (parameter to change variance of the perturbations' normal)
visibility_factor = 0.8         # Visibility factor (Binomial's probability)
memory_factor = 0.3       #Reach memory factor (parameter to change variance of the perturbations' normal)



seed = 207                    # Seed
set.seed(seed)

#Graph
dim = 1      # Graph dimension 
nei = 75     # Number of neighbours that each node is connected to. They are neighbors on each side of the node, so they are 2*nei connections
# before applying the randomization.
p   = 0.1    # Probability of randomize a connection. It is applied to all connections

Graph_population_matrix = getData(N, v_pop, v_pop_prob, hp_prob, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)

net_sw = Graph_population_matrix[[1]]       # Population´s graph
Population = Graph_population_matrix[[2]]   # Population
Mhp_vis = Graph_population_matrix[[3]]      # Population's visibility matrix

survey_hp = getSurvey(50, Population[Population$Hidden_Population == 1,])

################################################################################

# Estimates #

# Real value
visibility_factor

#Binomial aproach
vf_bin(survey_hp,Population, Mhp_vis, memory_factor)

# Double truncation
vf_dt_bin(survey_hp,Population, Mhp_vis, memory_factor)
vf_dt_fbin(survey_hp,Population, Mhp_vis, memory_factor) 
vf_dt_dt(survey_hp,Population, Mhp_vis, memory_factor) 
vf_dt_ut(survey_hp,Population, Mhp_vis, memory_factor) 
vf_dt_utc(survey_hp,Population, Mhp_vis, memory_factor) 


# Unilateral truncation 
vf_ut_bin(survey_hp,Population, Mhp_vis, memory_factor) 
vf_ut_fbin(survey_hp,Population, Mhp_vis, memory_factor)
vf_ut_dt(survey_hp,Population, Mhp_vis, memory_factor)
vf_ut_ut(survey_hp,Population, Mhp_vis, memory_factor)


# Unilateral truncation with correction
vf_utc_dt(survey_hp,Population, Mhp_vis, memory_factor) 
vf_utc_bin(survey_hp,Population, Mhp_vis, memory_factor) 
vf_utc_fbin(survey_hp,Population, Mhp_vis, memory_factor) 
vf_utc_utc(survey_hp,Population, Mhp_vis, memory_factor) 


