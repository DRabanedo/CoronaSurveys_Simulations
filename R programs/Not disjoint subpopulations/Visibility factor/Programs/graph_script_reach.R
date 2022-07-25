library(dplyr)
library(matrixStats)
library(ggplot2)

simulation_data = read.csv("~/GitHub/CoronaSurveys_Simulations/R programs/Disjoint subpopulations/Visibility factor/Analysis/Csv archives/seed_207/Simulations_visibilityfactor_estimate_reach")

######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.
simulation_data_v =  select(simulation_data, starts_with("vf_reach"))

Visibility_factor_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(simulation_data_v-1))),
                                        mse = rowMeans(as.matrix((simulation_data_v-1)^2)),
                                        bias = rowMeans(as.matrix(simulation_data_v)),
                                        sd = rowSds(as.matrix(simulation_data_v)))
################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### MSE #####

#Dataframe creation

graph_data_mse = data.frame(data = simulation_data$data, visibility_factor = Visibility_factor_analysis$mse)



ggplot(graph_data_mse) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate")) + 
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       x = "Memory factor",
       y = "Mean Squared Error (MSE)")


################################################################################

###### Bias analysis ######


graph_data_bias = data.frame(data = simulation_data$data, visibility_factor = Visibility_factor_analysis$bias)



ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate")) + 
  geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       x = "Memory factor",
       y = "Visibility factor value")



################################################################################

#### Standard deviation analysis ####

#Dataframe creation


graph_data_sd = data.frame(data = simulation_data$data, visibility_factor = Visibility_factor_analysis$sd)



ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate")) + 
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       x = "Memory factor",
       y = "Standard deviation")

#################################################################################

#### Absolute error ####


graph_data_abserror = data.frame(data = simulation_data$data, visibility_factor = Visibility_factor_analysis$abs_error)



ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate")) + 
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       x = "Memory factor",
       y = "Absolute error")

