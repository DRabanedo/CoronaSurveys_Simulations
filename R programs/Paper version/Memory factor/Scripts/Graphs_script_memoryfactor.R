#################### Memory factor graphs scripts ##############################

# This script generates the graphs automatically as a png archive #


###### Packages ######

library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)

######################
# Data import
simulation_data = read.csv("C:/Users/David Rabanedo/Documents/Simulations_memoryfactor_207.csv")

seed_number = "207"
getwd()

################################################################################
################################################################################

## Not disjoint ##

Nh_real_dataframe = select(simulation_data, starts_with("Nh_real"))

Nh_basic_sum_dataframe = select(simulation_data, starts_with("Nh_basic_sum"))
#Nh_basicvis_sum_dataframe = select(simulation_data, starts_with("Nh_basicvis_sum"))

Nh_basic_mean_dataframe = select(simulation_data, starts_with("Nh_basic_mean"))
#Nh_basicvis_mean_dataframe = select(simulation_data, starts_with("Nh_basicvis_mean"))


########################### Data analysis ######################################

# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

Nh_basic_sum_analysis      = data_analysis(Nh_basic_sum_dataframe, Nh_real_dataframe)
#Nh_basicvis_sum_analysis  = data_analysis(Nh_basicvis_sum_dataframe, Nh_real_dataframe)

Nh_basic_mean_analysis     = data_analysis(Nh_basic_mean_dataframe, Nh_real_dataframe)
#Nh_basicvis_mean_analysis = data_analysis(Nh_basicvis_mean_dataframe, Nh_real_dataframe)


################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### Absolute error #####

# Dataframe creation #

graph_data_abserror = data.frame(data = simulation_data$data)

graph_data_abserror = cbind(graph_data_abserror, Nh_basic_sum =  Nh_basic_sum_analysis$abserror)
#graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$abserror)
graph_data_abserror = cbind(graph_data_abserror, Nh_basic_mean =  Nh_basic_mean_analysis$abserror)
# graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$abserror)



plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_abserror.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       subtitle = sub_title,
       x = "Memory factor",
       y = "Mean Absolute Error")

dev.off()


################################################################################


##### Mean of Squares Error (MSE) #####

# Dataframe creation #

graph_data_mse = data.frame( data = simulation_data$data)

graph_data_mse = cbind(graph_data_mse, Nh_basic_sum =  Nh_basic_sum_analysis$mse)
#graph_data_mse = cbind(graph_data_mse, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$mse)
graph_data_mse = cbind(graph_data_mse, Nh_basic_mean =  Nh_basic_mean_analysis$mse)
#graph_data_mse = cbind(graph_data_mse, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$mse)


# Graph creation

plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_mse.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_mse) + 
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation memory factor",
       subtitle = sub_title,
       x = "Memory factor",
       y = "Mean Squared Error (MSE)")


dev.off()

################################################################################

###### Bias analysis ######

# Dataframe creation #

graph_data_bias = data.frame( data = simulation_data$data)

graph_data_bias = cbind(graph_data_bias, Nh_real =  simulation_data$Nh_real_1)
graph_data_bias = cbind(graph_data_bias, Nh_basic_sum =  Nh_basic_sum_analysis$bias)
#graph_data_bias = cbind(graph_data_bias, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$bias)
graph_data_bias = cbind(graph_data_bias, Nh_basic_mean =  Nh_basic_mean_analysis$bias)
#graph_data_bias = cbind(graph_data_bias, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$bias)


# Graph creation

plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_bias.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       subtitle = sub_title,
       x = "Memory factor",
       y = "Hidden population estimate")


dev.off()



################################################################################

#### Standard deviation analysis ####

# Dataframe creation #

graph_data_sd = data.frame( data = simulation_data$data)

graph_data_sd = cbind(graph_data_sd, Nh_basic_sum =  Nh_basic_sum_analysis$sd)
#graph_data_sd = cbind(graph_data_sd, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$sd)
graph_data_sd = cbind(graph_data_sd, Nh_basic_mean =  Nh_basic_mean_analysis$sd)
#graph_data_sd = cbind(graph_data_sd, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$sd)


# Graph creation

plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_sd.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       subtitle = sub_title,
       x = "Memory factor",
       y = "Standard deviation")


dev.off()



################################################################################

#### Median analysis ####

# Dataframe creation #

graph_data_median = data.frame( data = simulation_data$data)

graph_data_median = cbind(graph_data_median, Nh_basic_sum =  Nh_basic_sum_analysis$median)
#graph_data_median = cbind(graph_data_median, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$median)
graph_data_median = cbind(graph_data_median, Nh_basic_mean =  Nh_basic_mean_analysis$median)
#graph_data_median = cbind(graph_data_median, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$median)


# Graph creation

plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_median.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_median) + 
  geom_line(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       subtitle = sub_title,
       x = "Memory factor",
       y = "Median")

dev.off()

