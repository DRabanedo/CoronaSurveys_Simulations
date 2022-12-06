#################### Network probability graphs scripts ########################

# This script generates the graphs automatically as a png archive #

###### Packages ######

library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)

######################
# Data import

simulation_data = read.csv("~/GitHub/CoronaSurveys_Simulations/R programs/Paper version/Subpopulation size/CSV/Simulations_subpopulationsize_notdisjoint_2023.csv")
simulation_data_disjoint = read.csv("~/GitHub/CoronaSurveys_Simulations/R programs/Paper version/Subpopulation size/CSV/Simulations_subpopulationsize_disjoint_2023.csv")
getwd()
seed_number = 2023


################################################################################
################################################################################

## Not disjoint ##

Nh_real_dataframe = select(simulation_data, starts_with("Nh_real"))

Nh_PIMLE_dataframe    = select(simulation_data, starts_with("Nh_PIMLE_"))
#Nh_PIMLEvis_dataframe = select(simulation_data, starts_with("Nh_PIMLEvis_"))

Nh_MLE_dataframe     = select(simulation_data, starts_with("Nh_MLE_"))
#Nh_MLEvis_dataframe  = select(simulation_data, starts_with("Nh_MLEvis_"))

Nh_MoS_dataframe     = select(simulation_data, starts_with("Nh_MoS_"))
#Nh_MoSvis_dataframe  = select(simulation_data, starts_with("Nh_MoSvis_"))

Nh_GNSUM_dataframe   = select(simulation_data, starts_with("Nh_GNSUM"))



######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.


Nh_PIMLE_analysis     = data_analysis(Nh_PIMLE_dataframe, Nh_real_dataframe)
#Nh_PIMLEvis_analysis = data_analysis(Nh_PIMLEvis_dataframe, Nh_real_dataframe)

Nh_MLE_analysis     = data_analysis(Nh_MLE_dataframe, Nh_real_dataframe)
#Nh_MLEvis_analysis = data_analysis(Nh_MLEvis_dataframe, Nh_real_dataframe)

Nh_MoS_analysis     = data_analysis(Nh_MoS_dataframe, Nh_real_dataframe)
#Nh_MoSvis_analysis = data_analysis(Nh_MoSvis_dataframe, Nh_real_dataframe)

Nh_GNSUM_analysis  = data_analysis(Nh_GNSUM_dataframe, Nh_real_dataframe)


################################################################################
################################################################################

## Disjoint ##

Nh_real_dataframe_disjoint = select(simulation_data_disjoint, starts_with("Nh_real"))

Nh_PIMLE_dataframe_disjoint    = select(simulation_data_disjoint, starts_with("Nh_PIMLE_"))
#Nh_PIMLEvis_dataframe_disjoint = select(simulation_data_disjoint, starts_with("Nh_PIMLEvis_"))

Nh_MLE_dataframe_disjoint     = select(simulation_data_disjoint, starts_with("Nh_MLE_"))
#Nh_MLEvis_dataframe_disjoint  = select(simulation_data_disjoint, starts_with("Nh_MLEvis_"))

Nh_MoS_dataframe_disjoint     = select(simulation_data_disjoint, starts_with("Nh_MoS_"))
#Nh_MoSvis_dataframe_disjoint  = select(simulation_data_disjoint, starts_with("Nh_MoSvis_"))

Nh_GNSUM_dataframe_disjoint  = select(simulation_data_disjoint, starts_with("Nh_GNSUM"))


######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.


Nh_PIMLE_analysis_disjoint     = data_analysis(Nh_PIMLE_dataframe_disjoint, Nh_real_dataframe)
#Nh_PIMLEvis_analysis_disjoint = data_analysis(Nh_PIMLEvis_dataframe_disjoint, Nh_real_dataframe)

Nh_MLE_analysis_disjoint     = data_analysis(Nh_MLE_dataframe_disjoint, Nh_real_dataframe)
#Nh_MLEvis_analysis_disjoint = data_analysis(Nh_MLEvis_dataframe_disjoint, Nh_real_dataframe)

Nh_MoS_analysis_disjoint     = data_analysis(Nh_MoS_dataframe_disjoint, Nh_real_dataframe)
#Nh_MoSvis_analysis_disjoint = data_analysis(Nh_MoSvis_dataframe_disjoint, Nh_real_dataframe)

Nh_GNSUM_analysis_disjoint   = data_analysis(Nh_GNSUM_dataframe_disjoint, Nh_real_dataframe)




################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### Absolute error #####

## Not disjoint population ##

#Dataframe creation

graph_data_abserror = data.frame( data = simulation_data$data)

graph_data_abserror = cbind(graph_data_abserror, Nh_PIMLE =  Nh_PIMLE_analysis$abserror)
#graph_data_abserror = cbind(graph_data_abserror, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$abserror)

graph_data_abserror = cbind(graph_data_abserror, Nh_MLE =  Nh_MLE_analysis$abserror)
#graph_data_abserror = cbind(graph_data_abserror, Nh_MLEvis =  Nh_MLEvis_analysis$abserror)

graph_data_abserror = cbind(graph_data_abserror, Nh_MoS =  Nh_MoS_analysis$abserror)
#graph_data_abserror = cbind(graph_data_abserror, Nh_MoSvis =  Nh_MoSvis_analysis$abserror)

graph_data_abserror = cbind(graph_data_abserror, Nh_GNSUM  =  Nh_GNSUM_analysis$abserror)


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_notdisjoint_abserror.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_abserror) +
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Mean Absolute Error")

dev.off()



## Disjoint populations ##

#Dataframe creation

graph_data_abserror_disjoint = data.frame( data = simulation_data_disjoint$data)

graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$abserror)
#graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$abserror)

graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$abserror)
#graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$abserror)

graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$abserror)
#graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$abserror)

graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$abserror)


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_disjoint_abserror.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_abserror_disjoint) +
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Mean Absolute Error")

dev.off()



## Not disjoint & disjoint ##

# Dataframe creation

graph_data_abserror_total = cbind(graph_data_abserror, graph_data_abserror_disjoint[2:ncol(graph_data_abserror_disjoint)])

# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_total_abserror.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_abserror_total) +
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Mean Absolute Error")

dev.off()

################################################################################

##### Mean of Squares Error (MSE) #####

## Not disjoint populations ##

#Dataframe creation

graph_data_mse = data.frame(data = simulation_data$data)

graph_data_mse = cbind(graph_data_mse, Nh_PIMLE =  Nh_PIMLE_analysis$mse)
#graph_data_mse = cbind(graph_data_mse, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$mse)

graph_data_mse = cbind(graph_data_mse, Nh_MLE =  Nh_MLE_analysis$mse)
#graph_data_mse = cbind(graph_data_mse, Nh_MLEvis =  Nh_MLEvis_analysis$mse)

graph_data_mse = cbind(graph_data_mse, Nh_MoS =  Nh_MoS_analysis$mse)
#graph_data_mse = cbind(graph_data_mse, Nh_MoSvis =  Nh_MoSvis_analysis$mse)

graph_data_mse = cbind(graph_data_mse, Nh_GNSUM  =  Nh_GNSUM_analysis$mse)


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_notdisjoint_mse.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_mse) + 
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Mean Squared Error (MSE)")

dev.off()


## Disjoint population ##

# Dataframe creation

graph_data_mse_disjoint = data.frame(data = simulation_data_disjoint$data)

graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$mse)
#graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$mse)

graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$mse)
#graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$mse)

graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$mse)
#  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$mse)

graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$mse)


# Graph creation 
plot_name = str_c("Simulation_subpopsize_", seed_number, "_disjoint_mse.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_mse_disjoint) + 
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Mean Squared Error (MSE)")

dev.off()

## Not disjoint & disjoint ##

#Dataframe creation

graph_data_mse_total = cbind(graph_data_mse, graph_data_mse_disjoint[2:ncol(graph_data_mse_disjoint)])


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_total_mse.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_mse_total) +
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Mean Squared Error (MSE)")

dev.off()



################################################################################

###### Bias analysis ######

## Not disjoint ##

# Dataframe creation

graph_data_bias = data.frame( data = simulation_data$data)

graph_data_bias = cbind(graph_data_bias, Nh_real =  simulation_data$Nh_real_1)

graph_data_bias = cbind(graph_data_bias, Nh_PIMLE =  Nh_PIMLE_analysis$bias)
#graph_data_bias = cbind(graph_data_bias, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$bias)

graph_data_bias = cbind(graph_data_bias, Nh_MLE =  Nh_MLE_analysis$bias)
#graph_data_bias = cbind(graph_data_bias, Nh_MLEvis =  Nh_MLEvis_analysis$bias)

graph_data_bias = cbind(graph_data_bias, Nh_MoS =  Nh_MoS_analysis$bias)
#graph_data_bias = cbind(graph_data_bias, Nh_MoSvis =  Nh_MoSvis_analysis$bias)

graph_data_bias = cbind(graph_data_bias, Nh_GNSUM  =  Nh_GNSUM_analysis$bias)


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_notdisjoint_bias.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_bias) + 
  geom_point(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Hidden population estimate")

dev.off()


## Disjoint population ##

# Dataframe creation

graph_data_bias_disjoint = data.frame(data = simulation_data_disjoint$data)

graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_real =  simulation_data_disjoint$Nh_real_1)

graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$bias)
#graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$bias)

graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$bias)
#graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$bias)

graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$bias)
#graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$bias)

graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$bias)


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_disjoint_bias.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_bias_disjoint) + 
  geom_point(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Hidden population estimate")

dev.off()

## Not disjoint & disjoint ##

# Dataframe creation

graph_data_bias_total = cbind(graph_data_bias, graph_data_bias_disjoint[3:ncol(graph_data_bias_disjoint)])


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_total_bias.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_bias_total) +
  geom_point(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Hidden population estimate")

dev.off()



################################################################################

#### Standard deviation analysis ####

#Dataframe creation

graph_data_sd = data.frame( data = simulation_data$data)

graph_data_sd = cbind(graph_data_sd, Nh_PIMLE =  Nh_PIMLE_analysis$sd)
#graph_data_sd = cbind(graph_data_sd, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$sd)

graph_data_sd = cbind(graph_data_sd, Nh_MLE =  Nh_MLE_analysis$sd)
#graph_data_sd = cbind(graph_data_sd, Nh_MLEvis =  Nh_MLEvis_analysis$sd)

graph_data_sd = cbind(graph_data_sd, Nh_MoS =  Nh_MoS_analysis$sd)
#graph_data_sd = cbind(graph_data_sd, Nh_MoSvis =  Nh_MoSvis_analysis$sd)

graph_data_sd = cbind(graph_data_sd, Nh_GNSUM  =  Nh_GNSUM_analysis$sd)


# Data creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_notdisjoint_sd.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)


ggplot(graph_data_sd) + 
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Standard deviation")

dev.off()



## Disjoint population ##


graph_data_sd_disjoint = data.frame(data = simulation_data_disjoint$data)

graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$sd)
#graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$sd)

graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$sd)
#graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$sd)

graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$sd)
#graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$sd)

graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$sd)


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_disjoint_sd.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_sd_disjoint) + 
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Standard deviation")

dev.off()

## Not disjoint & disjoint ##

# Dataframe creation

graph_data_sd_total = cbind(graph_data_sd, graph_data_sd_disjoint[2:ncol(graph_data_sd_disjoint)])


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_total_sd.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_sd_total) +
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Standard error")

dev.off()

################################################################################

#### Median analysis ####

#Dataframe creation

graph_data_median = data.frame( data = simulation_data$data)

graph_data_median = cbind(graph_data_median, Nh_PIMLE =  Nh_PIMLE_analysis$median)
#graph_data_median = cbind(graph_data_median, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$median)

graph_data_median = cbind(graph_data_median, Nh_MLE =  Nh_MLE_analysis$median)
#graph_data_median = cbind(graph_data_median, Nh_MLEvis =  Nh_MLEvis_analysis$median)

graph_data_median = cbind(graph_data_median, Nh_MoS =  Nh_MoS_analysis$median)
#graph_data_median = cbind(graph_data_median, Nh_MoSvis =  Nh_MoSvis_analysis$median)

graph_data_median = cbind(graph_data_median, Nh_GNSUM  =  Nh_GNSUM_analysis$median)


# Data creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_notdisjoint_median.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)


ggplot(graph_data_median) + 
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Standard deviation")

dev.off()



## Disjoint population ##

# Dtaframe creation

graph_data_median_disjoint = data.frame(data = simulation_data_disjoint$data)

graph_data_median_disjoint = cbind(graph_data_median_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$median)
#graph_data_median_disjoint = cbind(graph_data_median_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$median)

graph_data_median_disjoint = cbind(graph_data_median_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$median)
#graph_data_median_disjoint = cbind(graph_data_median_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$median)

graph_data_median_disjoint = cbind(graph_data_median_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$median)
#graph_data_median_disjoint = cbind(graph_data_median_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$median)

graph_data_median_disjoint = cbind(graph_data_median_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$median)


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_disjoint_median.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_median_disjoint) + 
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Standard deviation")

dev.off()

## Not disjoint & disjoint ##

# Dataframe creation

graph_data_median_total = cbind(graph_data_median, graph_data_median_disjoint[2:ncol(graph_data_median_disjoint)])


# Graph creation

plot_name = str_c("Simulation_subpopsize_", seed_number, "_total_median.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_median_total) +
  #geom_point(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_point(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_point(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_point(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_point(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       subtitle = sub_title,
       x = "Subpopulation size",
       y = "Standard error")

dev.off()