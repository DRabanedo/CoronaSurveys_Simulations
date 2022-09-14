library(dplyr)
library(matrixStats)
library(ggplot2)

simulation_data = read.csv("~/Simulations_memoryfactor_disjoint")

Nh_real_dataframe = select(simulation_data, starts_with("Nh_real"))

Nh_basic_sum_dataframe = select(simulation_data, starts_with("Nh_basic_sum"))
#Nh_basicvis_sum_dataframe = select(simulation_data, starts_with("Nh_basicvis_sum"))

Nh_basic_mean_dataframe = select(simulation_data, starts_with("Nh_basic_mean"))
#Nh_basicvis_mean_dataframe = select(simulation_data, starts_with("Nh_basicvis_mean"))


######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

Nh_basic_sum_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basic_sum_dataframe-Nh_real_dataframe))),
                                   mse = rowMeans(as.matrix((Nh_basic_sum_dataframe-Nh_real_dataframe)^2)),
                                   bias = rowMeans(as.matrix(Nh_basic_sum_dataframe)),
                                   sd = rowSds(as.matrix(Nh_basic_sum_dataframe))) 

#Nh_basicvis_sum_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basicvis_sum_dataframe-Nh_real_dataframe))),
                                      #mse = rowMeans(as.matrix((Nh_basicvis_sum_dataframe-Nh_real_dataframe)^2)),
                                      #bias = rowMeans(as.matrix(Nh_basicvis_sum_dataframe)),
                                      #sd = rowSds(as.matrix(Nh_basicvis_sum_dataframe)))



Nh_basic_mean_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basic_mean_dataframe-Nh_real_dataframe))),
                                    mse = rowMeans(as.matrix((Nh_basic_mean_dataframe-Nh_real_dataframe)^2)),
                                    bias = rowMeans(as.matrix(Nh_basic_mean_dataframe)),
                                    sd = rowSds(as.matrix(Nh_basic_mean_dataframe)))

#Nh_basicvis_mean_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basicvis_mean_dataframe-Nh_real_dataframe))),
                                       #mse = rowMeans(as.matrix((Nh_basicvis_mean_dataframe-Nh_real_dataframe)^2)),
                                       #bias = rowMeans(as.matrix(Nh_basicvis_mean_dataframe)),
                                       #sd = rowSds(as.matrix(Nh_basicvis_mean_dataframe)))


################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### Absolute error #####

#Dataframe creation

graph_data_abserror = data.frame(data = simulation_data$data)

if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basic_sum =  Nh_basic_sum_analysis$abs_error)
}

#if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
#  graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$abs_error)
#}



if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basic_mean =  Nh_basic_mean_analysis$abs_error)
}

#if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
# # graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$abs_error)
#}




ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       x = "Memory factor",
       y = "Mean Absolute Error")


################################################################################


##### Mean of Squares Error (MSE) #####

#Dataframe creation

graph_data_mse = data.frame( data = simulation_data$data)


if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_basic_sum =  Nh_basic_sum_analysis$mse)
}

#if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
#  graph_data_mse = cbind(graph_data_mse, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$mse)
#}



if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_basic_mean =  Nh_basic_mean_analysis$mse)
}

#if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
#  graph_data_mse = cbind(graph_data_mse, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$mse)
#}




ggplot(graph_data_mse) + 
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation memory factor",
       x = "Memory factor",
       y = "Mean Squared Error (MSE)")


################################################################################

###### Bias analysis ######


graph_data_bias = data.frame( data = simulation_data$data)

graph_data_bias = cbind(graph_data_bias, Nh_real =  simulation_data$Nh_real_1)

if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_basic_sum =  Nh_basic_sum_analysis$bias)
}

#if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
#  graph_data_bias = cbind(graph_data_bias, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$bias)
#}



if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_basic_mean =  Nh_basic_mean_analysis$bias)
}

#if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
#  graph_data_bias = cbind(graph_data_bias, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$bias)
#}


ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       x = "Memory factor",
       y = "Hidden population estimate")



################################################################################

#### Standard deviation analysis ####

#Dataframe creation

graph_data_sd = data.frame( data = simulation_data$data)


if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_basic_sum =  Nh_basic_sum_analysis$sd)
}

#if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
#  graph_data_sd = cbind(graph_data_sd, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$sd)
#}



if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_basic_mean =  Nh_basic_mean_analysis$sd)
}

#if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
#  graph_data_sd = cbind(graph_data_sd, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$sd)
#}

ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       x = "Memory factor",
       y = "Standard deviation")
