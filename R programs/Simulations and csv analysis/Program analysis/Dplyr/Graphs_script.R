library(dplyr)
library(matrixStats)
library(ggplot2)

simulation_data = read.csv("~/GitHub/CoronaSurveys_Simulations/R programs/Simulations and csv analysis/Csv archives/Visibility factor/Simulations_visibilityfactor_h")

#Column selection

Nh_real_dataframe = select(simulation_data, starts_with("Nh_real"))

Nh_basic_sum_dataframe = select(simulation_data, starts_with("Nh_basic_sum"))
Nh_basicvis_sum_dataframe = select(simulation_data, starts_with("Nh_basicvis_sum"))

Nh_basic_mean_dataframe = select(simulation_data, starts_with("Nh_basic_mean"))
Nh_basicvis_mean_dataframe = select(simulation_data, starts_with("Nh_basicvis_mean"))

Nh_PIMLE_dataframe    = select(simulation_data, starts_with("Nh_PIMLE_"))
Nh_PIMLEvis_dataframe = select(simulation_data, starts_with("Nh_PIMLEvis_"))

Nh_MLE_dataframe     = select(simulation_data, starts_with("Nh_MLE_"))
Nh_MLEvis_dataframe  = select(simulation_data, starts_with("Nh_MLEvis_"))

Nh_MoS_dataframe     = select(simulation_data, starts_with("Nh_MoS_"))
Nh_MoSvis_dataframe  = select(simulation_data, starts_with("Nh_MoSvis_"))

Nh_GNSUM_dataframe   = select(simulation_data, starts_with("Nh_GNSUM"))

Nh_Direct_dataframe = select(simulation_data, starts_with("Nh_Direct"))



######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

Nh_basic_sum_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basic_sum_dataframe-Nh_real_dataframe))),
                                   mse = rowMeans(as.matrix((Nh_basic_sum_dataframe-Nh_real_dataframe)^2)),
                                   bias = rowMeans(as.matrix(Nh_basic_sum_dataframe-Nh_real_dataframe)),
                                   sd = rowSds(as.matrix(Nh_basic_sum_dataframe))) 

Nh_basicvis_sum_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basicvis_sum_dataframe-Nh_real_dataframe))),
                                      mse = rowMeans(as.matrix((Nh_basicvis_sum_dataframe-Nh_real_dataframe)^2)),
                                      bias = rowMeans(as.matrix(Nh_basicvis_sum_dataframe-Nh_real_dataframe)),
                                      sd = rowSds(as.matrix(Nh_basicvis_sum_dataframe)))



Nh_basic_mean_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basic_mean_dataframe-Nh_real_dataframe))),
                                    mse = rowMeans(as.matrix((Nh_basic_mean_dataframe-Nh_real_dataframe)^2)),
                                    bias = rowMeans(as.matrix(Nh_basic_mean_dataframe-Nh_real_dataframe)),
                                    sd = rowSds(as.matrix(Nh_basic_mean_dataframe)))

Nh_basicvis_mean_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_basicvis_mean_dataframe-Nh_real_dataframe))),
                                       mse = rowMeans(as.matrix((Nh_basicvis_mean_dataframe-Nh_real_dataframe)^2)),
                                       bias = rowMeans(as.matrix(Nh_basicvis_mean_dataframe-Nh_real_dataframe)),
                                       sd = rowSds(as.matrix(Nh_basicvis_mean_dataframe)))



Nh_PIMLE_analysis    = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_PIMLE_dataframe-Nh_real_dataframe))),
                                  mse = rowMeans(as.matrix((Nh_PIMLE_dataframe-Nh_real_dataframe)^2)),
                                  bias = rowMeans(as.matrix(Nh_PIMLE_dataframe-Nh_real_dataframe)),
                                  sd = rowSds(as.matrix(Nh_PIMLE_dataframe)))

Nh_PIMLEvis_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_PIMLEvis_dataframe-Nh_real_dataframe))),
                                  mse = rowMeans(as.matrix((Nh_PIMLEvis_dataframe-Nh_real_dataframe)^2)),
                                  bias = rowMeans(as.matrix(Nh_PIMLEvis_dataframe-Nh_real_dataframe)),
                                  sd = rowSds(as.matrix(Nh_PIMLEvis_dataframe)))



Nh_MLE_analysis     = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_MLE_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MLE_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MLE_dataframe-Nh_real_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MLE_dataframe)))

Nh_MLEvis_analysis  = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_MLEvis_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MLEvis_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MLEvis_dataframe-Nh_real_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MLEvis_dataframe)))



Nh_MoS_analysis     = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_MoS_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MoS_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MoS_dataframe-Nh_real_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MoS_dataframe)))

Nh_MoSvis_analysis  = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_MoSvis_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MoSvis_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MoSvis_dataframe-Nh_real_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MoSvis_dataframe)))



Nh_GNSUM_analysis  = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_GNSUM_dataframe-Nh_real_dataframe))),
                                mse = rowMeans(as.matrix((Nh_GNSUM_dataframe-Nh_real_dataframe)^2)),
                                bias = rowMeans(as.matrix(Nh_GNSUM_dataframe-Nh_real_dataframe)),
                                sd = rowSds(as.matrix(Nh_GNSUM_dataframe)))



Nh_Direct_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_Direct_dataframe-Nh_real_dataframe))),
                                mse = rowMeans(as.matrix((Nh_Direct_dataframe-Nh_real_dataframe)^2)),
                                bias = rowMeans(as.matrix(Nh_Direct_dataframe-Nh_real_dataframe)),
                                sd = rowSds(as.matrix(Nh_Direct_dataframe)))



################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### Absolute error #####

#Dataframe creation

graph_data_abserror = data.frame( data = simulation_data$data)


if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basic_sum =  Nh_basic_sum_analysis$abs_error)
}

if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$abs_error)
}



if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basic_mean =  Nh_basic_mean_analysis$abs_error)
}

if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$abs_error)
}




if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_PIMLE =  Nh_PIMLE_analysis$abs_error)
}

if(ncol(Nh_PIMLEvis_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$abs_error)
}




if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MLE =  Nh_MLE_analysis$abs_error)
}

if(ncol(Nh_MLEvis_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MLEvis =  Nh_MLEvis_analysis$abs_error)
}




if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MoS =  Nh_MoS_analysis$abs_error)
}

if(ncol(Nh_MoSvis_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MoSvis =  Nh_MoSvis_analysis$abs_error)
}



if(ncol(Nh_GNSUM_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_GNSUM  =  Nh_GNSUM_analysis$abs_error)
}

if(ncol(Nh_Direct_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_Direct  =  Nh_Direct_analysis$abs_error)
}



ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_line(aes(x = data, y =  Nh_Direct, col = "Nh_Direct")) +
  scale_color_discrete("Legend") + 
  labs(title = "",
       x = "",
       y = "Mean Absolute Error (MAE")


################################################################################


##### Mean of Squares Error (MSE) #####

#Dataframe creation

graph_data_mse = data.frame( data = simulation_data$data)


if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_basic_sum =  Nh_basic_sum_analysis$mse)
}

if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$mse)
}



if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_basic_mean =  Nh_basic_mean_analysis$mse)
}

if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$mse)
}




if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_PIMLE =  Nh_PIMLE_analysis$mse)
}

if(ncol(Nh_PIMLEvis_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$mse)
}




if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MLE =  Nh_MLE_analysis$mse)
}

if(ncol(Nh_MLEvis_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MLEvis =  Nh_MLEvis_analysis$mse)
}




if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MoS =  Nh_MoS_analysis$mse)
}

if(ncol(Nh_MoSvis_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MoSvis =  Nh_MoSvis_analysis$mse)
}



if(ncol(Nh_GNSUM_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_GNSUM  =  Nh_GNSUM_analysis$mse)
}

if(ncol(Nh_Direct_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_Direct  =  Nh_Direct_analysis$mse)
}



ggplot(graph_data_mse) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_line(aes(x = data, y =  Nh_Direct, col = "Nh_Direct")) +
  scale_color_discrete("Legend") + 
  labs(title = "",
       x = "",
       y = "Mean Squared Error (MSE)")


################################################################################

###### Bias analysis ######


graph_data_bias = data.frame( data = simulation_data$data)


if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_basic_sum =  Nh_basic_sum_analysis$bias)
}

if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$bias)
}



if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_basic_mean =  Nh_basic_mean_analysis$bias)
}

if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$bias)
}




if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_PIMLE =  Nh_PIMLE_analysis$bias)
}

if(ncol(Nh_PIMLEvis_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$bias)
}




if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MLE =  Nh_MLE_analysis$bias)
}

if(ncol(Nh_MLEvis_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MLEvis =  Nh_MLEvis_analysis$bias)
}




if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MoS =  Nh_MoS_analysis$bias)
}

if(ncol(Nh_MoSvis_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MoSvis =  Nh_MoSvis_analysis$bias)
}



if(ncol(Nh_GNSUM_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_GNSUM  =  Nh_GNSUM_analysis$bias)
}

if(ncol(Nh_Direct_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_Direct  =  Nh_Direct_analysis$bias)
}



ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_line(aes(x = data, y =  Nh_Direct, col = "Nh_Direct")) +
  scale_color_discrete("Legend") + 
  labs(title = "",
       x = "",
       y = "Hidden population estimate")



################################################################################





