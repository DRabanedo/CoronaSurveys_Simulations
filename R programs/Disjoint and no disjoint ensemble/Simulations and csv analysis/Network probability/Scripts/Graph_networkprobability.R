library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)

simulation_data = 
simulation_data_disjoint = 

seed_number = "207"
getwd()

simulation_data$data = simulation_data$data*2
simulation_data_disjoint$data = simulation_data_disjoint$data*2

##################
## Not disjoint ##

Nh_real_dataframe = dplyr::select(simulation_data, starts_with("Nh_real"))

Nh_basic_sum_dataframe = dplyr::select(simulation_data, starts_with("Nh_basic_sum"))
#Nh_basicvis_sum_dataframe = dplyr::select(simulation_data, starts_with("Nh_basicvis_sum"))

Nh_basic_mean_dataframe = dplyr::select(simulation_data, starts_with("Nh_basic_mean"))
#Nh_basicvis_mean_dataframe = dplyr::select(simulation_data, starts_with("Nh_basicvis_mean"))

Nh_PIMLE_dataframe    = dplyr::select(simulation_data, starts_with("Nh_PIMLE_"))
#Nh_PIMLEvis_dataframe = dplyr::select(simulation_data, starts_with("Nh_PIMLEvis_"))

Nh_MLE_dataframe     = dplyr::select(simulation_data, starts_with("Nh_MLE_"))
#Nh_MLEvis_dataframe  = dplyr::select(simulation_data, starts_with("Nh_MLEvis_"))

Nh_MoS_dataframe     = dplyr::select(simulation_data, starts_with("Nh_MoS_"))
#Nh_MoSvis_dataframe  = dplyr::select(simulation_data, starts_with("Nh_MoSvis_"))

Nh_GNSUM_dataframe   = dplyr::select(simulation_data, starts_with("Nh_GNSUM"))


######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

Nh_basic_sum_analysis = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basic_sum_dataframe-Nh_real_dataframe))),
                                   mse = rowMeans(as.matrix((Nh_basic_sum_dataframe-Nh_real_dataframe)^2)),
                                   bias = rowMeans(as.matrix(Nh_basic_sum_dataframe)),
                                   sd = rowSds(as.matrix(Nh_basic_sum_dataframe))) 

#Nh_basicvis_sum_analysis = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basicvis_sum_dataframe-Nh_real_dataframe))),
#mse = rowMeans(as.matrix((Nh_basicvis_sum_dataframe-Nh_real_dataframe)^2)),
#bias = rowMeans(as.matrix(Nh_basicvis_sum_dataframe)),
#sd = rowSds(as.matrix(Nh_basicvis_sum_dataframe)))



Nh_basic_mean_analysis = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basic_mean_dataframe-Nh_real_dataframe))),
                                    mse = rowMeans(as.matrix((Nh_basic_mean_dataframe-Nh_real_dataframe)^2)),
                                    bias = rowMeans(as.matrix(Nh_basic_mean_dataframe)),
                                    sd = rowSds(as.matrix(Nh_basic_mean_dataframe)))

#Nh_basicvis_mean_analysis = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basicvis_mean_dataframe-Nh_real_dataframe))),
#mse = rowMeans(as.matrix((Nh_basicvis_mean_dataframe-Nh_real_dataframe)^2)),
#bias = rowMeans(as.matrix(Nh_basicvis_mean_dataframe)),
#sd = rowSds(as.matrix(Nh_basicvis_mean_dataframe)))

Nh_PIMLE_analysis    = data.frame(abserror = rowMeans(as.matrix(abs(Nh_PIMLE_dataframe-Nh_real_dataframe))),
                                  mse = rowMeans(as.matrix((Nh_PIMLE_dataframe-Nh_real_dataframe)^2)),
                                  bias = rowMeans(as.matrix(Nh_PIMLE_dataframe)),
                                  sd = rowSds(as.matrix(Nh_PIMLE_dataframe)))

#Nh_PIMLEvis_analysis = data.frame(abserror = rowMeans(as.matrix(abs(Nh_PIMLEvis_dataframe-Nh_real_dataframe))),
#mse = rowMeans(as.matrix((Nh_PIMLEvis_dataframe-Nh_real_dataframe)^2)),
# bias = rowMeans(as.matrix(Nh_PIMLEvis_dataframe)),
#sd = rowSds(as.matrix(Nh_PIMLEvis_dataframe)))



Nh_MLE_analysis     = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MLE_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MLE_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MLE_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MLE_dataframe)))

#Nh_MLEvis_analysis  = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MLEvis_dataframe-Nh_real_dataframe))),
#mse = rowMeans(as.matrix((Nh_MLEvis_dataframe-Nh_real_dataframe)^2)),
#bias = rowMeans(as.matrix(Nh_MLEvis_dataframe)),
#sd = rowSds(as.matrix(Nh_MLEvis_dataframe)))



Nh_MoS_analysis     = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MoS_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MoS_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MoS_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MoS_dataframe)))

#Nh_MoSvis_analysis  = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MoSvis_dataframe-Nh_real_dataframe))),
#mse = rowMeans(as.matrix((Nh_MoSvis_dataframe-Nh_real_dataframe)^2)),
#bias = rowMeans(as.matrix(Nh_MoSvis_dataframe)),
#sd = rowSds(as.matrix(Nh_MoSvis_dataframe)))



Nh_GNSUM_analysis  = data.frame(abserror = rowMeans(as.matrix(abs(Nh_GNSUM_dataframe-Nh_real_dataframe))),
                                mse = rowMeans(as.matrix((Nh_GNSUM_dataframe-Nh_real_dataframe)^2)),
                                bias = rowMeans(as.matrix(Nh_GNSUM_dataframe)),
                                sd = rowSds(as.matrix(Nh_GNSUM_dataframe)))


##############
## Disjoint ##


Nh_real_dataframe_disjoint = dplyr::select(simulation_data_disjoint, starts_with("Nh_real"))

Nh_basic_sum_dataframe_disjoint = dplyr::select(simulation_data_disjoint, starts_with("Nh_basic_sum"))
#Nh_basicvis_sum_dataframe_disjoint = dplyr::select(simulation_data_disjoint, starts_with("Nh_basicvis_sum"))

Nh_basic_mean_dataframe_disjoint = dplyr::select(simulation_data_disjoint, starts_with("Nh_basic_mean"))
#Nh_basicvis_mean_dataframe_disjoint = dplyr::select(simulation_data_disjoint, starts_with("Nh_basicvis_mean"))

Nh_PIMLE_dataframe_disjoint    = dplyr::select(simulation_data_disjoint, starts_with("Nh_PIMLE_"))
#Nh_PIMLEvis_dataframe_disjoint = dplyr::select(simulation_data_disjoint, starts_with("Nh_PIMLEvis_"))

Nh_MLE_dataframe_disjoint     = dplyr::select(simulation_data_disjoint, starts_with("Nh_MLE_"))
#Nh_MLEvis_dataframe_disjoint  = dplyr::select(simulation_data_disjoint, starts_with("Nh_MLEvis_"))

Nh_MoS_dataframe_disjoint     = dplyr::select(simulation_data_disjoint, starts_with("Nh_MoS_"))
#Nh_MoSvis_dataframe_disjoint  = dplyr::select(simulation_data_disjoint, starts_with("Nh_MoSvis_"))

Nh_GNSUM_dataframe_disjoint  = dplyr::select(simulation_data_disjoint, starts_with("Nh_GNSUM"))


######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

Nh_basic_sum_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basic_sum_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                            mse = rowMeans(as.matrix((Nh_basic_sum_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                            bias = rowMeans(as.matrix(Nh_basic_sum_dataframe_disjoint)),
                                            sd = rowSds(as.matrix(Nh_basic_sum_dataframe_disjoint))) 

#Nh_basicvis_sum_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basicvis_sum_dataframe_disjoint-Nh_real_dataframe_disjoint))),
#                                               mse = rowMeans(as.matrix((Nh_basicvis_sum_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
#                                               bias = rowMeans(as.matrix(Nh_basicvis_sum_dataframe_disjoint)),
#                                               sd = rowSds(as.matrix(Nh_basicvis_sum_dataframe_disjoint)))



Nh_basic_mean_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basic_mean_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                             mse = rowMeans(as.matrix((Nh_basic_mean_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                             bias = rowMeans(as.matrix(Nh_basic_mean_dataframe_disjoint)),
                                             sd = rowSds(as.matrix(Nh_basic_mean_dataframe_disjoint)))

#Nh_basicvis_mean_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basicvis_mean_dataframe_disjoint-Nh_real_dataframe_disjoint))),
#mse = rowMeans(as.matrix((Nh_basicvis_mean_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
#bias = rowMeans(as.matrix(Nh_basicvis_mean_dataframe_disjoint)),
#sd = rowSds(as.matrix(Nh_basicvis_mean_dataframe_disjoint)))



Nh_PIMLE_analysis_disjoint    = data.frame(abserror = rowMeans(as.matrix(abs(Nh_PIMLE_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                           mse = rowMeans(as.matrix((Nh_PIMLE_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                           bias = rowMeans(as.matrix(Nh_PIMLE_dataframe_disjoint)),
                                           sd = rowSds(as.matrix(Nh_PIMLE_dataframe_disjoint)))

#Nh_PIMLEvis_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_PIMLEvis_dataframe_disjoint-Nh_real_dataframe_disjoint))),
#mse = rowMeans(as.matrix((Nh_PIMLEvis_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
# bias = rowMeans(as.matrix(Nh_PIMLEvis_dataframe_disjoint)),
#sd = rowSds(as.matrix(Nh_PIMLEvis_dataframe_disjoint)))



Nh_MLE_analysis_disjoint     = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MLE_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                          mse = rowMeans(as.matrix((Nh_MLE_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                          bias = rowMeans(as.matrix(Nh_MLE_dataframe_disjoint)),
                                          sd = rowSds(as.matrix(Nh_MLE_dataframe_disjoint)))

#Nh_MLEvis_analysis_disjoint  = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MLEvis_dataframe_disjoint-Nh_real_dataframe_disjoint))),
#mse = rowMeans(as.matrix((Nh_MLEvis_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
#bias = rowMeans(as.matrix(Nh_MLEvis_dataframe_disjoint)),
#sd = rowSds(as.matrix(Nh_MLEvis_dataframe_disjoint)))



Nh_MoS_analysis_disjoint     = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MoS_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                          mse = rowMeans(as.matrix((Nh_MoS_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                          bias = rowMeans(as.matrix(Nh_MoS_dataframe_disjoint)),
                                          sd = rowSds(as.matrix(Nh_MoS_dataframe_disjoint)))

#Nh_MoSvis_analysis_disjoint  = data.frame(abserror = rowMeans(as.matrix(abs(Nh_MoSvis_dataframe_disjoint-Nh_real_dataframe_disjoint))),
#mse = rowMeans(as.matrix((Nh_MoSvis_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
#bias = rowMeans(as.matrix(Nh_MoSvis_dataframe_disjoint)),
#sd = rowSds(as.matrix(Nh_MoSvis_dataframe_disjoint)))



Nh_GNSUM_analysis_disjoint  = data.frame(abserror = rowMeans(as.matrix(abs(Nh_GNSUM_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                         mse = rowMeans(as.matrix((Nh_GNSUM_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                         bias = rowMeans(as.matrix(Nh_GNSUM_dataframe_disjoint)),
                                         sd = rowSds(as.matrix(Nh_GNSUM_dataframe_disjoint)))



################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### Absolute error #####

## Not disjoint ##

#Dataframe creation

graph_data_abserror = data.frame(data = simulation_data$data)

if(ncol(Nh_basic_sum_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basic_sum =  Nh_basic_sum_analysis$abserror)
}

#if(ncol(Nh_basicvis_sum_dataframe) !=  0) {
#  graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_sum =  Nh_basicvis_sum_analysis$abserror)
#}


if(ncol(Nh_basic_mean_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_basic_mean =  Nh_basic_mean_analysis$abserror)
}

#if(ncol(Nh_basicvis_mean_dataframe) !=  0) {
# # graph_data_abserror = cbind(graph_data_abserror, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$abserror)
#}

if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_PIMLE =  Nh_PIMLE_analysis$abserror)
}

#if(ncol(Nh_PIMLEvis_dataframe) !=  0) {
#  graph_data_abserror = cbind(graph_data_abserror, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$abserror)
#}


if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MLE =  Nh_MLE_analysis$abserror)
}

#if(ncol(Nh_MLEvis_dataframe) !=  0) {
#  graph_data_abserror = cbind(graph_data_abserror, Nh_MLEvis =  Nh_MLEvis_analysis$abserror)
#}



if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MoS =  Nh_MoS_analysis$abserror)
}

#if(ncol(Nh_MoSvis_dataframe) !=  0) {
#  graph_data_abserror = cbind(graph_data_abserror, Nh_MoSvis =  Nh_MoSvis_analysis$abserror)
#}



if(ncol(Nh_GNSUM_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_GNSUM  =  Nh_GNSUM_analysis$abserror)
}




plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_notdisjoint_abserror.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)


ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Mean Absolute Error")

dev.off()



## Disjoint ##

#Dataframe creation

graph_data_abserror_disjoint = data.frame(data = simulation_data_disjoint$data)

if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$abserror)
}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$abserror)
#}


if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$abserror)
}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$abserror)
#}

if(ncol(Nh_PIMLE_dataframe_disjoint) !=  0) {
  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$abserror)
}

#if(ncol(Nh_PIMLEvis_dataframe_disjoint) !=  0) {
#  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$abserror)
#}


if(ncol(Nh_MLE_dataframe_disjoint) !=  0) {
  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$abserror)
}

#if(ncol(Nh_MLEvis_dataframe_disjoint) !=  0) {
#  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$abserror)
#}


if(ncol(Nh_MoS_dataframe_disjoint) !=  0) {
  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$abserror)
}

#if(ncol(Nh_MoSvis_dataframe_disjoint) !=  0) {
#  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$abserror)
#}


if(ncol(Nh_GNSUM_dataframe_disjoint) !=  0) {
  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$abserror)
}


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_disjoint_abserror.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)


png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_abserror_disjoint) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Mean Absolute Error")


dev.off()




## Not disjoint & disjoint ##

graph_data_abserror_total = cbind(graph_data_abserror, graph_data_abserror_disjoint[2:ncol(graph_data_abserror_disjoint)])


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_total_abserror.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_abserror_total) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Mean Absolute Error")


dev.off()



################################################################################


##### Mean of Squares Error (MSE) #####

## Not disjoint ##

#Dataframe creation

graph_data_mse = data.frame(data = simulation_data$data)

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
# # graph_data_mse = cbind(graph_data_mse, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$mse)
#}

if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_PIMLE =  Nh_PIMLE_analysis$mse)
}

#if(ncol(Nh_PIMLEvis_dataframe) !=  0) {
#  graph_data_mse = cbind(graph_data_mse, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$mse)
#}


if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MLE =  Nh_MLE_analysis$mse)
}

#if(ncol(Nh_MLEvis_dataframe) !=  0) {
#  graph_data_mse = cbind(graph_data_mse, Nh_MLEvis =  Nh_MLEvis_analysis$mse)
#}



if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MoS =  Nh_MoS_analysis$mse)
}

#if(ncol(Nh_MoSvis_dataframe) !=  0) {
#  graph_data_mse = cbind(graph_data_mse, Nh_MoSvis =  Nh_MoSvis_analysis$mse)
#}



if(ncol(Nh_GNSUM_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_GNSUM  =  Nh_GNSUM_analysis$mse)
}




plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_notdisjoint_mse.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)


ggplot(graph_data_mse) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Mean of Squares Error")

dev.off()



## Disjoint ##

#Dataframe creation

graph_data_mse_disjoint = data.frame(data = simulation_data_disjoint$data)

if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$mse)
}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$mse)
#}


if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$mse)
}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$mse)
#}

if(ncol(Nh_PIMLE_dataframe_disjoint) !=  0) {
  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$mse)
}

#if(ncol(Nh_PIMLEvis_dataframe_disjoint) !=  0) {
#  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$mse)
#}


if(ncol(Nh_MLE_dataframe_disjoint) !=  0) {
  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$mse)
}

#if(ncol(Nh_MLEvis_dataframe_disjoint) !=  0) {
#  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$mse)
#}


if(ncol(Nh_MoS_dataframe_disjoint) !=  0) {
  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$mse)
}

#if(ncol(Nh_MoSvis_dataframe_disjoint) !=  0) {
#  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$mse)
#}


if(ncol(Nh_GNSUM_dataframe_disjoint) !=  0) {
  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$mse)
}


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_disjoint_mse.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_mse_disjoint) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Mean of Squares Error")


dev.off()




## Not disjoint & disjoint ##

graph_data_mse_total = cbind(graph_data_mse, graph_data_mse_disjoint[2:ncol(graph_data_mse_disjoint)])


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_total_mse.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_mse_total) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Mean of Squares Error")


dev.off()





################################################################################

###### Bias analysis ######


## Not disjoint ##

#Dataframe creation

graph_data_bias = data.frame(data = simulation_data$data)

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
# # graph_data_bias = cbind(graph_data_bias, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$bias)
#}

if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_PIMLE =  Nh_PIMLE_analysis$bias)
}

#if(ncol(Nh_PIMLEvis_dataframe) !=  0) {
#  graph_data_bias = cbind(graph_data_bias, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$bias)
#}


if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MLE =  Nh_MLE_analysis$bias)
}

#if(ncol(Nh_MLEvis_dataframe) !=  0) {
#  graph_data_bias = cbind(graph_data_bias, Nh_MLEvis =  Nh_MLEvis_analysis$bias)
#}



if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MoS =  Nh_MoS_analysis$bias)
}

#if(ncol(Nh_MoSvis_dataframe) !=  0) {
#  graph_data_bias = cbind(graph_data_bias, Nh_MoSvis =  Nh_MoSvis_analysis$bias)
#}



if(ncol(Nh_GNSUM_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_GNSUM  =  Nh_GNSUM_analysis$bias)
}




plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_notdisjoint_bias.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)


ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Hidden population estimate")

dev.off()



## Disjoint ##

#Dataframe creation

graph_data_bias_disjoint = data.frame(data = simulation_data_disjoint$data)

if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$bias)
}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$bias)
#}


if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$bias)
}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$bias)
#}

if(ncol(Nh_PIMLE_dataframe_disjoint) !=  0) {
  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$bias)
}

#if(ncol(Nh_PIMLEvis_dataframe_disjoint) !=  0) {
#  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$bias)
#}


if(ncol(Nh_MLE_dataframe_disjoint) !=  0) {
  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$bias)
}

#if(ncol(Nh_MLEvis_dataframe_disjoint) !=  0) {
#  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$bias)
#}


if(ncol(Nh_MoS_dataframe_disjoint) !=  0) {
  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$bias)
}

#if(ncol(Nh_MoSvis_dataframe_disjoint) !=  0) {
#  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$bias)
#}


if(ncol(Nh_GNSUM_dataframe_disjoint) !=  0) {
  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$bias)
}


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_disjoint_bias.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_bias_disjoint) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Hidden population estimate")


dev.off()




## Not disjoint & disjoint ##

graph_data_bias_total = cbind(graph_data_bias, graph_data_bias_disjoint[2:ncol(graph_data_bias_disjoint)])


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_total_bias.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_bias_total) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Hidden population estimate")


dev.off()





################################################################################

#### Standard deviation analysis ####

#Dataframe creation


## Not disjoint ##

#Dataframe creation

graph_data_sd = data.frame(data = simulation_data$data)

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
# # graph_data_sd = cbind(graph_data_sd, Nh_basicvis_mean =  Nh_basicvis_mean_analysis$sd)
#}

if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_PIMLE =  Nh_PIMLE_analysis$sd)
}

#if(ncol(Nh_PIMLEvis_dataframe) !=  0) {
#  graph_data_sd = cbind(graph_data_sd, Nh_PIMLEvis =  Nh_PIMLEvis_analysis$sd)
#}


if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_MLE =  Nh_MLE_analysis$sd)
}

#if(ncol(Nh_MLEvis_dataframe) !=  0) {
#  graph_data_sd = cbind(graph_data_sd, Nh_MLEvis =  Nh_MLEvis_analysis$sd)
#}



if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_MoS =  Nh_MoS_analysis$sd)
}

#if(ncol(Nh_MoSvis_dataframe) !=  0) {
#  graph_data_sd = cbind(graph_data_sd, Nh_MoSvis =  Nh_MoSvis_analysis$sd)
#}



if(ncol(Nh_GNSUM_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_GNSUM  =  Nh_GNSUM_analysis$sd)
}


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_notdisjoint_sd.png")
sub_title = str_c("Not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)


ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Standard deviation")

dev.off()



## Disjoint ##

#Dataframe creation

graph_data_sd_disjoint = data.frame(data = simulation_data_disjoint$data)

if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$sd)
}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$sd)
#}


if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$sd)
}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$sd)
#}

if(ncol(Nh_PIMLE_dataframe_disjoint) !=  0) {
  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_PIMLE_disjoint =  Nh_PIMLE_analysis_disjoint$sd)
}

#if(ncol(Nh_PIMLEvis_dataframe_disjoint) !=  0) {
#  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_PIMLEvis_disjoint =  Nh_PIMLEvis_analysis_disjoint$sd)
#}


if(ncol(Nh_MLE_dataframe_disjoint) !=  0) {
  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MLE_disjoint =  Nh_MLE_analysis_disjoint$sd)
}

#if(ncol(Nh_MLEvis_dataframe_disjoint) !=  0) {
#  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MLEvis_disjoint =  Nh_MLEvis_analysis_disjoint$sd)
#}


if(ncol(Nh_MoS_dataframe_disjoint) !=  0) {
  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MoS_disjoint =  Nh_MoS_analysis_disjoint$sd)
}

#if(ncol(Nh_MoSvis_dataframe_disjoint) !=  0) {
#  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_MoSvis_disjoint =  Nh_MoSvis_analysis_disjoint$sd)
#}


if(ncol(Nh_GNSUM_dataframe_disjoint) !=  0) {
  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_GNSUM_disjoint  =  Nh_GNSUM_analysis_disjoint$sd)
}


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_disjoint_sd.png")
sub_title = str_c("Disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_sd_disjoint) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Standard deviation")


dev.off()




## Not disjoint & disjoint ##

graph_data_sd_total = cbind(graph_data_sd, graph_data_sd_disjoint[2:ncol(graph_data_sd_disjoint)])


plot_name = str_c("Simulation_probabilitynetwork_", seed_number, "_total_sd.png")
sub_title = str_c("Disjoint & not disjoint populations plot, seed ", seed_number)

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_sd_total) + 
  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  #geom_line(aes(x = data, y =  Nh_PIMLEvis_disjoint, col = "Nh_PIMLEvis_disjoint")) + 
  geom_line(aes(x = data, y =  Nh_PIMLE_disjoint, col = "Nh_PIMLE_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MLE_disjoint, col = "Nh_MLE_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MLEvis_disjoint, col = "Nh_MLEvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_MoS_disjoint, col = "Nh_MoS_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_MoSvis_disjoint, col = "Nh_MoSvis_disjoint")) + 
  
  geom_line(aes(x = data, y =  Nh_GNSUM_disjoint, col = "Nh_GNSUM_disjoint")) +   
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the network probability",
       subtitle = sub_title,
       x = "Network probability",
       y = "Standard deviation")


dev.off()

