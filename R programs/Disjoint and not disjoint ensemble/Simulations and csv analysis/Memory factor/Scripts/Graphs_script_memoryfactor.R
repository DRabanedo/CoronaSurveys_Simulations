library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)

simulation_data = read.csv("C:/Users/David Rabanedo/Documents/GitHub/CoronaSurveys_Simulations/R programs/Disjoint and no disjoint ensemble/Simulations and csv analysis/Memory factor/CSV/Simulations_memoryfactor_dt_notdisjoint_2022.csv")
simulation_data_disjoint = read.csv("C:/Users/David Rabanedo/Documents/GitHub/CoronaSurveys_Simulations/R programs/Disjoint and no disjoint ensemble/Simulations and csv analysis/Memory factor/CSV/Simulations_memoryfactor_dt_disjoint_2022.csv")

seed_number = "2022"
getwd()

##################
## Not disjoint ##

Nh_real_dataframe = select(simulation_data, starts_with("Nh_real"))

Nh_basic_sum_dataframe = select(simulation_data, starts_with("Nh_basic_sum"))
#Nh_basicvis_sum_dataframe = select(simulation_data, starts_with("Nh_basicvis_sum"))

Nh_basic_mean_dataframe = select(simulation_data, starts_with("Nh_basic_mean"))
#Nh_basicvis_mean_dataframe = select(simulation_data, starts_with("Nh_basicvis_mean"))


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


##############
## Disjoint ##


#Nh_real_dataframe_disjoint = select(simulation_data_disjoint, starts_with("Nh_real"))

#Nh_basic_sum_dataframe_disjoint = select(simulation_data_disjoint, starts_with("Nh_basic_sum"))
#Nh_basicvis_sum_dataframe_disjoint = select(simulation_data_disjoint, starts_with("Nh_basicvis_sum"))

#Nh_basic_mean_dataframe_disjoint = select(simulation_data_disjoint, starts_with("Nh_basic_mean"))
#Nh_basicvis_mean_dataframe_disjoint = select(simulation_data_disjoint, starts_with("Nh_basicvis_mean"))


######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

#Nh_basic_sum_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basic_sum_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                  # mse = rowMeans(as.matrix((Nh_basic_sum_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                  # bias = rowMeans(as.matrix(Nh_basic_sum_dataframe_disjoint)),
                                  # sd = rowSds(as.matrix(Nh_basic_sum_dataframe_disjoint))) 

#Nh_basicvis_sum_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basicvis_sum_dataframe_disjoint-Nh_real_dataframe_disjoint))),
#mse = rowMeans(as.matrix((Nh_basicvis_sum_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
#bias = rowMeans(as.matrix(Nh_basicvis_sum_dataframe_disjoint)),
#sd = rowSds(as.matrix(Nh_basicvis_sum_dataframe_disjoint)))



#Nh_basic_mean_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basic_mean_dataframe_disjoint-Nh_real_dataframe_disjoint))),
                                   # mse = rowMeans(as.matrix((Nh_basic_mean_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
                                   # bias = rowMeans(as.matrix(Nh_basic_mean_dataframe_disjoint)),
                                   # sd = rowSds(as.matrix(Nh_basic_mean_dataframe_disjoint)))

#Nh_basicvis_mean_analysis_disjoint = data.frame(abserror = rowMeans(as.matrix(abs(Nh_basicvis_mean_dataframe_disjoint-Nh_real_dataframe_disjoint))),
#mse = rowMeans(as.matrix((Nh_basicvis_mean_dataframe_disjoint-Nh_real_dataframe_disjoint)^2)),
#bias = rowMeans(as.matrix(Nh_basicvis_mean_dataframe_disjoint)),
#sd = rowSds(as.matrix(Nh_basicvis_mean_dataframe_disjoint)))



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



plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_abserror.png")


png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       x = "Memory factor",
       y = "Mean Absolute Error")

dev.off()



## Disjoint ##

#Dataframe creation

#graph_data_abserror_disjoint = data.frame(data = simulation_data_disjoint$data)

#if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
#graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$abserror)
#}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$abserror)
#}



#if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
#  graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$abserror)
#}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_abserror_disjoint = cbind(graph_data_abserror_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$abserror)
#}



#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_disjoint_abserror.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_abserror_disjoint) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
 #      x = "Memory factor",
 #      y = "Mean Absolute Error")


#dev.off()




## Not disjoint & disjoint ##

#graph_data_abserror_total = cbind(graph_data_abserror, graph_data_abserror_disjoint[2:ncol(graph_data_abserror_disjoint)])


#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_total_abserror.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_abserror_total) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 
  
#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +
  
#  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
#  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
       #x = "Memory factor",
       #y = "Mean Absolute Error")


#dev.off()






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



plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_mse.png")

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_mse) + 
  
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation memory factor",
       x = "Memory factor",
       y = "Mean Squared Error (MSE)")


dev.off()


## Disjoint ##

#Dataframe creation

#graph_data_mse_disjoint = data.frame(data = simulation_data_disjoint$data)

#if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
#graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$mse)
#}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$mse)
#}



#if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
#  graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$mse)
#}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_mse_disjoint = cbind(graph_data_mse_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$mse)
#}



#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_disjoint_mse.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_mse_disjoint) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +

#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
#      x = "Memory factor",
#      y = "Mean Absolute Error")


#dev.off()




## Not disjoint & disjoint ##

#graph_data_mse_total = cbind(graph_data_mse, graph_data_mse_disjoint[2:ncol(graph_data_mse_disjoint)])


#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_total_mse.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_mse_total) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +

#  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +

#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
#x = "Memory factor",
#y = "Mean Absolute Error")


#dev.off()







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


plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_bias.png")

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
       x = "Memory factor",
       y = "Hidden population estimate")


dev.off()



## Disjoint ##

#Dataframe creation

#graph_data_bias_disjoint = data.frame(data = simulation_data_disjoint$data)

#if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
#graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$bias)
#}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$bias)
#}



#if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
#  graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$bias)
#}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_bias_disjoint = cbind(graph_data_bias_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$bias)
#}



#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_disjoint_bias.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_bias_disjoint) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +

#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
#      x = "Memory factor",
#      y = "Mean Absolute Error")


#dev.off()




## Not disjoint & disjoint ##

#graph_data_bias_total = cbind(graph_data_bias, graph_data_bias_disjoint[2:ncol(graph_data_bias_disjoint)])


#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_total_bias.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_bias_total) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +

#  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +

#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
#x = "Memory factor",
#y = "Mean Absolute Error")


#dev.off()






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






plot_name = str_c("Simulation_memoryfactor_", seed_number, "_notdisjoint_sd.png")

png(filename = plot_name,
    width = 1000, height = 600)

ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the memory factor",
       x = "Memory factor",
       y = "Standard deviation")


dev.off()



## Disjoint ##

#Dataframe creation

#graph_data_sd_disjoint = data.frame(data = simulation_data_disjoint$data)

#if(ncol(Nh_basic_sum_dataframe_disjoint) !=  0) {
#graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basic_sum_disjoint =  Nh_basic_sum_analysis_disjoint$sd)
#}

#if(ncol(Nh_basicvis_sum_dataframe_disjoint) !=  0) {
#  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basicvis_sum_disjoint =  Nh_basicvis_sum_analysis_disjoint$sd)
#}



#if(ncol(Nh_basic_mean_dataframe_disjoint) !=  0) {
#  graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basic_mean_disjoint =  Nh_basic_mean_analysis_disjoint$sd)
#}

#if(ncol(Nh_basicvis_mean_dataframe_disjoint) !=  0) {
# # graph_data_sd_disjoint = cbind(graph_data_sd_disjoint, Nh_basicvis_mean_disjoint =  Nh_basicvis_mean_analysis_disjoint$sd)
#}



#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_disjoint_sd.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_sd_disjoint) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +

#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
#      x = "Memory factor",
#      y = "Mean Absolute Error")


#dev.off()




## Not disjoint & disjoint ##

#graph_data_sd_total = cbind(graph_data_sd, graph_data_sd_disjoint[2:ncol(graph_data_sd_disjoint)])


#plot_name = str_c("Simulation_memoryfactor_", seed_number, "_total_sd.png")

#png(filename = plot_name,
#    width = 1000, height = 600)

#ggplot(graph_data_sd_total) + 
#  geom_line(aes(x = data, y =  Nh_basic_sum_disjoint, col = "Nh_basic_sum_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum_disjoint, col = "Nh_basicvis_sum_disjoint")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean_disjoint, col = "Nh_basic_mean_disjoint")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean_disjoint, col = "Nh_basicvis_mean_disjoint")) +

#  geom_line(aes(x = data, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 

#  geom_line(aes(x = data, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
#geom_line(aes(x = data, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +

#  scale_color_discrete("Legend") + 
#  labs(title = "Simulations based on the memory factor",
#x = "Memory factor",
#y = "Mean Absolute Error")


#dev.off()
