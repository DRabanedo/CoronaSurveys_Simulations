library(dplyr)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(grid)

################################################################################

# Reading the data

simulation_data = read.csv("C:/Users/David Rabanedo/Desktop/Cs_subpop size/Simulations_subpopulationsize")

table_PIMLE = read.csv("C:/Users/David Rabanedo/Desktop/Cs_subpop size/Simulations_subpopulationsize_summary_PIMLE")
table_PIMLE_Nh = read.csv("C:/Users/David Rabanedo/Desktop/Cs_subpop size/Simulations_subpopulationsize_summary_PIMLE_Nh")
table_MoS = read.csv("C:/Users/David Rabanedo/Desktop/Cs_subpop size/Simulations_subpopulationsize_summary_MoS")
table_MoS_Nh = read.csv("C:/Users/David Rabanedo/Desktop/Cs_subpop size/Simulations_subpopulationsize_summary_MoS_Nh")

################################################################################

#####################
## Table creations ##
#####################

# Table PIMLE's di creation
table_PIMLE = tableGrob(round(table_PIMLE[c(-1,-7)], digits = 3))
grid.newpage()
h <- grobHeight(table_PIMLE)
w <- grobWidth(table_PIMLE)
title <- textGrob("PIMLE/MLE di analysis", y=unit(0.5,"npc") + 1.2*h, 
                  vjust=0, gp=gpar(fontsize=20))
gt <- gTree(children=gList(table_PIMLE, title))
grid.draw(gt)


# Table PIMLE's Nh analysis
table_PIMLE_Nh = tableGrob(round(table_PIMLE_Nh[c(-1,-7)], digits = 3))
grid.newpage()
h <- grobHeight(table_PIMLE_Nh)
w <- grobWidth(table_PIMLE_Nh)
title <- textGrob("PIMLE Nh analysis", y=unit(0.5,"npc") + 1.2*h, 
                  vjust=0, gp=gpar(fontsize=20))
gt <- gTree(children=gList(table_PIMLE_Nh, title))
grid.draw(gt)



# Table MoS' di creation
table_MoS = tableGrob(round(table_MoS[c(-1,-7)], digits = 3))
grid.newpage()
h <- grobHeight(table_MoS)
w <- grobWidth(table_MoS)
title <- textGrob("MoS di analysis", y=unit(0.5,"npc") + 1.2*h, 
                  vjust=0, gp=gpar(fontsize=20))
gt <- gTree(children=gList(table_MoS, title))
grid.draw(gt)


# Table MoS Nh analysis

table_MoS_Nh = tableGrob(round(table_MoS_Nh[c(-1,-7)], digits = 3))
grid.newpage()
h <- grobHeight(table_MoS_Nh)
w <- grobWidth(table_MoS_Nh)
title <- textGrob("MoS Nh analysis", y=unit(0.5,"npc") + 1.2*h, 
                  vjust=0, gp=gpar(fontsize=20))
gt <- gTree(children=gList(table_MoS_Nh, title))
grid.draw(gt)

################################################################################

# Data selection #

Nh_real_dataframe = dplyr::select(simulation_data, starts_with("Nh_real"))

Nh_PIMLE_dataframe    = dplyr::select(simulation_data, starts_with("Nh_PIMLE_"))

Nh_MLE_dataframe     = dplyr::select(simulation_data, starts_with("Nh_MLE_"))

Nh_MoS_dataframe     = dplyr::select(simulation_data, starts_with("Nh_MoS_"))

################################################################################

######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.



Nh_PIMLE_analysis    = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_PIMLE_dataframe-Nh_real_dataframe))),
                                  mse = rowMeans(as.matrix((Nh_PIMLE_dataframe-Nh_real_dataframe)^2)),
                                  bias = rowMeans(as.matrix(Nh_PIMLE_dataframe)),
                                  sd = rowSds(as.matrix(Nh_PIMLE_dataframe)))



Nh_MLE_analysis     = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_MLE_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MLE_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MLE_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MLE_dataframe)))



Nh_MoS_analysis     = data.frame(abs_error = rowMeans(as.matrix(abs(Nh_MoS_dataframe-Nh_real_dataframe))),
                                 mse = rowMeans(as.matrix((Nh_MoS_dataframe-Nh_real_dataframe)^2)),
                                 bias = rowMeans(as.matrix(Nh_MoS_dataframe)),
                                 sd = rowSds(as.matrix(Nh_MoS_dataframe)))



################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### Absolute error #####

#Dataframe creation

graph_data_abserror = data.frame( data = simulation_data$data)



if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_PIMLE =  Nh_PIMLE_analysis$abs_error)
}



if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MLE =  Nh_MLE_analysis$abs_error)
}


if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_abserror = cbind(graph_data_abserror, Nh_MoS =  Nh_MoS_analysis$abs_error)
}



ggplot(graph_data_abserror) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       x = "Type of subpopulation",
       y = "Mean Absolute Error")


################################################################################


##### Mean of Squares Error (MSE) #####

#Dataframe creation

graph_data_mse = data.frame( data = simulation_data$data)


if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_PIMLE =  Nh_PIMLE_analysis$mse)
}


if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MLE =  Nh_MLE_analysis$mse)
}


if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_mse = cbind(graph_data_mse, Nh_MoS =  Nh_MoS_analysis$mse)
}



ggplot(graph_data_mse) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       x = "Type of subpopulation",
       y = "Mean Squared Error (MSE)")


################################################################################

###### Bias analysis ######


graph_data_bias = data.frame( data = simulation_data$data)

graph_data_bias = cbind(graph_data_bias, Nh_real =  simulation_data$Nh_real_1)


if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_PIMLE =  Nh_PIMLE_analysis$bias)
}


if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MLE =  Nh_MLE_analysis$bias)
}


if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_bias = cbind(graph_data_bias, Nh_MoS =  Nh_MoS_analysis$bias)
}


ggplot(graph_data_bias) +
  geom_point(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       x = "Type of subpopulation",
       y = "Hidden population estimate")



################################################################################

#### Standard deviation analysis ####

#Dataframe creation

graph_data_sd = data.frame( data = simulation_data$data)


if(ncol(Nh_PIMLE_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_PIMLE =  Nh_PIMLE_analysis$sd)
}


if(ncol(Nh_MLE_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_MLE =  Nh_MLE_analysis$sd)
}



if(ncol(Nh_MoS_dataframe) !=  0) {
  graph_data_sd = cbind(graph_data_sd, Nh_MoS =  Nh_MoS_analysis$sd)
}



ggplot(graph_data_sd) + 
  geom_point(aes(x = data, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MLE, col = "Nh_MLE")) + 
  
  geom_point(aes(x = data, y =  Nh_MoS, col = "Nh_MoS")) + 
  
  scale_color_discrete("Legend") + 
  labs(title = "Simulations based on the subpopulation size",
       x = "Type of subpopulation",
       y = "Standard deviation")
