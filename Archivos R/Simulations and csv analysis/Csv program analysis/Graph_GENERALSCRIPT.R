simulation_data = read.csv("~/GitHub/CoronasurveysSimulations/Archivos R/Simulaciones y tratamiento de csv/Archivos csv/Tamaño encuesta/Poblacion fuera del bucle/Simulaciones_tamañoencuestas_ig.txt")

# WARNING: The archive is imported with the rows name #

k = 3 # Column number in which the simulations starts

# Auxiliar parameters for the calcs
n_row = nrow(simulation_data) # Number of rows
n_estimators = 13 # Number of estimators used
n_simulations = (ncol(simulation_data)-k+1)/n_estimators #Represented as b in other programs


#################### BIAS ######################

#Variable reset

Nh_real = rep(0, length(n_row))

Nh_basic_sum = rep(0, length(n_row))
Nh_basicvis_sum = rep(0, length(n_row))

Nh_basic_mean = rep(0, length(n_row))
Nh_basicvis_mean = rep(0, length(n_row))

Nh_PIMLE = rep(0, length(n_row))
Nh_PIMLEvis = rep(0, length(n_row))

Nh_MLE = rep(0, length(n_row))
Nh_MLEvis = rep(0, length(n_row))

Nh_MoS = rep(0, length(n_row))
Nh_MoSvis = rep(0, length(n_row))

Nh_GNSUM = rep(0, length(n_row))

Nh_Direct = rep(0, length(n_row))


#Bias calculations

for (i in 1:n_simulations) {
  Nh_real =  simulation_data[,k+(i-1)*n_estimators] + Nh_real
  
  Nh_basic_sum =      simulation_data[,(k+1)+(i-1)*n_estimators] +  Nh_basic_sum
  Nh_basicvis_sum =   simulation_data[,(k+2)+(i-1)*n_estimators] +  Nh_basicvis_sum
  
  Nh_basic_mean =      simulation_data[,(k+3)+(i-1)*n_estimators] +  Nh_basic_mean
  Nh_basicvis_mean =   simulation_data[,(k+4)+(i-1)*n_estimators] +  Nh_basicvis_mean
  
  Nh_PIMLE =      simulation_data[,(k+5)+(i-1)*n_estimators] + Nh_PIMLE
  Nh_PIMLEvis =   simulation_data[,(k+6)+(i-1)*n_estimators] + Nh_PIMLEvis
  
  Nh_MLE =        simulation_data[,(k+7)+(i-1)*n_estimators] + Nh_MLE
  Nh_MLEvis =     simulation_data[,(k+8)+(i-1)*n_estimators] + Nh_MLEvis
  
  Nh_MoS =        simulation_data[,(k+9)+(i-1)*n_estimators] + Nh_MoS
  Nh_MoSvis =     simulation_data[,(k+10)+(i-1)*n_estimators] + Nh_MoSvis
  
  Nh_GNSUM =      simulation_data[,(k+11)+(i-1)*n_estimators] + Nh_GNSUM
  
  Nh_Direct =    simulation_data[,(k+12)+(i-1)*n_estimators] + Nh_Direct
}


graph_data = data.frame( data = simulation_data$data, 
                         Nh_real     =  Nh_real/n_simulations,
                         
                         Nh_basic_sum    =  Nh_basic_sum/n_simulations,
                         Nh_basicvis_sum =  Nh_basicvis_sum/n_simulations,
                         
                         Nh_basic_mean    =  Nh_basic_mean/n_simulations,
                         Nh_basicvis_mean =  Nh_basicvis_mean/n_simulations,
                         
                         Nh_PIMLE    =  Nh_PIMLE/n_simulations,
                         Nh_PIMLEvis =  Nh_PIMLEvis/n_simulations,
                         
                         Nh_MLE      =  Nh_MLE/n_simulations,
                         Nh_MLEvis   =  Nh_MLEvis/n_simulations,
                         
                         Nh_MoS      =  Nh_MoS/n_simulations,
                         Nh_MoSvis   =  Nh_MoSvis/n_simulations,
                         
                         Nh_Direct  =  Nh_Direct/n_simulations,
                         
                         Nh_GNSUM    =  Nh_GNSUM/n_simulations)


#Graph 

library(ggplot2)
ggplot(graph_data) + 
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
  
  geom_line(aes(x = data, y =  Nh_real, col = "Nh_real")) +
  
  geom_line(aes(x = data, y =  Nh_Direct, col = "Nh_Direct"))+
  scale_color_discrete("Legend") + 
  labs(title = "",
       x = "",
       y = "Hidden population estimate")


#write.csv(graph_data,   # Data frame
#file = "_sd",           # CSV name
#row.names = TRUE )      # Row names: TRUE or FALSE

################################################################################

############## MEAN SQUARED ERROR ################

#Variable reset

Nh_real = rep(0, length(n_row))

Nh_basic_sum = rep(0, length(n_row))
Nh_basicvis_sum = rep(0, length(n_row))

Nh_basic_mean = rep(0, length(n_row))
Nh_basicvis_mean = rep(0, length(n_row))

Nh_PIMLE = rep(0, length(n_row))
Nh_PIMLEvis = rep(0, length(n_row))

Nh_MLE = rep(0, length(n_row))
Nh_MLEvis = rep(0, length(n_row))

Nh_MoS = rep(0, length(n_row))
Nh_MoSvis = rep(0, length(n_row))

Nh_GNSUM = rep(0, length(n_row))

Nh_Direct = rep(0, length(n_row))


#Mean squared error calculations

for (i in 1:n_simulations) {
  
  Nh_basic_sum =      (simulation_data[,(k+1)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators])^2 +  Nh_basic
  Nh_basicvis_sum =   (simulation_data[,(k+2)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators])^2 +  Nh_basicvis
  
  Nh_basic_mean =      (simulation_data[,(k+3)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators])^2 +  Nh_basic
  Nh_basicvis_mean =   (simulation_data[,(k+4)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators])^2 +  Nh_basicvis
  
  Nh_PIMLE =      (simulation_data[,(k+5)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_PIMLE
  Nh_PIMLEvis =   (simulation_data[,(k+6)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_PIMLEvis
  
  Nh_MLE =        (simulation_data[,(k+7)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_MLE
  Nh_MLEvis =     (simulation_data[,(k+8)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_MLEvis
  
  Nh_MoS =        (simulation_data[,(k+9)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_MoS
  Nh_MoSvis =     (simulation_data[,(k+10)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_MoSvis
  
  Nh_GNSUM =      (simulation_data[,(k+11)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_GNSUM
  
  Nh_Direct =    (simulation_data[,(k+12)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators])^2 + Nh_Direct
}

graph_data = data.frame( data = simulation_data$data, 
                         Nh_basic_sum    =  Nh_basic_sum/n_simulations,
                         Nh_basicvis_sum =  Nh_basicvis_sum/n_simulations,
                         
                         Nh_basic_mean    =  Nh_basic_mean/n_simulations,
                         Nh_basicvis_mean =  Nh_basicvis_mean/n_simulations,
                         
                         Nh_PIMLE    =  Nh_PIMLE/n_simulations,
                         Nh_PIMLEvis =  Nh_PIMLEvis/n_simulations,
                         
                         Nh_MLE      =  Nh_MLE/n_simulations,
                         Nh_MLEvis   =  Nh_MLEvis/n_simulations,
                         
                         Nh_MoS      =  Nh_MoS/n_simulations,
                         Nh_MoSvis   =  Nh_MoSvis/n_simulations,
                         
                         Nh_GNSUM    =  Nh_GNSUM/n_simulations,
                         
                         Nh_Direct  =  Nh_Direct/n_simulations)


library(ggplot2)
ggplot(graph_data) + 
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
  
  geom_line(aes(x = data, y =  Nh_Direct, col = "Nh_Direct"))
scale_color_discrete("Legend") + 
  labs(title = "",
       x = "",
       y = "Mean Squared Error (MSE)")



#write.csv(graph_data,                    # Data frame to be exported
#file = "RELLENAR_ECM",         # CSV name
#row.names = TRUE )             # Row names: TRUE o no FALSE

################################################################################

########### ABSOLUTE ERROR ##############

#Variable reset

Nh_real = rep(0, length(n_row))

Nh_basic_sum = rep(0, length(n_row))
Nh_basicvis_sum = rep(0, length(n_row))

Nh_basic_mean = rep(0, length(n_row))
Nh_basicvis_mean = rep(0, length(n_row))

Nh_PIMLE = rep(0, length(n_row))
Nh_PIMLEvis = rep(0, length(n_row))

Nh_MLE = rep(0, length(n_row))
Nh_MLEvis = rep(0, length(n_row))

Nh_MoS = rep(0, length(n_row))
Nh_MoSvis = rep(0, length(n_row))

Nh_GNSUM = rep(0, length(n_row))

Nh_Direct = rep(0, length(n_row))


# Absolute error calculation

for (i in 1:n_simulations) {
  
  Nh_basic_sum =      abs(simulation_data[,(k+1)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators]) +  Nh_basic_sum
  Nh_basicvis_sum =   abs(simulation_data[,(k+2)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators]) +  Nh_basicvis_sum
  
  Nh_basic_mean =      abs(simulation_data[,(k+3)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators]) +  Nh_basic_mean
  Nh_basicvis_mean =   abs(simulation_data[,(k+4)+(i-1)*n_estimators] - simulation_data[,k+(i-1)*n_estimators]) +  Nh_basicvis_mean
  
  Nh_PIMLE =      abs(simulation_data[,(k+5)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_PIMLE
  Nh_PIMLEvis =   abs(simulation_data[,(k+6)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_PIMLEvis
  
  Nh_MLE =        abs(simulation_data[,(k+7)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_MLE
  Nh_MLEvis =     abs(simulation_data[,(k+8)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_MLEvis
  
  Nh_MoS =        abs(simulation_data[,(k+9)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_MoS
  Nh_MoSvis =     abs(simulation_data[,(k+10)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_MoSvis
  
  Nh_GNSUM =      abs(simulation_data[,(k+11)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_GNSUM
  
  Nh_Direct =    abs(simulation_data[,(k+12)+(i-1)*n_estimators] -  simulation_data[,k+(i-1)*n_estimators]) + Nh_Direct
  
}


graph_data = data.frame( data = simulation_data$data, 
                         Nh_basic_sum    =  Nh_basic_sum/n_simulations,
                         Nh_basicvis_sum =  Nh_basicvis_sum/n_simulations,
                         
                         Nh_basic_mean    =  Nh_basic_mean/n_simulations,
                         Nh_basicvis_mean =  Nh_basicvis_mean/n_simulations,
                         
                         Nh_PIMLE    =  Nh_PIMLE/n_simulations,
                         Nh_PIMLEvis =  Nh_PIMLEvis/n_simulations,
                         
                         Nh_MLE      =  Nh_MLE/n_simulations,
                         Nh_MLEvis   =  Nh_MLEvis/n_simulations,
                         
                         Nh_MoS      =  Nh_MoS/n_simulations,
                         Nh_MoSvis   =  Nh_MoSvis/n_simulations,
                         
                         Nh_GNSUM    =  Nh_GNSUM/n_simulations,
                         
                         Nh_Direct  =  Nh_Direct/n_simulations)

#Graph

library(ggplot2)
ggplot(graph_data) + 
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
  
  geom_line(aes(x = data, y =  Nh_Direct, col = "Nh_Direct")) +
  scale_color_discrete("Legend") + 
  labs(title = "",
       x = "",
       y = "Mean Absolute Error (MAE")

#write.csv(graph_data,   # Data frame
#file = "_sd",           # CSV name
#row.names = TRUE )      # Row names: TRUE or FALSE

################################################################################

############# STANDARD DEVIATION ###############

#Variable reset

Nh_real = rep(0, length(n_row))

Nh_basic_sum = rep(0, length(n_row))
Nh_basicvis_sum = rep(0, length(n_row))

Nh_basic_mean = rep(0, length(n_row))
Nh_basicvis_mean = rep(0, length(n_row))

Nh_PIMLE = rep(0, length(n_row))
Nh_PIMLEvis = rep(0, length(n_row))

Nh_MLE = rep(0, length(n_row))
Nh_MLEvis = rep(0, length(n_row))

Nh_MoS = rep(0, length(n_row))
Nh_MoSvis = rep(0, length(n_row)
)
Nh_GNSUM = rep(0, length(n_row))

Nh_Direct = rep(0, length(n_row))


#Standard deviation calculation

for (i in 1:n_simulations) {
  Nh_real[,i] =       simulation_data[,k+(i-1)*n_estimators] 
  
  Nh_basic_sum[,i] =      simulation_data[,(k+1)+(i-1)*n_estimators] 
  Nh_basicvis_sum[,i] =   simulation_data[,(k+2)+(i-1)*n_estimators]
  
  Nh_basic_mean[,i] =      simulation_data[,(k+3)+(i-1)*n_estimators] 
  Nh_basicvis_mean[,i] =   simulation_data[,(k+4)+(i-1)*n_estimators]
  
  Nh_PIMLE[,i] =      simulation_data[,(k+5)+(i-1)*n_estimators]
  Nh_PIMLEvis[,i] =   simulation_data[,(k+6)+(i-1)*n_estimators] 
  Nh_MLE[,i] =        simulation_data[,(k+7)+(i-1)*n_estimators] 
  Nh_MLEvis[,i] =     simulation_data[,(k+8)+(i-1)*n_estimators]
  
  Nh_MoS[,i] =        simulation_data[,(k+9)+(i-1)*n_estimators] 
  Nh_MoSvis[,i] =     simulation_data[,(k+10)+(i-1)*n_estimators] 
  
  Nh_GNSUM[,i] =      simulation_data[,(k+11)+(i-1)*n_estimators]
  
  Nh_Direct[,i] =    simulation_data[,(k+12)+(i-1)*n_estimators]
}

# Variable reset

Nh_real_sd = rep(0, length(n_row))

Nh_basic_sd_sum = rep(0, length(n_row))
Nh_basicvis_sd_sum = rep(0, length(n_row))

Nh_basic_sd = rep(0, length(n_row))
Nh_basicvis_sd = rep(0, length(n_row))

Nh_PIMLE_sd = rep(0, length(n_row))
Nh_PIMLEvis_sd = rep(0, length(n_row))

Nh_MLE_sd = rep(0, length(n_row))
Nh_MLEvis_sd = rep(0, length(n_row))

Nh_MoS_sd = rep(0, length(n_row))
Nh_MoSvis_sd = rep(0, length(n_row))

Nh_GNSUM_sd = rep(0, length(n_row))

Nh_Direct_sd = rep(0, length(n_row))

for (i in 1:n_row) {
  Nh_real_sd[i]     = sd(Nh_real[i,])
  
  Nh_basic_sum_sd[i]    = sd(Nh_basic_sum[i,])
  Nh_basicvis_sum_sd[i] = sd(Nh_basicvis_sum[i,])
  
  Nh_basic_mean_sd[i]    = sd(Nh_basic_mean[i,])
  Nh_basicvis_mean_sd[i] = sd(Nh_basicvis_mean[i,])
  
  Nh_PIMLE_sd[i]    = sd(Nh_PIMLE[i,])
  Nh_PIMLEvis_sd[i] = sd(Nh_PIMLEvis[i,])
  
  Nh_MLE_sd[i]      = sd(Nh_MLE[i,])
  Nh_MLEvis_sd[i]   = sd(Nh_MLEvis[i,])
  
  Nh_MoS_sd[i]      = sd(Nh_MoS[i,])
  Nh_MoSvis_sd[i]   = sd(Nh_MoSvis[i,])
  
  Nh_GNSUM_sd[i]    = sd(Nh_GNSUM[i,])
  
  Nh_Direct_sd[i]  = sd(Nh_Direct[i,])
  
}

#Graph

graph_data = data.frame( datos = simulation_data$datos, 
                         Nh_basic_sum     =  Nh_basic_sum_sd,
                         Nh_basicvis_sum  =  Nh_basicvis_sum_sd,
                         
                         Nh_basic_mean    =  Nh_basic_mean_sd,
                         Nh_basicvis_mean =  Nh_basicvis_mean_sd,
                         
                         Nh_PIMLE    =  Nh_PIMLE_sd,
                         Nh_PIMLEvis =  Nh_PIMLEvis_sd,
                         
                         Nh_MLE      =  Nh_MLE_sd,
                         Nh_MLEvis   =  Nh_MLEvis_sd,
                         
                         Nh_MoS      =  Nh_MoS_sd,
                         Nh_MoSvis   =  Nh_MoSvis_sd,
                         
                         Nh_GNSUM    =  Nh_GNSUM_sd,
                         
                         Nh_Direct   =  Nh_Direct_sd)


#Graph

library(ggplot2)
ggplot(graph_data) + 
  geom_line(aes(x = datos, y =  Nh_basic_sum, col = "Nh_basic_sum")) + 
  #geom_line(aes(x = datos, y =  Nh_basicvis_sum, col = "Nh_basicvis_sum")) + 
  
  geom_line(aes(x = datos, y =  Nh_basic_mean, col = "Nh_basic_mean")) + 
  #geom_line(aes(x = datos, y =  Nh_basicvis_mean, col = "Nh_basicvis_mean")) +
  
  #geom_line(aes(x = datos, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = datos, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  
  geom_line(aes(x = datos, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = datos, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  
  geom_line(aes(x = datos, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = datos, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  
  geom_line(aes(x = datos, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  
  geom_line(aes(x = datos, y =  Nh_Direct, col = "Nh_Direct"))+
  
  scale_color_discrete("Legend") + 
  labs(title = "",
       x = "",
       y = "Standard deviation")


#write.csv(graph_data,    # Data frame
#file = "RELLENAR_sd",    # CSV name
#row.names = TRUE )       # Row names: TRUE or FALSE


