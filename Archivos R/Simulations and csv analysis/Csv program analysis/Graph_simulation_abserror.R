simulation_data = read.csv("~/GitHub/CoronasurveysSimulations/Archivos R/Simulaciones y tratamiento de csv/Archivos csv/Tamaño encuesta/Poblacion fuera del bucle/Simulaciones_tamañoencuestas_ig.txt")

# WARNING: The archive is imported with the rows name #

# Calculations and representation of the absolute error Error

k = 3 # Column number in which the simulations starts

# Auxiliar parameters for the calcs
n_row = nrow(simulation_data) # Number of rows
n_estimators = 11 # Number of estimators used
n_simulations = (ncol(simulation_data)-k+1)/n_estimators #Represented as b in other programs


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
