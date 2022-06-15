datos_simulacion = read.csv("~/GitHub/CoronasurveysSimulations/Archivos R/Simulaciones y tratamiento de csv/Archivos csv/Tamaño encuesta/Poblacion fuera del bucle/Simulaciones_tamañoencuestas_ig.txt")
#CUIDADO: Importamos el archivo con el nombre de las columnas

#Hacemos las medias correspondientes a cada uno de los estimadores, y lo guardamos en un vector
#que es el que vamos a representar

k = 3 #Numero de la columna en la que empiezan los datos de las simulaciones

#Para comprobar que estamos tomando los datos correctos
n_filas = nrow(datos_simulacion)
n_estimadores = 11 #Numero de estimadores diferentes usados
n_simulaciones = (ncol(datos_simulacion)-k+1)/n_estimadores


#Hacemos las medias de todos los valores correspondientes

Nh_real = rep(0, length(n_filas))
Nh_basic = rep(0, length(n_filas))
Nh_basicvis = rep(0, length(n_filas))
Nh_PIMLE = rep(0, length(n_filas))
Nh_PIMLEvis = rep(0, length(n_filas))
Nh_MLE = rep(0, length(n_filas))
Nh_MLEvis = rep(0, length(n_filas))
Nh_MoS = rep(0, length(n_filas))
Nh_MoSvis = rep(0, length(n_filas))
Nh_GNSUM = rep(0, length(n_filas))
Nh_Directo = rep(0, length(n_filas))

for (i in 1:n_simulaciones) {
  
  Nh_basic =      abs(datos_simulacion[,(k+1)+(i-1)*n_estimadores] - datos_simulacion[,k+(i-1)*n_estimadores]) +  Nh_basic
  print(Nh_basic)
  Nh_basicvis =   abs(datos_simulacion[,(k+2)+(i-1)*n_estimadores] - datos_simulacion[,k+(i-1)*n_estimadores]) +  Nh_basicvis
  
  Nh_PIMLE =      abs(datos_simulacion[,(k+3)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_PIMLE
  Nh_PIMLEvis =   abs(datos_simulacion[,(k+4)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_PIMLEvis
  
  Nh_MLE =        abs(datos_simulacion[,(k+5)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_MLE
  Nh_MLEvis =     abs(datos_simulacion[,(k+6)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_MLEvis
  
  Nh_MoS =        abs(datos_simulacion[,(k+7)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_MoS
  Nh_MoSvis =     abs(datos_simulacion[,(k+8)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_MoSvis
  
  Nh_GNSUM =      abs(datos_simulacion[,(k+9)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_GNSUM
  
  Nh_Directo =    abs(datos_simulacion[,(k+10)+(i-1)*n_estimadores] -  datos_simulacion[,k+(i-1)*n_estimadores]) + Nh_Directo

}


datos_grafica = data.frame( datos = datos_simulacion$datos, 
                            Nh_basic    =  Nh_basic/n_simulaciones,
                            Nh_basicvis =  Nh_basicvis/n_simulaciones,
                            Nh_PIMLE    =  Nh_PIMLE/n_simulaciones,
                            Nh_PIMLEvis =  Nh_PIMLEvis/n_simulaciones,
                            Nh_MLE      =  Nh_MLE/n_simulaciones,
                            Nh_MLEvis   =  Nh_MLEvis/n_simulaciones,
                            Nh_MoS      =  Nh_MoS/n_simulaciones,
                            Nh_MoSvis   =  Nh_MoSvis/n_simulaciones,
                            Nh_GNSUM    =  Nh_GNSUM/n_simulaciones,
                            Nh_Directo  =  Nh_Directo/n_simulaciones)

library(ggplot2)
ggplot(datos_grafica) + 
  geom_line(aes(x = datos, y =  Nh_basic, col = "Nh_basic")) + 
  #geom_line(aes(x = datos, y =  Nh_basicvis, col = "Nh_basicvis")) + 
  #geom_line(aes(x = datos, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = datos, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  geom_line(aes(x = datos, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = datos, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  geom_line(aes(x = datos, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = datos, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  geom_line(aes(x = datos, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  geom_line(aes(x = datos, y =  Nh_Directo, col = "Nh_Directo")) +
  scale_color_discrete("Leyenda") + 
  labs(title = "Variabilidad en función del tamaño de la encuesta",
       x = "Tamaño de la encuesta",
       y = "Error absoluto medio")

#write.csv(datos_grafica,                    # Data frame a ser exportado
#file = "RELLENAR_var", # Nombre de la hoja de Excel
#row.names = TRUE )    # Incluir los nombres de las filas (TRUE) o no (FALSE))
