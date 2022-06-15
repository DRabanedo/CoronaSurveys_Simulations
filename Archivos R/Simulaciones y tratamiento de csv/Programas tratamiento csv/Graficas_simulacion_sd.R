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

Nh_real = matrix(0,n_filas,n_simulaciones)
Nh_basic = matrix(0,n_filas,n_simulaciones)
Nh_basicvis = matrix(0,n_filas,n_simulaciones)
Nh_PIMLE = matrix(0,n_filas,n_simulaciones)
Nh_PIMLEvis = matrix(0,n_filas,n_simulaciones)
Nh_MLE = matrix(0,n_filas,n_simulaciones)
Nh_MLEvis = matrix(0,n_filas,n_simulaciones)
Nh_MoS = matrix(0,n_filas,n_simulaciones)
Nh_MoSvis = matrix(0,n_filas,n_simulaciones)
Nh_GNSUM = matrix(0,n_filas,n_simulaciones)
Nh_Directo = matrix(0,n_filas,n_simulaciones)


for (i in 1:n_simulaciones) {
  Nh_real[,i] =       datos_simulacion[,k+(i-1)*n_estimadores] 
  
  Nh_basic[,i] =      datos_simulacion[,(k+1)+(i-1)*n_estimadores] 
  Nh_basicvis[,i] =   datos_simulacion[,(k+2)+(i-1)*n_estimadores]
  
  Nh_PIMLE[,i] =      datos_simulacion[,(k+3)+(i-1)*n_estimadores]
  Nh_PIMLEvis[,i] =   datos_simulacion[,(k+4)+(i-1)*n_estimadores] 
  Nh_MLE[,i] =        datos_simulacion[,(k+5)+(i-1)*n_estimadores] 
  Nh_MLEvis[,i] =     datos_simulacion[,(k+6)+(i-1)*n_estimadores]
  
  Nh_MoS[,i] =        datos_simulacion[,(k+7)+(i-1)*n_estimadores] 
  Nh_MoSvis[,i] =     datos_simulacion[,(k+8)+(i-1)*n_estimadores] 
  
  Nh_GNSUM[,i] =      datos_simulacion[,(k+9)+(i-1)*n_estimadores]
  
  Nh_Directo[,i] =    datos_simulacion[,(k+10)+(i-1)*n_estimadores]
}

Nh_basic

#Hacemos las medias de todos los valores correspondientes

Nh_real_sd = rep(0, length(n_filas))
Nh_basic_sd = rep(0, length(n_filas))
Nh_basicvis_sd = rep(0, length(n_filas))
Nh_PIMLE_sd = rep(0, length(n_filas))
Nh_PIMLEvis_sd = rep(0, length(n_filas))
Nh_MLE_sd = rep(0, length(n_filas))
Nh_MLEvis_sd = rep(0, length(n_filas))
Nh_MoS_sd = rep(0, length(n_filas))
Nh_MoSvis_sd = rep(0, length(n_filas))
Nh_GNSUM_sd = rep(0, length(n_filas))
Nh_Directo_sd = rep(0, length(n_filas))

for (i in 1:n_filas) {
  Nh_real_sd[i]     = sd(Nh_real[i,])
  Nh_basic_sd[i]    = sd(Nh_basic[i,])
  Nh_basicvis_sd[i] = sd(Nh_basicvis[i,])
  Nh_PIMLE_sd[i]    = sd(Nh_PIMLE[i,])
  Nh_PIMLEvis_sd[i] = sd(Nh_PIMLEvis[i,])
  Nh_MLE_sd[i]      = sd(Nh_MLE[i,])
  Nh_MLEvis_sd[i]   = sd(Nh_MLEvis[i,])
  Nh_MoS_sd[i]      = sd(Nh_MoS[i,])
  Nh_MoSvis_sd[i]   = sd(Nh_MoSvis[i,])
  Nh_GNSUM_sd[i]    = sd(Nh_GNSUM[i,])
  Nh_Directo_sd[i]  = sd(Nh_Directo[i,])
  
}
sd(Nh_basic[1,])
mean(Nh_basic[1,])
datos_grafica = data.frame( datos = datos_simulacion$datos, 
                            Nh_basic    =  Nh_basic_sd,
                            Nh_basicvis =  Nh_basicvis_sd,
                            Nh_PIMLE    =  Nh_PIMLE_sd,
                            Nh_PIMLEvis =  Nh_PIMLEvis_sd,
                            Nh_MLE      =  Nh_MLE_sd,
                            Nh_MLEvis   =  Nh_MLEvis_sd,
                            Nh_MoS      =  Nh_MoS_sd,
                            Nh_MoSvis   =  Nh_MoSvis_sd,
                            Nh_GNSUM    =  Nh_GNSUM_sd,
                            Nh_Directo  =  Nh_Directo_sd)


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
  geom_line(aes(x = datos, y =  Nh_Directo, col = "Nh_Directo"))+
  scale_color_discrete("Leyenda") + 
  labs(title = "Variabilidad en función del tamaño de la encuesta",
       x = "Tamaño de la encuesta",
       y = "Desviación típica")


#write.csv(datos_grafica,                    # Data frame a ser exportado
#file = "RELLENAR_sd", # Nombre de la hoja de Excel
#row.names = TRUE )    # Incluir los nombres de las filas (TRUE) o no (FALSE))