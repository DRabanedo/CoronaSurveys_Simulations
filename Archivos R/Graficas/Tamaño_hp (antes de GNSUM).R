#################################################################################################
# Gráfica en función del valor del tamaño de la población oculta, dejando el resto de cosas fijas
#################################################################################################

N = 10000 #Tamaño de la población
poblaciones_dis = c(0:10) #Poblaciones a distinguir. Son disjuntas y el 0 corresponde a no clasificar al individuo en ninguna
v_probabilidades = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probabilidad de cada una de las poblaciones anteriores
PropHiddenPop = 0.1 #Probabilidad de un individuo de pertenecer a la población oculta
n_encuesta = 300 #Número de muestras que extraemos en la encuesta
n_encuesta_hp = 50
n_poblaciones = length(poblaciones_dis)-1 #Número de poblaciones en las que clasificamos
par_poison = 20  #Parámetro de la Poisson con la que vamos a sobreestimar la variable Reach
vect_memory_subpob = c(0.8,0.9,1,1.1,1.2) #Vector del factor memoria que vamos a aplicar a las subpoblaciones
vect_prob_memory_subpob = c(0.1,0.2,0.4,0.2,0.1) #Probabilidad de aplicar cada una de las componentes del vector
memory_factor = 0.7      #Factor de memoria 
visibility_factor = 0.8  #Factor visibilidad
semilla = 207            #Está guay cambiarla pq aunque se mantenga la estructura, los resultados cambian
set.seed(semilla)

#Grafo
dim = 1    #Dimensión del grafo
nei = 75   #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1  #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio


#Vector de iteraciones de la probabilidad de población oculta, haciendo tramos desde el 0.1 al 1
iteraciones = seq(from = 0.1, to = 1, length.out = 19)

#Reiniciamos las variables
Nh_real =  rep(NA,length(iteraciones)) 

Nh_basic = rep(NA,length(iteraciones)) 
Nh_basicvis = rep(NA,length(iteraciones)) 

Nh_PIMLE = rep(NA,length(iteraciones)) 
Nh_PIMLEvis = rep(NA,length(iteraciones)) 

Nh_MLE = rep(NA,length(iteraciones)) 
Nh_MLEvis = rep(NA,length(iteraciones)) 

Nh_MoS = rep(NA,length(iteraciones)) 
Nh_MoSvis = rep(NA,length(iteraciones)) 

Nh_GNSUM = rep(NA,length(iteraciones)) 

#Generamos la población general que vamos a utilizar y la distribución de dicha población
Grafo_y_poblacion =  getDatos(N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor,vect_memory_subpob, vect_prob_memory_subpob,par_poison)
Poblacion =  Grafo_y_poblacion[[2]]
net_sw = Grafo_y_poblacion[[1]]

for (j in 1:nrow(Poblacion)) {
  vect_hp[j] = sum(Poblacion[net_sw[[j]][[1]],]$Poblacion_Oculta)
  vect_hp_vis[j] = rbinom(1,vect_hp[j],prob = visibility_factor)
}


#Contamos el tamaño de las subpoblaciones conocidas, que no dependen de la cantidad de gente en la población oculta
v_pob_totales = rep(NA, n_poblaciones)
for (j in 1:n_poblaciones) {
  v_pob_totales[j] = sum(Poblacion$Poblacion == j)
}


#Realizamos el bucle para hacer la simulación de la variación de este factor
for (i in 1:length(iteraciones)) {
  
  PropHiddenPop = iteraciones[i]   #Metemos el dato correspondiente a la población oculta
  #Generamos los nuevos vectores de población oculta y los sustituimos en el dataframe
  Poblacion$Poblacion_Oculta = sample(c(0,1), nrow(Poblacion), replace = TRUE, p = c(1-PropHiddenPop,PropHiddenPop))
  vect_hp = rep(NA,N)       #Vector sin visibility
  vect_hp_vis = rep(NA,N)   #Vector con visibility aplicado
  for (j in 1:nrow(Poblacion)) {
    vect_hp[j] = sum(Poblacion[net_sw[[j]][[1]],]$Poblacion_Oculta)
    vect_hp_vis[j] = rbinom(1,vect_hp[j],prob = visibility_factor)
  }
  Poblacion$HP_total_conocida = vect_hp
  Poblacion$HP_total_apvis = vect_hp_vis
  
  
  #Teniendo la población como queremos, hacemos una encuesta a esa población y a la oculta correspondiente
  encuesta = getEncuesta(n_encuesta,Poblacion)
  encuesta_hp = getEncuesta(n_encuesta_hp, Poblacion[Poblacion$Poblacion_Oculta==1,])
  
  #Añadimos las estimaciones de cada uno de los diferentes estimadores
  Nh_real[i] = sum(Poblacion$Poblacion_Oculta) 
  
  Nh_basic[i] = getNh_basic(encuesta,N) 
  Nh_basicvis[i] = getNh_basicvis(encuesta,N,visibility_factor) 
  
  Nh_PIMLE[i] = getNh_PIMLE(encuesta, v_pob_totales, N)
  Nh_PIMLEvis[i] = getNh_PIMLEvis(encuesta, v_pob_totales, N, visibility_factor)
  
  Nh_MLE[i] = getNh_MLE(encuesta, v_pob_totales)
  Nh_MLEvis[i] = getNh_MLEvis(encuesta, v_pob_totales, visibility_factor)
  
  Nh_MoS[i] = getNh_MoS(encuesta, v_pob_totales, N)
  Nh_MoSvis[i] = getNh_MoSvis(encuesta, v_pob_totales, N, visibility_factor)
  
  Nh_GNSUM[i] =  getNh_GNSUM(encuesta, encuesta_hp, v_pob_totales, N, visibility_factor)
}
################################################################################

#Hacemos las gráficas correspondientes a esta simulación
x_1 = iteraciones
ggplot() + 
  geom_line(aes(x = x_1, y =  Nh_basic, col = "Nh_basic")) + 
  geom_line(aes(x = x_1, y =  Nh_basicvis, col = "Nh_basicvis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLE, col = "Nh_MLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_MoS, col = "Nh_MoS")) + 
  geom_line(aes(x = x_1, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  geom_line(aes(x = x_1, y =  Nh_GNSUM, col = "Nh_GNSUM")) +
  
  geom_line(aes(x = x_1, y =  Nh_real, col = "Nh_real")) +
  scale_color_discrete("Leyenda") + 
  labs(title = "Variabilidad en función del tamaño de la población oculta",
       x = "Porcentaje de población oculta",
       y = "Estimación de la población oculta")