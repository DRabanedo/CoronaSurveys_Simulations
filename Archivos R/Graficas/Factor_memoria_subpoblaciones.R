###########################################################################################################
# Gráfica en función del valor del factor de memoria de las subpoblaciones, dejando el resto de cosas fijas
###########################################################################################################

N = 10000                  #Tamaño de la población general
poblaciones_dis = c(0:10)  #Vector de poblaciones disjuntas
#Vector de probabilidades de pertenencia a cada una de las poblaciones
v_probabilidades = rep(1/length(poblaciones_dis), length(poblaciones_dis)) 
PropHiddenPop = 0.1        #Probabilidad de pertenecer a la población oculta
n_encuesta = 500           #Tamaño de la encuesta a la población general
n_encuesta_hp = 50         #Tamaño de la encuesta a la población oculta
#Número de poblaciones utilizadas para la simulación
n_poblaciones =  length(poblaciones_dis)-1 

################CUIDADO#########################################################
# ES IMPORTANTE TENER ESTE FACTOR A 0, YA QUE CON EL SE CONSTRUYE LA POBLACIÓN QUE
# VAMOS A TOMAR COMO REFERENCIA
sub_memory_factor = 0      #Factor de memoria de las subpoblaciones
################################################################################

memory_factor = 0      #Factor memoria 
visibility_factor = 1  #Factor visibilidad

semilla = 207          #Está guay cambiarla pq aunque se mantenga la estructura, los resultados cambian
set.seed(semilla)


#Grafo
dim = 1    #Dimensión del grafo
nei = 75   #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1  #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio


#Vector de iteraciones del factor memoria
iteraciones = seq(from = 0, to = 0.25, length.out = 21)


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

#Generamos la población y la matriz
Grafo_poblacion_matriz = getDatos(N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)
Poblacion = Grafo_poblacion_matriz[[2]]
Mhp_vis = Grafo_poblacion_matriz[[3]]

#Encuesta a la población y población oculta
encuesta = getEncuesta(n_encuesta,Poblacion)
encuesta_hp = getEncuesta(n_encuesta_hp, Poblacion[Poblacion$Poblacion_Oculta==1,])

#Calculamos el tamaño de las subpoblaciones conocidas
v_pob_totales = rep(NA, n_poblaciones)
for (j in 1:n_poblaciones) {
  v_pob_totales[j] = sum(Poblacion$Poblacion == j)
}

#Cambiamos solo las columnas necesarias para que la simulación funcione más rápido
vect_reach = encuesta$Reach
vect_reach_re =  rep(NA, nrow(encuesta))


for (i in 1:length(iteraciones)) {
  
  sub_memory_factor = iteraciones[i]   #Metemos el dato correspondiente a la memoria
  
  for(j in 0:n_poblaciones){
    v_1 = rep(NA,nrow(encuesta))
    vis_pob = encuesta[,ncol(encuesta)-(n_poblaciones-j)]
    for(k in 1:nrow(encuesta)) {
      v_1[k] = max(round(rnorm(1, mean = vis_pob[k], sd = sub_memory_factor*vis_pob[k])),0)
    }
    encuesta[,ncol(encuesta)-(n_poblaciones-j)] = v_1
  }
  
  
  Nh_real[i] = sum(Poblacion$Poblacion_Oculta) 
  
  Nh_basic[i] = getNh_basic(encuesta,N) 
  Nh_basicvis[i] = getNh_basicvis(encuesta,N,visibility_factor) 
  
  Nh_PIMLE[i] = getNh_PIMLE(encuesta, v_pob_totales, N)
  Nh_PIMLEvis[i] = getNh_PIMLEvis(encuesta, v_pob_totales, N, visibility_factor)
  
  Nh_MLE[i] = getNh_MLE(encuesta, v_pob_totales)
  Nh_MLEvis[i] = getNh_MLEvis(encuesta, v_pob_totales, visibility_factor)
  
  Nh_MoS[i] = getNh_MoS(encuesta, v_pob_totales, N)
  Nh_MoSvis[i] = getNh_MoSvis(encuesta, v_pob_totales, N, visibility_factor)
  
  Nh_GNSUM[i] =  getNh_GNSUM(Poblacion, encuesta, encuesta_hp, Mhp_vis, v_pob_totales, N)
}



################################################################################

x_1 = iteraciones
ggplot() + 
  #geom_line(aes(x = x_1, y =  Nh_basic, col = "Nh_basic")) + 
  #geom_line(aes(x = x_1, y =  Nh_basicvis, col = "Nh_basicvis")) + 
  #geom_line(aes(x = x_1, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLE, col = "Nh_MLE")) + 
  #geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_MoS, col = "Nh_MoS")) + 
  #geom_line(aes(x = x_1, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  #geom_line(aes(x = x_1, y =  Nh_GNSUM, col = "Nh_GNSUM")) +
  
  geom_line(aes(x = x_1, y =  Nh_real, col = "Nh_real")) +
  scale_color_discrete("Leyenda") + 
  labs(title = "Variabilidad en función del factor memoria de las subpoblaciones",
       x = "Factor de memoria",
       y = "Estimación de la población oculta")

