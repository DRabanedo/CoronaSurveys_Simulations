######################################################################################
# Gráfica en función del tamaño de las subpoblaciones, dejando el resto de cosas fijas
######################################################################################

N = 10000 #Tamaño de la población
poblaciones_dis = c(0:10) #Poblaciones a distinguir. Son disjuntas y el 0 corresponde a no clasificar al individuo en ninguna
v_probabilidades = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probabilidad de cada una de las poblaciones anteriores
PropHiddenPop = 0.1 #Probabilidad de un individuo de pertenecer a la población oculta
n_encuesta = 300 #Número de muestras que extraemos en la encuesta
n_encuesta_hp = 50
n_poblaciones = length(poblaciones_dis)-1 #Número de poblaciones en las que clasificamos
par_poison = 0  #Parámetro de la Poisson con la que vamos a sobreestimar la variable Reach
vect_memory_subpob = 1 #Vector del factor memoria que vamos a aplicar a las subpoblaciones
vect_prob_memory_subpob = 1 #Probabilidad de aplicar cada una de las componentes del vector
memory_factor = 1      #Factor de memoria 
visibility_factor = 1  #Factor visibilidad
semilla = 207            #Está guay cambiarla pq aunque se mantenga la estructura, los resultados cambian
set.seed(semilla)

#Grafo
dim = 1    #Dimensión del grafo
nei = 52   #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1  #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio



#Vector de iteraciones
iteraciones = list(rep(1/length(poblaciones_dis), length(poblaciones_dis)), c(0.5,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05),
                   c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025), c(0.2, 0.5, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.05,0.05))


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


Grafo_y_poblacion = getDatos(N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor,vect_memory_subpob, vect_prob_memory_subpob,par_poison)
Poblacion = Grafo_y_poblacion[[2]]
net_sw = Grafo_y_poblacion[[1]]
#Factor de memoria para las subpoblaciones
memory_factor_subpob = sample(vect_memory_subpob, N*(n_poblaciones+1), replace = TRUE, vect_prob_memory_subpob)
#Variable para poder seleccionar las componentes adecuadas
k = length(poblaciones_dis)

################################################################################
#Bucle para calcular la simulación
for (i in 1:length(iteraciones)) {
  
  v_probabilidades = iteraciones[[i]]
  m_pob = length(iteraciones[[i]])-1
  n_columnas = ncol(Poblacion)
  poblaciones_dis = c(0:m_pob)
  n_poblaciones = length(poblaciones_dis)-1 
  
  Poblacion$Poblacion = sample(poblaciones_dis, N, replace = TRUE, p = v_probabilidades)

  
  for(j in 0:n_poblaciones){
    v_1 = rep(NA,N)
    for(v in 1:N) {
      # Visibilidad de la poblaci?n j por v aplicando un factor de visibilidad para las subpoblaciones
      v_1[v] = sum(Poblacion[net_sw[[v]][[1]],]$Poblacion == j)
    }
    
    Poblacion[,n_columnas-2*k+1+j] = v_1
    names(Poblacion)[n_columnas-2*k+1+j] = str_c("KP_total_",j)
  }
  
  for(j in 0:n_poblaciones){
    v_1 = rep(NA,N)
    for(v in 1:N) {
      # Visibilidad de la poblaci?n j por v aplicando un factor de visibilidad para las subpoblaciones
      v_1[v] = round(sum(Poblacion[net_sw[[v]][[1]],]$Poblacion == j)*memory_factor_subpob[(n_poblaciones+1)*(v-1)+j+1])
    }
    
    Poblacion[,n_columnas-(-length(poblaciones_dis) + 2*k) +1 +j] = v_1
    names(Poblacion)[n_columnas-(-length(poblaciones_dis) + 2*k) +1 +j] = str_c("KP_total_apvis_",j)
  }
  
  k = length(poblaciones_dis)
  
  v_pob_totales = rep(NA, n_poblaciones)
  for (j in 1:n_poblaciones) {
    v_pob_totales[j] = sum(Poblacion$Poblacion == j)
  }
  #Hacemos las encuestas
  encuesta = getEncuesta(n_encuesta,Poblacion)
  encuesta_hp = getEncuesta(n_encuesta_hp, Poblacion[Poblacion$Poblacion_Oculta==1,])
  
  #Estimamos la poblacion oculta
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

#################################################################################

x_0 = c("1-Todas iguales", "2-Iguales pero pequeñas", "3-Desiguales", "4-Muy desiguales")
x_1 = 1:4
ggplot() +
  geom_line(aes(x = x_0, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_basicvis, col = "Nh_basicvis")) + 
  geom_line(aes(x = x_1,  y =  Nh_basic, col = "Nh_basic")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLEvis, col = "Nh_PIMLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_PIMLE, col = "Nh_PIMLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLE, col = "Nh_MLE")) + 
  geom_line(aes(x = x_1, y =  Nh_MLEvis, col = "Nh_MLEvis")) + 
  geom_line(aes(x = x_1, y =  Nh_MoS, col = "Nh_MoS")) + 
  geom_line(aes(x = x_1, y =  Nh_MoSvis, col = "Nh_MoSvis")) + 
  geom_line(aes(x = x_1, y =  Nh_GNSUM, col = "Nh_GNSUM")) + 
  geom_line(aes(x = x_1, y =  Nh_real, col = "Nh_real")) +
  scale_color_discrete("Leyenda") + 
  labs(title = "Variabilidad en función de las subpoblaciones de tamaño conocido",
       x = "Número de subpoblaciones",
       y = "Estimación de la población oculta")