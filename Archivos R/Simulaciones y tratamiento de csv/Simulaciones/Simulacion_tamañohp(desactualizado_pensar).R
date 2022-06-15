#################################################################
# Simulación del tamaño de la hp, dejando el resto de cosas fijas
#################################################################


#Datos para la simulación

N = 10000 #Tamaño de la población
poblaciones_dis = c(0:10) #Poblaciones a distinguir. Son disjuntas y el 0 corresponde a no clasificar al individuo en ninguna
v_probabilidades = rep(1/length(poblaciones_dis),length(poblaciones_dis)) #Probabilidad de cada una de las poblaciones anteriores
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
nei = 75   #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1  #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio

#Vector de iteraciones de la probabilidad de población oculta, haciendo tramos desde el 0.1 al 1
iteraciones = seq(from = 0.05, to = 1, length.out = 20)

#Creamos un dataframe para guardar los datos
simulaciones = data.frame(datos =  iteraciones)


set.seed(semilla)
for (l in 1:25) {
  
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
  
  
  simulaciones = cbind(simulaciones,Nh_real = Nh_real)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_real_",l)
  
  simulaciones = cbind(simulaciones,Nh_basic = Nh_basic)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basic",l)
  
  simulaciones = cbind(simulaciones,Nh_basicvis = Nh_basicvis)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_basicvis_",l)
  
  simulaciones = cbind(simulaciones,Nh_PIMLE = Nh_PIMLE)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_PIMLE_",l)
  
  simulaciones = cbind(simulaciones,Nh_PIMLEvis = Nh_PIMLEvis)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_PIMLEvis_",l)
  
  simulaciones = cbind(simulaciones,Nh_MLE = Nh_MLE)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MLE_",l)
  
  simulaciones = cbind(simulaciones,Nh_MLEvis = Nh_MLEvis)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MLEvis_",l)
  
  simulaciones = cbind(simulaciones,Nh_MoS = Nh_MoS)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MoS_",l)
  
  simulaciones = cbind(simulaciones,Nh_MoSvis = Nh_MoSvis)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_MoSvis_",l)
  
  simulaciones = cbind(simulaciones,Nh_GNSUM = Nh_GNSUM)
  names(simulaciones)[dim(simulaciones)[2]] = str_c("Nh_GNSUM_",l)
  
  
  
  print(l)
  
}

simulaciones
write.csv(simulaciones,                    # Data frame a ser exportado
          file = "Simulaciones_tamañohp", # Nombre de la hoja de Excel
          row.names = TRUE )    # Incluir los nombres de las filas (TRUE) o no (FALSE)) 
