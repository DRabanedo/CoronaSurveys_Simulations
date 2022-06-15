###################################
# Estimador de suma de medias (MoS)
###################################

#Cargamos todos los datos necesarios para que funcionen bien los programas
N = 10000 #Tamaño de la población
poblaciones_dis = c(0:10) #Poblaciones a distinguir. Son disjuntas y el 0 corresponde a no clasificar al individuo en ninguna
v_probabilidades = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probabilidad de cada una de las poblaciones anteriores
PropHiddenPop = 0.1 #Probabilidad de un individuo de pertenecer a la población oculta
n_encuesta = 300 #Número de muestras que extraemos en la encuesta
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
nei = 52   #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1  #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio


#Generamos la poblacion y la encuesta
Poblacion = getPoblacionTotal(N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor,vect_memory_subpob, vect_prob_memory_subpob,par_poison)
encuesta = getEncuesta(n_encuesta,Poblacion)


#Vector de poblaciones, calcula el número de personas en cada subpoblación

v_pob_totales = rep(NA, n_poblaciones)
for (k in 1:n_poblaciones) {
  v_pob_totales[k] = sum(Poblacion$Poblacion == k) # N_k
  
}

################################################################################



getNh_MoS = function(enc, v_pob, N){
  #enc es la encuesta, ?ltimos valores para Known Population
  #N es el tamaño de la población TOTAL
  #v_pob es el vector con el tamaño de cada subpoblación
  #Primero calculamos d_i_est, que es la estimación de la red personal a través
  #de las L = n_poblaciones poblaciones diferentes
  #Despues usamos la fórmula utilizada en MoS para calcular el estimador. En este 
  #caso no estamos teniendo en cuenta la visibilidad, pero si estamos corrigiendo el error de memoria
  
  # \hat{d_i} = N/L \sum_k (y_{ik}/N_k)
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$HP_total_apvis/d_i_est)
  Nh_f
}

getNh_MoSvis = function(enc, v_pob, N, vis){
  #enc es la encuesta, ?ltimos valores para Known Population
  #N es el tamaño de la población TOTAL
  #v_pob es el vector con el tamaño de cada subpoblación
  #vis es el factor de visibilidad
  #Primero calculamos d_i_est, que es la estimación de la red personal a través
  #de las L = n_poblaciones poblaciones diferentes
  #Despues usamos la fórmula utilizada en MoS para calcular el estimador. En este 
  #caso estamos corrigiendo la visibilidad y el error de memoria.
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$HP_total_apvis/d_i_est) * (1/vis)
  Nh_f
}

getNh_MoSsinvis  = function(enc, v_pob, N){
  #enc es la encuesta, ?ltimos valores para Known Population
  #N es el tamaño de la población TOTAL
  #v_pob es el vector con el tamaño de cada subpoblación
  #Primero calculamos d_i_est, que es la estimación de la red personal a través
  #de las L = n_poblaciones poblaciones diferentes
  #Despues usamos la fórmula utilizada en MoS para calcular el estimador. En este 
  #caso  estamos corrigiendo el error de memoria y suponiendo que no existe el error
  #de visibilidad
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  
  Nh_f = N * mean(enc$HP_total_conocida/d_i_est)
  Nh_f
}

################################################################################


#Valores de los estimadores
Nh_MoS = getNh_MoS(encuesta, v_pob_totales, N)
Nh_MoS
Nh_MoSvis = getNh_MoSvis(encuesta, v_pob_totales, N, visibility_factor)
Nh_MoSvis
Nh_MoSsinvis = getNh_MoSsinvis(encuesta, v_pob_totales, N)
Nh_MoSsinvis


#Valor real del tamaño de la población oculta
sum(Poblacion$Poblacion_Oculta)


