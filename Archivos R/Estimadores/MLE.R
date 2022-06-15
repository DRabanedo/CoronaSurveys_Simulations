######################################
# Metodo de máxima verosimilitud (MLE)
######################################

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


getNh_MLE = function(enc,v_pob) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #En esta fórmula tenemos en cuenta el error de memoria pero no el de visibilidad
  #enc es la encuesta, ?ltimos valores para las preguntas de la poblaci?n conocida
  #v_pob es un vector con los valores de las poblaciones conocidas

  suma_KP = sum(enc[tail(names(enc),length(v_pob))]) # suma de la Known Population
  # (\sum y_{iu})/(\frac{\sum N_k}{\sum \sum y_{ik}} )
  Nh_f = sum(enc$HP_total_apvis)*(sum(v_pob)/suma_KP)
  Nh_f
}

getNh_MLEvis = function(enc,v_pob,vis) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #En esta fórmula tenemos en cuenta el error de memoria y el de visibilidad, y lo corregimos.
  #enc es la encuesta, ?ltimos valores para las preguntas de la poblaci?n conocida
  #v_pob es un vector con los valores de las poblaciones conocidas
  #vis es el factor de visibilidad
  
  suma_KP = sum(enc[tail(names(enc),length(v_pob))])
  Nh_f = (sum(enc$HP_total_apvis))*(sum(v_pob)/suma_KP)*(1/vis)
  Nh_f
}

getNh_MLEsinvis =  function(enc,v_pob) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #En este caso no hay error de visibilidad, y tenemos en cuenta el error de memoria si lo hubiera
  #enc es la encuesta, ?ltimos valores para las preguntas de la poblaci?n conocida
  #v_pob es un vector con los valores de las poblaciones conocidas
  
  suma_KP = sum(enc[tail(names(enc),length(v_pob))])
  Nh_f = (sum(enc$HP_total_conocida))*(sum(v_pob)/suma_KP) 
  Nh_f
}

################################################################################

#Como vemos, se obtiene una muy buena aproximación de la variable Reach media de
#la población, luego el error de memoria se elimina gracias a la incorporación de
#las poblaciones auxiliares conocidas
sum(encuesta[8:dim(encuesta)[2]])/sum(v_pob_totales)*(N/n_encuesta)  #Estimación de la variable Reach de media
sum(Poblacion$Reach)/N                                             #Valor real


#Valores de los estimadores
Nh_MLE = getNh_MLE(encuesta, v_pob_totales)
Nh_MLE
Nh_MLEvis = getNh_MLEvis(encuesta, v_pob_totales, visibility_factor)
Nh_MLEvis
Nh_MLESsinvis = getNh_MLEsinvis(encuesta, v_pob_totales)
Nh_MLESsinvis


#Valor real del tamaño de la población oculta
sum(Poblacion$Poblacion_Oculta) 
