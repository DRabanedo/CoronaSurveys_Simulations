###############################
# Plug-in-MLE estimator (PIMLE)
###############################

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


getNh_PIMLE = function(enc,v_pob,N) {
  #enc es la encuesta, ?ltimos valores para la Known Population
  #v_pob es el vector con el tamaño de cada subpoblación
  #N es el tamaño de la población TOTAL
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #En esta fórmula tenemos en cuenta el error de memoria pero no el de visibilidad
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  Nh_f = N * mean(enc$HP_total_apvis/d_iest)
  Nh_f
}

getNh_PIMLEvis = function(enc,v_pob,N,vis) {
  #enc es la encuesta, ?ltimos valores para la Known Population
  #v_pob es el vector con el tamaño de cada subpoblación
  #n es el tamaño de la población TOTAL
  #vis es el factor de visibilidad
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #En esta fórmula tenemos en cuenta el error de memoria y el de visibilidad, y lo corregimos.
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  Nh_f = N * mean(enc$HP_total_apvis/d_iest) * (1/vis) # media de N \frac{y_{iu}}{\hat{d_i}}
  Nh_f
}

getNh_PIMLEsinvis =  function(enc,v_pob,n) {
  #enc es la encuesta, ?ltimos valores para la Known Population
  #v_pob es el vector con el tamaño de cada subpoblación
  #n es el tamaño de la población TOTAL
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #En este caso no hay error de visibilidad, y tenemos en cuenta el error de memoria si lo hubiera
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  Nh_f = N * mean(enc$HP_total_conocida/d_iest)
  Nh_f
}

################################################################################

#Valores de los estimadores
Nh_PIMLE = getNh_PIMLE(encuesta, v_pob_totales, N)
Nh_PIMLE
Nh_PIMLEvis = getNh_PIMLEvis(encuesta, v_pob_totales, N, visibility_factor)
Nh_PIMLEvis
Nh_PIMLESsinvis = getNh_PIMLEsinvis(encuesta, v_pob_totales, N)
Nh_PIMLESsinvis

#Valor real del tamaño de la población oculta
sum(Poblacion$Poblacion_Oculta) 


