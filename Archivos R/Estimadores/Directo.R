#Estimación directa de la población oculta

getNh_Directo = function(enc,N){
  Nh = sum(enc$Poblacion_Oculta)/nrow(enc) * N
  return(Nh)
}

################################################################################
N = 1000 #Tamaño de la población
poblaciones_dis = c(0:10) #Poblaciones a distinguir. Son disjuntas y el 0 corresponde a no clasificar al individuo en ninguna
v_probabilidades = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probabilidad de cada una de las poblaciones anteriores
PropHiddenPop = 0.1 #Probabilidad de un individuo de pertenecer a la población oculta
n_encuesta = 600 #Número de muestras que extraemos en la encuesta
n_poblaciones = length(poblaciones_dis)-1 #Número de poblaciones en las que clasificamos
par_poison = 20  #Parámetro de la Poisson con la que vamos a sobreestimar la variable Reach
vect_memory_subpob = c(0.8,0.9,1,1.1,1.2) #Vector del factor memoria que vamos a aplicar a las subpoblaciones
vect_prob_memory_subpob = c(0.1,0.2,0.4,0.2,0.1) #Probabilidad de aplicar cada una de las componentes del vector
memory_factor = 0.6      #Factor de memoria 
visibility_factor = 0.7  #Factor visibilidad
semilla = 207            #Está guay cambiarla pq aunque se mantenga la estructura, los resultados cambian

#Grafo
dim = 1    #Dimensión del grafo
nei = 52   #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1  #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio


Poblacion = getPoblacionTotal(N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor,vect_memory_subpob, vect_prob_memory_subpob,par_poison)
encuesta = getEncuesta(500,Poblacion)
encuesta_hp = getEncuesta(50, Poblacion[Poblacion$Poblacion_Oculta==1,])
sum(Poblacion$Poblacion_Oculta)

getNh_Directo(encuesta,N)

################################################################################


