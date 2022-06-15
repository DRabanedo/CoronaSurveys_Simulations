# Prueba GNSUM
library(igraph)
library(igraph)
library(tidyverse)
library(stringr)
library(ggplot2)

######################
getNh_GNSUM = function(enc, enc_hp, v_pob, N, visibility_factor){
  #enc es la encuesta realizada
  #N es el tamaño de la poblacion general
  #v_pob el vector de probabilidades de las poblaciones
  
  #Calculamos la estimador del numerador de la fórmula, correspondiente al número de personas que se conocen de la población oculta
  n_enc = nrow(enc)
  prob_inc = n_enc/N
  numerador = (1/prob_inc) * sum(enc$HP_total_apvis) #En nuestro caso la probabilidad de inclusión es la misma para todos los elementos, la sacamos del sumatorio
  
  
  #Calculamos la estimación del denominador
  #Para ello necesitamos una encuesta a la población oculta
  #enc[enc$Poblacion_Oculta ==1,] -> Usamos la encuesta que ya teníamos y entrevistamos a esos elementos
  denominador = (N/sum(v_pob)) * rbinom(1,sum(enc_hp[,(ncol(enc_hp)-length(v_pob)+1):(ncol(enc_hp))-(length(v_pob)+1)]),visibility_factor)/nrow(enc_hp)
  
  Nh = numerador/denominador
  return(Nh)
}

#############################
N = 10000 #Tamaño de la población
poblaciones_dis = c(0:10) #Poblaciones a distinguir. Son disjuntas y el 0 corresponde a no clasificar al individuo en ninguna
v_probabilidades = c(0.3, 0.1,0.05,0.005,0.005,0.04, 0.2, 0.1, 0.15, 0.025, 0.025) #Probabilidad de cada una de las poblaciones anteriores
PropHiddenPop = 0.1 #Probabilidad de un individuo de pertenecer a la población oculta
n_encuesta = 300 #Número de muestras que extraemos en la encuesta
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
encuesta = getEncuesta(300,Poblacion)
encuesta_hp = getEncuesta(50, Poblacion[Poblacion$Poblacion_Oculta==1,])

v_pob_totales = rep(NA, n_poblaciones)
for (j in 1:n_poblaciones) {
  v_pob_totales[j] = sum(Poblacion$Poblacion == j)
}

getNh_GNSUM(encuesta, encuesta_hp, v_pob_totales, N, visibility_factor)
sum(Poblacion$Poblacion_Oculta)
