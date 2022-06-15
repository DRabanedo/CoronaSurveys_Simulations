###########################
# Estimador básico del NSUM
###########################

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

################################################################################

getNh_basic = function(enc,N) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Es el estimador básico, luego no tenemos nada en cuenta
  #enc es la encuesta
  #N es el tamaño de la población
  Nh_f =  N*sum(enc$HP_total_apvis)/sum(enc$Reach_memory) 
  Nh_f
}

getNh_basicvis = function(enc,N,vis) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Es el estimador básico, aplicando el factor corrector de la visibilidad
  #enc es la encuesta
  #N es el tamaño de la población TOTAL
  Nh_f =  N*sum(enc$HP_total_apvis)/sum(enc$Reach_memory) * (1/vis)
  Nh_f
}
getNh_basicsinvis = function(enc,N) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Estimador básico si los datos no tuvieran un error de visibilidad
  #enc es la encuesta
  #N es el tamaño de la población TOTAL
  Nh_f =  N*sum(enc$HP_total_conocida)/sum(enc$Reach_memory)
  Nh_f
}

getNh_basicsinvismem = function(enc,N) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Estimador básico si los datos no tuvieran un error de memoria ni de visibilidad
  #enc es la encuesta
  #N es el tamaño de la población TOTAL
  Nh_f =  N*sum(enc$HP_total_conocida)/sum(enc$Reach)
  Nh_f
}

################################################################################

#Valor de los estimadores
Nh_basic =getNh_basic(encuesta,N) 
Nh_basic
Nh_basicvis =getNh_basicvis(encuesta,N,visibility_factor) 
Nh_basicvis
Nh_basicsinvis = getNh_basicsinvis(encuesta,N)
Nh_basicsinvis
Nh_basicsinvismem = getNh_basicsinvismem(encuesta,N)
Nh_basicsinvismem

#Valor real del tamaño de la población oculta
sum(Poblacion$Poblacion_Oculta) 
