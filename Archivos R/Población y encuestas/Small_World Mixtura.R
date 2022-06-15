library(igraph)
library(tidyverse)
library(stringr)
#Este programa va a asignar con probabilidad uniforme una clase a cada individuo.
#En este caso trabajaremos con poblaciones disjuntas, como podría ser clasificar 
#a la población en distitos estratos en función de su nombre. También vamos a utilizar
#una probabilidad uniforme para decidir si un individuo pertenece a la población
#oculta o no. De esta forma podemos trabajar con los datos fuera de la red de 
#contactos y utilizar la red generada sólo para ver que individuos están conectados.



#Número de estratos
n_poblaciones = 10
semilla = 207
set.seed(semilla)
#Vector con los nombres de las poblaciones, uso el 0 como no pertenecer a 
#ninguna subpoblación de las descritas
poblaciones_dis = c(0,1:n_poblaciones) 
#Vector con las probabilidades de cada subpoblación conocida
v_probabilidades = c(rep(1/(1+n_poblaciones), 1+n_poblaciones)) 
#Número de individuos
N = 5000
#Probabilidad de pertenecer a la población oculta
PropHiddenPop = 0.05



#Definimos los factores que van a influir en nuestra predicción y que hay que
#estimar de alguna forma
visibility_factor = 0.6
#barrier_effect = 0.95 (No creo que exista en este caso)
#frame_ratio = 0.95    (Si arreglamos post-strat todo bien)
memory_factor = 0.8
par_poison = 20
vect_memory_subpob = c(0.8,0.9,1,1.1,1.2)
vect_prob_memory_subpob = c(0.1,0.2,0.4,0.2,0.1)

memory_factor_subpob = sample(vect_memory_subpob, N*(n_poblaciones+1), replace = TRUE, vect_prob_memory_subpob)
#Aplicar un factor modificador que infraestime/sobreestime la 
#                variable Reach para obtener esa variabilidad que existe al no 
#                conocer cada persona exactamente su red de contactos. Yo la 
#                introduciría en Reach pero no en el estudio de las subpoblaciones
#                dado que ahí las personas son más precisas al decir el número de
#                personas que conocen, para así darle sentido al MLE y MoS.



#Asignamos a cada uno de los individuos una subpoblación (FUNCIÓN)
getPobDis <- function(k, pob = poblaciones_dis, p = v_probabilidades) {
  # Generador de poblaciones disjuntas
  # Devuelve un vector con la poblaci?n correspondiente a cada individuo
  # k: n?mero de individuos
  # pob:  vector con enteros representando las distintas poblaciones
  # p: vector con enteros con las probabilidades de las subpoblaciones
  sample(pob, k, replace = TRUE, p = p)
}

# Ejemplo 
getPobDis(100)
getPobDis(100,c(0,1:5),c(rep(1/(1+5), 1+5)))




#Asignamos uniformemente si pertenecen a la población oculta (FUNCIÓN)
getHiddenPop <- function(k, prob = PropHiddenPop) {
  # Generador de la poblaci?n oculta
  # Devuelve un vector de ceros y unos, los unos representan los elementos de la poblaci?n oculta
  # prob: probabilidad de la poblaci?n oculta
  sample(c(0,1), k, replace = TRUE, p = c(1-prob,prob))
}
getHiddenPop(100)




#Generamos una población a partir de los datos que queremos/necesitamos (FUNCIÓN)
genPoblacion <- function(n, pob_dis, v_pob,prob_hidden) {
  # Genera un data frame con la poblaci?n y la pertenencia a la poblaci?n oculta
  enc = data.frame(Poblacion = getPobDis(n,pob_dis,v_pob))
  enc = cbind(enc, Poblacion_Oculta = getHiddenPop(n,prob_hidden))
  rownames(enc) <- c(1:n)
  return(enc)
} 



#Generamos la población de estudio
Pob_general = genPoblacion(N, poblaciones_dis,v_probabilidades,PropHiddenPop)



#Determinamos la distribución de nuestra población, en este caso usando una
#network basada en el Small-World. Habría que buscar cúal se adapta mejor a 
#la realidad, variando el modelo o los parámetros.
net_sw = sample_smallworld(dim = 1,size =  N, nei = 50, p = 0.34, loops = FALSE, multiple = FALSE)



#Cuanta gente conoce cada individuo de la población oculta. Para tener en cuenta
#el factor de visibilidad se le puede aplicar directamente sobre la suma final, de 
#forma que siempre sea constante o aplicar en cada uno de los individuos. He implementado 
#lo segundo, ya que me ha parecido lo más indicado, manteniendo asi la aleatoriedad.

# Simulamos las conexiones con la poblaci?n oculta, las conexiones con la poblaci?n oculta
# que realmente se conocen y el grado de cada nodo

#Vector sin visibility
vect_hp = rep(NA,N) 
#Vector con visibility aplicado
vect_hp_vis = rep(NA,N)
#Vector del número de nodos que se conecta cada nodo
vect_reach = rep(NA,N)
#Vector Reach con error de memoria
vect_reach_re = rep(NA,N)
#Vector de decisiones de la mixtura
vect_mix = sample(c(1,0), N*n_poblaciones ,replace = TRUE, p = c(0.5,0.5))
#Vector de lo que añadimos en el caso de sobreestimar
vect_poison = rpois(N, par_poison)
for (i in 1:N) {
  # En la clase de igraph, [[]] busca vértices adyacentes
  # net_sw[[i]] es un lista con un elemento, la lista con los vértices adyacentes
  
  vect_hp[i] = sum(Pob_general[net_sw[[i]][[1]],]$Poblacion_Oculta)
  vect_reach[i] = length(net_sw[[i]][[1]])
  #  binomial de tamaño las conexiones con la población oculta y probabilidad la visibilidad
  vect_hp_vis[i] = rbinom(1,vect_hp[i],prob = visibility_factor)
  
 
  if (vect_mix[i] == 1){
    #  binomial de tamaño las conexiones con la población  y probabilidad el factor de memoria
    vect_reach_re[i] = rbinom(1,vect_reach[i],prob = memory_factor)
  } else {
    #  poison para simular la sobreestimación
    vect_reach_re[i] = vect_reach[i] + vect_poison[i]
  }
  
}
 
Pob_general = cbind(Pob_general, Reach = vect_reach)
Pob_general = cbind(Pob_general, Reach_memory = vect_reach_re)
Pob_general = cbind(Pob_general, HP_total_conocida = vect_hp) 
Pob_general = cbind(Pob_general, HP_total_apvis = vect_hp_vis)


#Ahora vamos a encontrar el número de personas en cada subpoblación que conoce
#el individuo. No me dejaba emplear nombres dinámicos dadas las resticiones al 
#definir las columnas del dataframe, lo he hecho un poco a mano con un vector de
#nombres aprovechando el bucle e implementandolo despúes.
for(j in 0:n_poblaciones){
  v_1 = rep(NA,N)
  for(i in 1:N) {
    # Visibilidad de la poblaci?n j por i aplicando un factor de visibilidad para las subpoblaciones
    v_1[i] = round(sum(Pob_general[net_sw[[i]][[1]],]$Poblacion == j)*memory_factor_subpob[(n_poblaciones+1)*(i-1)+j+1])
  }
  
  Pob_general = cbind(Pob_general,Subpoblacion_total = v_1)
  names(Pob_general)[dim(Pob_general)[2]] = str_c("KP_total_",j)
}

getPoblacionTotal = function(N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor,vect_memory_subpob,vect_prob_memory_subpob,par_poison){
  #Esta funcion genera una poblacion de estudio para realizar las encuestas previstas
  #Genera un dataframe con la población
  #N es el tamaño de la población que queremos
  #poblaciones_dis vector con todas las poblaciones
  #v_probabilidades vector con las probabilidades de las poblaciones
  #PropHiddenPop proporción de la población oculta
  #dim Integer constant, the dimension of the starting lattice.
  #nei Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p Real constant between zero and one, the rewiring probability.
  #memory_factor es el factor de memoria
  #visibility_factor es l factor de visibilidad
  Pob_general = genPoblacion(N, poblaciones_dis,v_probabilidades,PropHiddenPop)
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  n_poblaciones = length(poblaciones_dis)-1
  memory_factor_subpob = sample(vect_memory_subpob, N*(n_poblaciones+1), replace = TRUE, vect_prob_memory_subpob)
  vect_hp = rep(NA,N)       #Vector sin visibility
  vect_hp_vis = rep(NA,N)   #Vector con visibility aplicado
  vect_reach = rep(NA,N)    #Vector del número de nodos que se conecta cada nodo
  vect_reach_re = rep(NA,N) #Vector Reach con error de memoria
  vect_poison = rpois(N, par_poison) #Vector poison para la sobre estimación
  vect_mix = sample(c(1,0), N*n_poblaciones ,replace = TRUE, p = c(0.5,0.5)) #Vector de decisiones de la mixtura
  for (i in 1:N) {
    # En la clase de igraph, [[]] busca vértices adyacentes
    # net_sw[[i]] es un lista con un elemento, la lista con los vértices adyacentes
    
    vect_hp[i] = sum(Pob_general[net_sw[[i]][[1]],]$Poblacion_Oculta)
    vect_reach[i] = length(net_sw[[i]][[1]])
    #  binomial de tamaño las conexiones con la población oculta y probabilidad la visibilidad
    vect_hp_vis[i] = rbinom(1,vect_hp[i],prob = visibility_factor)
    
    
    if (vect_mix[i] == 1){
      #  binomial de tamaño las conexiones con la población  y probabilidad el factor de memoria
      vect_reach_re[i] = rbinom(1,vect_reach[i],prob = memory_factor)
    } else {
      #  poison para simular la sobreestimación
      vect_reach_re[i] = vect_reach[i] + vect_poison[i]
    }
    
  }
  
  Pob_general = cbind(Pob_general, Reach = vect_reach)
  Pob_general = cbind(Pob_general, Reach_memory = vect_reach_re)
  Pob_general = cbind(Pob_general, HP_total_conocida = vect_hp) 
  Pob_general = cbind(Pob_general, HP_total_apvis = vect_hp_vis)
  
  for(j in 0:n_poblaciones){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      # Visibilidad de la poblaci?n j por i aplicando un factor de visibilidad para las subpoblaciones
      v_1[i] = round(sum(Pob_general[net_sw[[i]][[1]],]$Poblacion == j)*memory_factor_subpob[(n_poblaciones+1)*(i-1)+j+1])
    }
    
    Pob_general = cbind(Pob_general,Subpoblacion_total = v_1)
    names(Pob_general)[dim(Pob_general)[2]] = str_c("KP_total_",j)
  }

  return(Pob_general)
}


#Falta introducir el resto de los coefs y elegir la forma más adecuada de hacer el grafo

#Problemas:
#1. La variable Reach tiene una varianza muy pequeña mientras que en los datos 
#    reales encontramos una variación bastante fuerte. Esto supongo que podrá 
#    cambiarse modificando el tamaño/forma del grafo, pero no he encontrado una forma
#    de hacerlo con el Small World

#2. Hay que implementar todas las pequeñas perturbaciones generadas por los dife-
#    rentes sesgos, ya sea en las propias conexiones o en el total, yo opino que es 
#    mejor en las conexiones personales dado que dan una aleatoriedad más real, pero
#    son más difíciles de implementar 
#--> He implementado Visibility, faltaría ver como implementar el error de memoria
#    de una forma óptima/correcta para que lo represente bien.

#3. Mi código es un desastre, pero funciona. Cualquier mejora es bien recibida.




rpois(19,20)


#Hacer la mixtura de distribuciones