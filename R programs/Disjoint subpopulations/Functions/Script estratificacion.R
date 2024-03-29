library(igraph)
library(igraph)
library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)
library(sampler)
library(tibble)


#Asignamos a cada uno de los individuos una subpoblación (FUNCIÓN)
getPobDis <- function(k, pob, p) {
  # Generador de poblaciones disjuntas
  # Devuelve un vector con la poblaci?n correspondiente a cada individuo
  # k: n?mero de individuos
  # pob:  vector con enteros representando las distintas poblaciones
  # p: vector con enteros con las probabilidades de las subpoblaciones
  sample(pob, k, replace = TRUE, p = p)
}

#Asignamos uniformemente si pertenecen a la población oculta (FUNCIÓN)
getHiddenPop <- function(k, prob) {
  # Generador de la poblaci?n oculta
  # Devuelve un vector de ceros y unos, los unos representan los elementos de la poblaci?n oculta
  # prob: probabilidad de la poblaci?n oculta
  sample(c(0,1), k, replace = TRUE, p = c(1-prob,prob))
}

#Asignamos uniformemente si pertenecen a la población oculta (FUNCIÓN)
getHiddenPop <- function(k, prob) {
  # Generador de la poblaci?n oculta
  # Devuelve un vector de ceros y unos, los unos representan los elementos de la poblaci?n oculta
  # prob: probabilidad de la poblaci?n oculta
  sample(c(0,1), k, replace = TRUE, p = c(1-prob,prob))
}

# Asigna los valores del estrato en funci�n de la proporcion
getEstrato <- function(n,StratumProp){
  sample(1:length(StratumProp),n,replace = TRUE, p= StratumProp)
}

#Generamos las poblaciones
genPoblacion <- function(n, pob_dis, v_pob,prob_hidden) {
  # Genera un data frame con la poblaci?n y la pertenencia a la poblaci?n oculta
  enc = data.frame(Poblacion = getPobDis(n,pob_dis,v_pob))
  enc = cbind(enc, Poblacion_Oculta = getHiddenPop(n,prob_hidden))
  rownames(enc) <- c(1:n)
  return(enc)
}

######################
# Matriz para el GNSUM
######################

matrixHP = function(grafo,Pob){
  # Matriz de adyaciencia del grafo dirigido de conexiones con la HP
  ad = as_adj(grafo) # Matriz de adyacencia
  for (j in 1:ncol(ad)) {
    #if (V(grafo)$label[j] == 0){
    if (Pob$Poblacion_Oculta[j]==0){
      ad[,j] = 0
    }
  }
  return(ad)
}

#######################

berHP = function(x,p){
  if(x!=0){
    return(x*rbinom(1,1,p))
  }
  else {
    return(0)
  }
}

######################

###############################
# Generación de las poblaciones
###############################


getDatos_conEstrat = function(N, poblaciones_dis,v_probabilidades,PropHiddenPop,StratumProp, dim, nei, p, visibility_factor, memory_factor, sub_memory_factor){
  # Lista que contine el grafo con v?rtices etiquetados, y los datos de la poblaci?n
  
  #N es el tamaño de la población que queremos
  #poblaciones_dis vector con todas las poblaciones
  #v_probabilidades vector con las probabilidades de las poblaciones
  #PropHiddenPop proporción de la población oculta
  #StratumProp , vector con las proporciones del estrato en la Poblaci�n
  #dim Integer constant, the dimension of the starting lattice.
  #nei Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p Real constant between zero and one, the rewiring probability.
  #memory_factor es el factor de memoria
  #visibility_factor es el factor de visibilidad
  #memory_factor es la proporción que vamos a tomar la varianza de la normal que vamos a aplicar al Reach
  #sub_memory_factor es la proporción que vamos a tomar como varianza de la normal que vamos a aplicar a las subpoblaciones
  
  Pob_general = genPoblacion(N, poblaciones_dis,v_probabilidades,PropHiddenPop)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_poblaciones = length(poblaciones_dis)-1
  
  # Vector con la pertenencia a los diferentes estratos
  Estratos = getEstrato(N,StratumProp)
  # Muestra con reemplazamiento de los factores de correcci?n con probabilidades vect_prob_memory_subpob
  vect_hp = rep(NA,N)       #Vector sin visibility
  vect_hp_vis = rep(NA,N)   #Vector con visibility aplicado
  vect_reach = rep(NA,N)    #Vector del número de nodos que se conecta cada nodo
  vect_reach_re = rep(NA,N) #Vector Reach con error de memoria
  
  # Creamos la matriz del grafo dirigido de personas que conocen a la poblaci?n oculta
  Mhp = matrixHP(net_sw,Pob_general)
  Mhp_vis =  apply(Mhp,c(1,2), berHP,p = visibility_factor)
  
  
  for (i in 1:N) {
    # En la clase de igraph, [[]] busca vértices adyacentes
    # net_sw[[i]] es un lista con un elemento, la lista con los vértices adyacentes
    
    vect_hp[i] = sum(Mhp[i,])
    vect_reach[i] = length(net_sw[[i]][[1]])
    #  binomial de tamaño las conexiones con la población oculta y probabilidad la visibilidad
    #vect_hp_vis[i] = rbinom(1,vect_hp[i],prob = visibility_factor)
    vect_hp_vis[i] = sum(Mhp_vis[i,])
    
    vect_reach_re[i] = round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i]))
  }
  
  Pob_general = cbind(Pob_general, Estrato = Estratos)
  Pob_general = cbind(Pob_general, Reach = vect_reach)
  Pob_general = cbind(Pob_general, Reach_memory = vect_reach_re)
  Pob_general = cbind(Pob_general, HP_total_conocida = vect_hp) 
  Pob_general = cbind(Pob_general, HP_total_apvis = vect_hp_vis)
  
  for(j in 0:n_poblaciones){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Pob_general[net_sw[[i]][[1]],]$Poblacion == j)
      # Visibilidad de la poblaci?n j por i aplicando un factor de visibilidad para las subpoblaciones
      v_1[i] = round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob))
    }
    
    Pob_general = cbind(Pob_general,Subpoblacion_total = v_1)
    names(Pob_general)[dim(Pob_general)[2]] = str_c("KP_total_apvis_",j)
  }
  
  
  
  
  
  # A?adimos las etiquetas de la pertenencia a la Hidden Population al grafo
  #V(net_sw)$label = Pob_general$Poblacion_Oculta
  
  return(list(net_sw, Pob_general, Mhp_vis))
}


getEncuesta = function(n_enc, dataframe){
  #Se introduce un número n_enc que es el tamaño de la misma
  #Esta función genera aleatoriamente un dataframe con parte de las columnas del Pob_general
  encuesta = dataframe[sample(nrow(dataframe), n_enc, replace = FALSE),]
  encuesta
}

getMuestraEstrat = function(n_enc,dataframe){
  # Muestro estratificado por asignaci�n proporcional
  # ssamp(dataframe,n_enc,Estrato)
  dataframe %>% rownames_to_column('index') %>% ssamp(n_enc,Estrato) %>% column_to_rownames('index')
}
############################

##################
# Estimador b?sico
##################

getNh_basic = function(enc,N) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Es el estimador básico, luego no tenemos nada en cuenta
  #enc es la encuesta
  #N es el tamaño de la población
  Nh_f =  N*sum(enc$HP_total_apvis)/sum(enc$Reach_memory) 
  Nh_f
}

NSUM.basic.s = function(enc,N,StratumProp){
  Reach.medio = mean(enc$Reach_memory)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_apvis/Reach.medio), sd = sd(N*HP_total_apvis/Reach.medio), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  #varianza = sum((1/tabla$n-1/N)*StratumProp*tabla$S^2)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.basic.ps = function(enc,N,StratumProp){
  Reach.medio = mean(enc$Reach_memory)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_apvis/Reach.medio), sd = sd(N*HP_total_apvis/Reach.medio), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

getNh_basicvis = function(enc,N,vis) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Es el estimador básico, aplicando el factor corrector de la visibilidad
  #enc es la encuesta
  #N es el tamaño de la población TOTAL
  Nh_f =  N*sum(enc$HP_total_apvis)/sum(enc$Reach_memory) * (1/vis)
  Nh_f
}

NSUM.basicvis.s = function(enc,N,vis,StratumProp){
  Reach.medio = mean(enc$Reach_memory)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_apvis/Reach.medio/vis), sd = sd(N*HP_total_apvis/Reach.medio/vis), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.basicvis.ps = function(enc,N,vis,StratumProp){
  Reach.medio = mean(enc$Reach_memory)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_apvis/Reach.medio/vis), sd = sd(N*HP_total_apvis/Reach.medio/vis), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}
getNh_basicsinvis = function(enc,N) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Estimador básico si los datos no tuvieran un error de visibilidad
  #enc es la encuesta
  #N es el tamaño de la población TOTAL
  Nh_f =  N*sum(enc$HP_total_conocida)/sum(enc$Reach_memory)
  Nh_f
}

NSUM.basicsinvis.s = function(enc,N,StratumProp){
  Reach.medio = mean(enc$Reach_memory)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_conocida/Reach.medio), sd = sd(N*HP_total_conocida/Reach.medio), n = length(Poblacion) )
  #summarise(Nh = mean (N*HP_total_apvis/Reach.medio), S = sd(N*HP_total_apvis/Reach.medio))
  Nh = sum(tabla$Nh*StratumProp)
  #varianza = sum((1/tabla$n-1/N)*StratumProp*tabla$S^2)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.basicsinvis.ps = function(enc,N,StratumProp){
  Reach.medio = mean(enc$Reach_memory)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_conocida/Reach.medio), sd = sd(N*HP_total_conocida/Reach.medio), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

getNh_basicsinvismem = function(enc,N) {
  #Aplicando la fórmula del paper Thirty Years of NSUM obtenemos la estimacion de la encuesta realizada
  #Estimador básico si los datos no tuvieran un error de memoria ni de visibilidad
  #enc es la encuesta
  #N es el tamaño de la población TOTAL
  Nh_f =  N*sum(enc$HP_total_conocida)/sum(enc$Reach)
  Nh_f
}

NSUM.basicsinvismem.s = function(enc,N,StratumProp){
  Reach.medio = mean(enc$Reach)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_conocida/Reach.medio), sd = sd(N*HP_total_conocida/Reach.medio), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.basicsinvismem.ps = function(enc,N,StratumProp){
  Reach.medio = mean(enc$Reach)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*HP_total_conocida/Reach.medio), sd = sd(N*HP_total_conocida/Reach.medio), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

###############
# Estimador MLE
###############

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

NSUM.MLE.s = function(enc,v_pob,N,StratumProp){
  Reach.medio = sum(enc[tail(names(enc),length(v_pob))])/nrow(enc)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (sum(v_pob)*HP_total_apvis/Reach.medio), sd = sd(sum(v_pob)*HP_total_apvis/Reach.medio), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  #varianza = sum((1/tabla$n-1/N)*StratumProp*tabla$S^2)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.MLE.ps = function(enc,v_pob,N,StratumProp){
  Reach.medio = sum(enc[tail(names(enc),length(v_pob))])/nrow(enc)
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (sum(v_pob)*HP_total_apvis/Reach.medio), sd = sd(sum(v_pob)*HP_total_apvis/Reach.medio), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
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

NSUM.MLEvis.s = function(enc,v_pob,N,vis,StratumProp){
  Reach.medio = sum(enc[tail(names(enc),length(v_pob))])/nrow(enc)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (sum(v_pob)*HP_total_apvis/Reach.medio/vis), sd = sd(sum(v_pob)*HP_total_apvis/Reach.medio/vis), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.MLEvis.ps = function(enc,v_pob,N,vis,StratumProp){
  Reach.medio = sum(enc[tail(names(enc),length(v_pob))])/nrow(enc)
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (sum(v_pob)*HP_total_apvis/Reach.medio/vis), sd = sd(sum(v_pob)*HP_total_apvis/Reach.medio/vis), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
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

NSUM.MLEsinvis.s = function(enc,v_pob,N,StratumProp){
  Reach.medio = sum(enc[tail(names(enc),length(v_pob))])/nrow(enc)
  #enc$prob = enc$HP_total_apvis/Reach.medio
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (sum(v_pob)*HP_total_conocida/Reach.medio), sd = sd(sum(v_pob)*HP_total_conocida/Reach.medio), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.MLEsinvis.ps = function(enc,v_pob,N,StratumProp){
  Reach.medio = sum(enc[tail(names(enc),length(v_pob))])/nrow(enc)
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (sum(v_pob)*HP_total_conocida/Reach.medio), sd = sd(sum(v_pob)*HP_total_conocida/Reach.medio), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

################
#Estimador PIMLE
################
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

NSUM.PIMLE.s = function(enc,v_pob,N,StratumProp){
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  enc$prob = enc$HP_total_apvis/d_iest
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), sd = sd(N*prob), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.PIMLE.ps = function(enc,v_pob,N,StratumProp){
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  enc$prob = enc$HP_total_apvis/d_iest
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), sd = sd(N*prob), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
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

NSUM.PIMLEvis.s = function(enc,v_pob,N,vis,StratumProp){
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  enc$prob = enc$HP_total_apvis/d_iest
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob/vis), sd = sd(N*prob/vis), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.PIMLEvis.ps = function(enc,v_pob,N,vis,StratumProp){
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  enc$prob = enc$HP_total_apvis/d_iest
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob/vis), sd = sd(N*prob/vis), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
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

NSUM.PIMLEsinvis.s = function(enc,v_pob,N,StratumProp){
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  enc$prob = enc$HP_total_conocida/d_iest
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), sd = sd(N*prob), n = length(Poblacion) )
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.PIMLEsinvis.ps = function(enc,v_pob,N,StratumProp){
  d_iest = c()
  for (i in 1:nrow(enc)) {
    d_iest[i] = N * sum(enc[i,tail(names(enc),length(v_pob))])/sum(v_pob)
  }
  enc$prob = enc$HP_total_conocida/d_iest
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), sd = sd(N*prob), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}
###############
# Estimador MoS
###############

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

NSUM.MoS.s = function(enc,v_pob,N,StratumProp){
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  enc$prob = enc$HP_total_apvis/d_i_est
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), s = sd(N*prob), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$s^2)
  return(list(tabla,Nh,varianza))
}

NSUM.MoS.ps = function(enc,v_pob,N,StratumProp){
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  enc$prob = enc$HP_total_apvis/d_i_est
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), sd = sd(N*prob), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
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

NSUM.MoSvis.s = function(enc,v_pob,N,vis,StratumProp){
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  enc$prob = enc$HP_total_apvis/d_i_est
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob/vis), sd = sd(N*prob/vis), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.MoSvis.ps = function(enc,v_pob,N,vis,StratumProp){
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  enc$prob = enc$HP_total_apvis/d_i_est
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob/vis), sd = sd(N*prob/vis), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
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

NSUM.MoSsinvis.s = function(enc,v_pob,N,StratumProp){
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  enc$prob = enc$HP_total_conocida/d_i_est
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), sd = sd(N*prob), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

NSUM.MoSsinvis.ps = function(enc,v_pob,N,StratumProp){
  d_i_est = rep(NA, nrow(enc))
  for (i in 1:nrow(enc)) {
    d_i_est[i] = (sum((enc[i,tail(names(enc),length(v_pob))])/v_pob))/length(v_pob) * N
  }
  enc$prob = enc$HP_total_conocida/d_i_est
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (N*prob), sd = sd(N*prob), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

#######
# GNSUM
#######

getNh_GNSUM_antiguo = function(enc, enc_hp, v_pob, N, visibility_factor){
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


getNh_GNSUM = function(Pob, enc, enc_hp, Mhp_vis, v_pob, N){
  #enc es la encuesta realizada
  #enc_hp es la encuesta de la pob oculta
  #M, matriz de las conexiones de la visibilidad con la HP
  #N es el tamaño de la poblacion general
  #v_pob el vector de probabilidades de las poblaciones
  
  #Calculamos la estimador del numerador de la fórmula, correspondiente al número de personas que se conocen de la población oculta
  n_enc = nrow(enc)
  prob_inc = n_enc/N
  numerador = (1/prob_inc) * sum(enc$HP_total_apvis) #En nuestro caso la probabilidad de inclusión es la misma para todos los elementos, la sacamos del sumatorio
  
  
  #Calculamos la estimación del denominador
  #Para ello necesitamos una encuesta a la población oculta
  #enc[enc$Poblacion_Oculta ==1,] -> Usamos la encuesta que ya teníamos y entrevistamos a esos elementos
  #denominador = (N/sum(v_pob)) * rbinom(1,sum(enc_hp[,(ncol(enc_hp)-length(v_pob)+1):(ncol(enc_hp))-(length(v_pob)+1)]),visibility_factor)/nrow(enc_hp)
  ind1 = as.numeric(rownames(enc_hp))
  ind2 = as.numeric(rownames(Pob[Pob$Poblacion != 0,]))
  suma = sum(Mhp_vis[ind2,ind1])
  
  denominador = N/sum(v_pob)*suma/nrow(enc_hp)
  
  Nh = numerador/denominador
  return(Nh)
}

GNSUM.s = function(Pob, enc, enc_hp, Mhp_vis, v_pob, N, StratumProp){
  #Calculamos la estimacion del denominador
  ind1 = as.numeric(rownames(enc_hp))
  ind2 = as.numeric(rownames(Pob[Pob$Poblacion != 0,]))
  suma = sum(Mhp_vis[ind2,ind1])
  
  # Calculamos el denominador del GNSUM
  denominador = N/sum(v_pob)*suma/nrow(enc_hp)
  
  
  enc$prop =N* enc$HP_total_apvis/denominador
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (prop), sd = sd(prop), n = length(Poblacion))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}

GNSUM.ps = function(Pob, enc, enc_hp, Mhp_vis, v_pob, N, StratumProp){
  #Calculamos la estimacion del denominador
  ind1 = as.numeric(rownames(enc_hp))
  ind2 = as.numeric(rownames(Pob[Pob$Poblacion != 0,]))
  suma = sum(Mhp_vis[ind2,ind1])
  
  # Calculamos el denominador del GNSUM
  denominador = N/sum(v_pob)*suma/nrow(enc_hp)
  
 
  enc$prop =N* enc$HP_total_apvis/denominador
  tabla = enc %>%
    group_by(Estrato) %>%
    summarise(Nh = mean (prop), sd = sd(prop), n = length(Poblacion))
  #summarise(Nh = mean (N*HP_total_apvis/Reach.medio), S = sd(N*HP_total_apvis/Reach.medio))
  Nh = sum(tabla$Nh*StratumProp)
  varianza = (1/nrow(enc)-1/N)*sum(StratumProp*tabla$sd^2)+1/nrow(enc)^2*sum((1-StratumProp)*tabla$sd^2)
  return(list(tabla,Nh,varianza))
}



#########
# Directo
#########

getNh_Directo = function(enc,N){
  Nh = sum(enc$Poblacion_Oculta)/nrow(enc) * N
  return(Nh)
}





