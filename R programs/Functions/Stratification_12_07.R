getStratum <- function(n,StratumProp){
  # Generator of stratum in the population
  # n: the number of individuals
  # StratumProp: vector with the subpopulations probabilities
  sample(1:length(StratumProp),n,replace = TRUE, p= StratumProp)
}


getDatos_conEstrat = function(N, poblaciones_dis,v_probabilidades,PropHiddenPop,StratumProp, dim, nei, p, visibility_factor, memory_factor, sub_memory_factor){
  # Lista que contine el grafo con v?rtices etiquetados, y los datos de la poblaci?n
  
  #N es el tamaño de la población que queremos
  #poblaciones_dis vector con todas las poblaciones
  #v_probabilidades vector con las probabilidades de las poblaciones
  #PropHiddenPop proporción de la población oculta
  #StratumProp , vector con las proporciones totales del estrato
  #dim Integer constant, the dimension of the starting lattice.
  #nei Integer constant, the neighborhood within which the vertices of the lattice will be connected.
  #p Real constant between zero and one, the rewiring probability.
  #memory_factor es el factor de memoria
  #visibility_factor es el factor de visibilidad
  #memory_factor es la proporción que vamos a tomar la varianza de la normal que vamos a aplicar al Reach
  #sub_memory_factor es la proporción que vamos a tomar como varianza de la normal que vamos a aplicar a las subpoblaciones
  
  Pob_general = genPopulation(N, poblaciones_dis,v_probabilidades,PropHiddenPop)
  
  net_sw = sample_smallworld(dim, N, nei, p, loops = FALSE, multiple = FALSE)
  
  n_populations = length(poblaciones_dis)-1
  
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
    
    vect_reach_re[i] = max(0,round(rnorm(1, mean = vect_reach[i], sd = memory_factor*vect_reach[i])))
  }
  
  Pob_general = cbind(Pob_general, Estrato = Estratos)
  Pob_general = cbind(Pob_general, Reach = vect_reach)
  Pob_general = cbind(Pob_general, Reach_memory = vect_reach_re)
  Pob_general = cbind(Pob_general, HP_total_conocida = vect_hp) 
  Pob_general = cbind(Pob_general, HP_total_apvis = vect_hp_vis)
  
  for(j in 0:n_populations){
    v_1 = rep(NA,N)
    for(i in 1:N) {
      vis_pob = sum(Pob_general[net_sw[[i]][[1]],]$Poblacion == j)
      # Visibilidad de la poblaci?n j por i aplicando un factor de visibilidad para las subpoblaciones
      v_1[i] = max(0,round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob)))
    }
    
    Pob_general = cbind(Pob_general,Subpoblacion_total = v_1)
    names(Pob_general)[dim(Pob_general)[2]] = str_c("KP_total_apvis_",j)
  }
  
  
  
  
  
  # A?adimos las etiquetas de la pertenencia a la Hidden Population al grafo
  #V(net_sw)$label = Pob_general$Poblacion_Oculta
  
  return(list(net_sw, Pob_general, Mhp_vis))
}

getMuestraEstrat = function(n_enc,dataframe, Estrato){
  # Muestro estratificado por asignaci?n proporcional
  ssamp(dataframe,n_enc,Estrato)
}