######################################################################################
# Simulación del número subpoblaciones con mclapply, dejando el resto de cosas fijas
######################################################################################

#Vector de iteraciones
iteraciones = round(seq(from = 1, to = 20, length.out = 11))


#Definimos la semilla
semilla = 207            
set.seed(semilla)



  
#Introducimos los datos
  
N = 10000                  #Tamaño de la población general
poblaciones_dis = c(0:10)  #Vector de poblaciones disjuntas
#Vector de probabilidades de pertenencia a cada una de las poblaciones
v_probabilidades = rep(1/length(poblaciones_dis), length(poblaciones_dis)) 
PropHiddenPop = 0.1        #Probabilidad de pertenecer a la población oculta
n_encuesta = 500           #Tamaño de la encuesta a la población general
n_encuesta_hp = 50         #Tamaño de la encuesta a la población oculta
#Número de poblaciones utilizadas para la simulación
n_poblaciones =  length(poblaciones_dis)-1 
sub_memory_factor = 0      #Factor de memoria de las subpoblaciones
memory_factor = 0          #Factor memoria 
visibility_factor = 1      #Factor visibilidad
  
dim = 1     #Dimensión del grafo
nei = 75    #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1   #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio
  
numerosubpoblaciones.simulaciones <- function(l,N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor,
                                 memory_factor,sub_memory_factor,n_encuesta,n_encuesta_hp,iteraciones) {  
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
  
  
  
  #Generamos la población y la matriz 
  Grafo_poblacion_matriz = getDatos(N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor, memory_factor,sub_memory_factor)
  
  net_sw = Grafo_poblacion_matriz[[1]]
  Poblacion = Grafo_poblacion_matriz[[2]]
  Mhp = Grafo_poblacion_matriz[[3]]
  
  #Variable para poder seleccionar las componentes adecuadas
  k = length(poblaciones_dis)
  
  
  for (i in 1:length(iteraciones)) {
    print(i)
    m_pob = iteraciones[i]
    n_columnas = ncol(Poblacion)
    poblaciones_dis = c(0:m_pob)
    n_poblaciones = length(poblaciones_dis)-1 
    v_probabilidades = rep(1/length(poblaciones_dis), length(poblaciones_dis))
    
    Poblacion$Poblacion = sample(poblaciones_dis, N, replace = TRUE, p = v_probabilidades)
    Poblacion = Poblacion[,1:(ncol(Poblacion)-k)]
    
    for(j in 0:n_poblaciones){
      v_1 = rep(NA,N)
      for(v in 1:N) {
        vis_pob = sum(Poblacion[net_sw[[v]][[1]],]$Poblacion == j)
        # Visibilidad de la poblaci?n j por i aplicando un factor de visibilidad para las subpoblaciones
        v_1[v] = round(rnorm(1, mean = vis_pob, sd = sub_memory_factor*vis_pob))
      }
      
      Poblacion = cbind(Poblacion,Subpoblacion_total = v_1)
      names(Poblacion)[dim(Poblacion)[2]] = str_c("KP_total_apvis_",j)
    }
    
    k = length(poblaciones_dis)
    
    v_pob_totales = rep(NA, n_poblaciones)
    for (j in 1:n_poblaciones) {
      v_pob_totales[j] = sum(Poblacion$Poblacion == j)
    }
    
    #Hacemos las encuestas
    encuesta = getEncuesta(n_encuesta,Poblacion)
    encuesta_hp = getEncuesta(n_encuesta_hp, Poblacion[Poblacion$Poblacion_Oculta==1,])
    
    
    #Estimamos la poblacion oculta
    Nh_real[i] = sum(Poblacion$Poblacion_Oculta) 
    
    Nh_basic[i] = getNh_basic(encuesta,N) 
    Nh_basicvis[i] = getNh_basicvis(encuesta,N,visibility_factor) 
    
    Nh_PIMLE[i] = getNh_PIMLE(encuesta, v_pob_totales, N)
    Nh_PIMLEvis[i] = getNh_PIMLEvis(encuesta, v_pob_totales, N, visibility_factor)
    
    Nh_MLE[i] = getNh_MLE(encuesta, v_pob_totales)
    Nh_MLEvis[i] = getNh_MLEvis(encuesta, v_pob_totales, visibility_factor)
    
    Nh_MoS[i] = getNh_MoS(encuesta, v_pob_totales, N)
    Nh_MoSvis[i] = getNh_MoSvis(encuesta, v_pob_totales, N, visibility_factor)
    
    Nh_GNSUM[i] =  getNh_GNSUM(Poblacion, encuesta, encuesta_hp, Mhp, v_pob_totales, N)
  }
  
  
  #simulaciones = cbind(simulaciones,Nh_real = Nh_real)
  simulaciones = data.frame(Nh_real = Nh_real)
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
  return(simulaciones)
  
}

start_time = Sys.time()
lista_simulaciones =mclapply(1:25,numerosubpoblaciones.simulaciones,N, poblaciones_dis,v_probabilidades,PropHiddenPop, dim, nei, p, visibility_factor,
                             memory_factor,sub_memory_factor,n_encuesta,n_encuesta_hp,iteraciones)
end_time = Sys.time()
print(end_time-start_time)

simulaciones = data.frame(datos =  1:length(iteraciones))
for (i in lista_simulaciones) {
  simulaciones = cbind(simulaciones,i)
}

simulaciones
write.csv(simulaciones,                    # Data frame a ser exportado
          file = "Simulaciones_numerosubpoblaciones", # Nombre de la hoja de Excel
          row.names = TRUE )    # Incluir los nombres de las filas (TRUE) o no (FALSE)) 