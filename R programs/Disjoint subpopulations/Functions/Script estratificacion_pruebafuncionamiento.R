#Simulaciones de prueba de los estimadores con estratificación.

N = 1000                  #Tamaño de la población general
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
StratumProp = c(1/5,1/5,1/10,1/10,1/5,1/5)

semilla = 207              #Semilla utilizada para la simulación 
set.seed(semilla)

dim = 1     #Dimensión del grafo
nei = 75    #Número de vecinos con los que vamos a conectar. Son vecinos a cada lado del nodo, luego son 2*nei
p   = 0.1   #Probabilidad de aleatorizar una conexión. Cuanto más grande más aleatorio


#Generamos la población sobre la que vamos a realizar el estudio
Datos = getDatos_conEstrat(N, poblaciones_dis,v_probabilidades,PropHiddenPop,StratumProp, dim, nei, p, visibility_factor, memory_factor, sub_memory_factor)

net_sw = Datos[[1]]
Poblacion = Datos[[2]]
Mhp_vis = Datos[[3]]


#Sacamos las encuestas
encuesta = getEncuesta(n_encuesta, Poblacion)
encuesta_s = getMuestraEstrat(n_encuesta,Poblacion)
encuesta_hp = getEncuesta(n_encuesta_hp, Poblacion[Poblacion$Poblacion_Oculta==1,])

#V_pob
v_pob = rep(NA, n_poblaciones)
for (j in 1:n_poblaciones) {
  v_pob[j] = sum(Poblacion$Poblacion == j)
}


#Estimadores con sus respectivas estimaciones

#Básico
getNh_basic(encuesta,N)

NSUM.basic.s(encuesta_s,N,StratumProp)[[2]]
NSUM.basic.ps(encuesta,N,StratumProp)[[2]]

getNh_basicvis(encuesta,N,visibility_factor)

NSUM.basicvis.s(encuesta_s,N,visibility_factor,StratumProp)[[2]]
NSUM.basicvis.ps(encuesta,N,visibility_factor,StratumProp)[[2]]

#MLE
getNh_MLE(encuesta,v_pob)

NSUM.MLE.s(encuesta_s,v_pob,N,StratumProp)[[2]]
NSUM.MLE.ps(encuesta,v_pob,N,StratumProp)[[2]]

getNh_MLEvis(encuesta,v_pob,visibility_factor)

NSUM.MLEvis.s(encuesta_s,v_pob,N,visibility_factor,StratumProp)[[2]]
NSUM.MLEvis.ps(encuesta,v_pob,N,visibility_factor,StratumProp)[[2]]

#PIMLE
getNh_PIMLE(encuesta,v_pob,N)

NSUM.PIMLE.s(encuesta_s,v_pob,N,StratumProp)[[2]]
NSUM.PIMLE.ps(encuesta,v_pob,N,StratumProp)[[2]]

getNh_PIMLEvis(encuesta,v_pob,N,visibility_factor)

NSUM.PIMLEvis.s(encuesta_s,v_pob,N,visibility_factor,StratumProp)[[2]]
NSUM.PIMLEvis.ps(encuesta,v_pob,N,visibility_factor,StratumProp)[[2]]

#MoS
getNh_MoS(encuesta,v_pob,N)

NSUM.MoS.s(encuesta_s,v_pob,N,StratumProp)[[2]]
NSUM.MoS.ps(encuesta,v_pob,N,StratumProp)[[2]]

getNh_MoSvis(encuesta,v_pob,N,visibility_factor)

NSUM.MoSvis.s(encuesta_s,v_pob,N,visibility_factor,StratumProp)[[2]]
NSUM.MoSvis.ps(encuesta,v_pob,N,visibility_factor,StratumProp)[[2]]

#GNSUM

getNh_GNSUM(Poblacion, encuesta, encuesta_hp, Mhp_vis, v_pob, N)

GNSUM.s(Poblacion, encuesta_s, encuesta_hp, Mhp_vis, v_pob, N, StratumProp)[[2]]
GNSUM.ps(Poblacion, encuesta, encuesta_hp, Mhp_vis, v_pob, N, StratumProp)[[2]]





