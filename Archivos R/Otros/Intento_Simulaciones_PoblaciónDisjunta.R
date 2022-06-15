n_poblaciones = 30
semilla = 207
poblaciones = c(1:n_poblaciones)
poblaciones_dis = c(0,1:n_poblaciones)
v_probabilidades = c(rep(1/(1+n_poblaciones), 1+n_poblaciones))
N = 5000
PropHiddenPop = 0.2

#Definimos los factores que van a influir en nuestra predicción y que hay que
#estimar de alguna forma.
visibility_factor = 0.85
barrier_effect = 0.95
frame_ratio = 0.95
error_transmision = 0.90


# Si cada uno tiene que tener una y son disjuntas es fácil, incluso se le puede
# dar el valor a la última población como ninguna población, sería el caso de
# utilizar nombres y sólo nombres, o profesiones y sólo profesiones

getPobDis <- function(k, semilla) {
  set.seed(semilla)
  sample(poblaciones_dis, k, replace = TRUE, p = v_probabilidades)
}


# Asignamos una variable Reach a cada uno de ellos
getReach <- function(k, semilla) {
  set.seed(semilla)
  round(rexp(5000, rate = 1/105),0) #También usan 160
}

# Asignamos si pertenece a la población oculta o no cierta persona de la población
getHiddenPop <- function(k, semilla) {
  set.seed(semilla)
  sample(c(0,1), k, replace = TRUE, p = c(1-PropHiddenPop,PropHiddenPop))
}

#Sacamos las personas de la red de contactos, a partir de las cuales podremos 
#contar para sacar nuestras subpoblaciones

getReach2 <- function(n,r,semilla) {
  red <- list() #hay que hacer que acepte vectores como componentes
  set.seed(semilla)
  for (i in 1:5000) {
    ind_reach = r_ej[i]
    red_ind = sample(1:N, ind_reach, replace = FALSE)
    red[[i]] = red_ind
  }
  red
}

# ARD Enriquecido, nos dice cuanta gente conoce la gente de la subpoblación, lo
# que nos va a permitir estimar el factor de visibilidad y usar el estimador
# GNSUM, de forma que 

#cbind mete por columna y rbind por fila.
genEncuesta <- function(n, semilla) {
  encuesta <- data.frame(Reach = getReach(n,semilla))
  encuesta <- cbind(encuesta, Poblacion = getPobDis(n,semilla))
  encuesta <- cbind(encuesta, Poblacion_Oculta = getHiddenPop(n,semilla))
  rownames(encuesta) <- c(1:n)
  return(encuesta)
}
encuesta = genEncuesta(N,semilla)
encuesta_pob = list(encuesta, getReach2(N,encuesta$Reach,semilla) )
encuesta_pob[[1]]
encuesta_pob[[1]][encuesta_pob[[2]][[1]],3]

#Con esto ya se puede hacer recuento de poblaciones, cuantos conoce de cada una,
#contar en caso de pertenecer a la hidden population cuanta gente conoce que
#pertenece a la hidden aplicando el factor de corrección etc etc

#Con esto sería suficiente para implementar los estimadores que queríamos utilizar,
#el MoS, MLE y PIMLE. Falta discutir y aplicar los factores a cada uno de los pasos 
#en los cuales se puede producir un tipo de error.

