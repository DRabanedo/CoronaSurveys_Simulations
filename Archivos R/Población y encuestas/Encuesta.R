#Número de personas sobre las que realizamos la encuesta
n_encuesta = 500
#Número de encuestas que realizamos
m_encuesta = 10


getEncuesta = function(n_enc, dataframe){
  #Se introduce un número n_enc que es el tamaño de la misma
  #Esta función genera aleatoriamente un dataframe con parte de las columnas del Pob_general
  encuesta = dataframe[sample(nrow(dataframe), n_enc, replace = FALSE),]
  encuesta
}


#Si queremos hacer m encuestas diferentes, simplemente usamos este bucle
lista_encuestas = list()
for(j in 1:m_encuesta){
  lista_encuestas[[j]] = getEncuesta(n_encuesta,Pob_general)
}
lista_encuestas
