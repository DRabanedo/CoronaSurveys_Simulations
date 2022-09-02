library(dplyr)
library(matrixStats)
library(ggplot2)

simulation_data = read.csv("C:/Users/David Rabanedo/Desktop/Simulations_visibilityfactorestimate_vf")

######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.


simulation_data_reach =  select(simulation_data, -starts_with("vf_reach_out") & -starts_with("X") & -starts_with("data"))
simulation_data_reach_out = select(simulation_data, starts_with("vf_reach_out"))

simulation_data_subpop =  select(simulation_data, -starts_with("vf_subpop_") & -starts_with("X") & -starts_with("data"))
simulation_data_subpop_out = select(simulation_data, starts_with("vf_subpop_out"))


vf = simulation_data$data

vf_analysis_reach = data.frame(abs_error = rowMeans(as.matrix(abs(simulation_data_reach-vf))),
                                        mse = rowMeans(as.matrix((simulation_data_reach-vf)^2)),
                                        bias = rowMeans(as.matrix(simulation_data_reach)),
                                        sd = rowSds(as.matrix(simulation_data_reach)))

vf_analysis_reach_out = data.frame(abs_error = rowMeans(as.matrix(abs(simulation_data_reach_out-vf))),
                                            mse = rowMeans(as.matrix((simulation_data_reach_out-vf)^2)),
                                            bias = rowMeans(as.matrix(simulation_data_reach_out)),
                                            sd = rowSds(as.matrix(simulation_data_reach_out)))


vf_analysis_subpop = data.frame(abs_error = rowMeans(as.matrix(abs(simulation_data_subpop-vf))),
                                        mse = rowMeans(as.matrix((simulation_data_subpop-vf)^2)),
                                        bias = rowMeans(as.matrix(simulation_data_subpop)),
                                        sd = rowSds(as.matrix(simulation_data_subpop)))

vf_analysis_subpop_out = data.frame(abs_error = rowMeans(as.matrix(abs(simulation_data_subpop_out-vf))),
                                            mse = rowMeans(as.matrix((simulation_data_subpop_out-vf)^2)),
                                            bias = rowMeans(as.matrix(simulation_data_subpop_out)),
                                            sd = rowSds(as.matrix(simulation_data_subpop_out)))

################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### MSE #####

#Dataframe creation

graph_data_mse = data.frame(data = simulation_data$data, 
                            vf_reach = vf_analysis_reach$mse, 
                            vf_reach_out = vf_analysis_reach_out$mse,
                            vf_subpop = vf_analysis_subpop$mse, 
                            vf_subpop_out = vf_analysis_subpop_out$mse)



ggplot(graph_data_mse) + 
  geom_line(aes(x = data, y =  vf_reach, col = "Visibility factor estimate with remplacement (reach)")) + 
  geom_line(aes(x = data, y =  vf_reach_out, col = "Visibility factor estimate without outliers (reach)")) + 
  geom_line(aes(x = data, y =  vf_subpop, col = "Visibility factor estimate with remplacement (subpop)")) + 
  geom_line(aes(x = data, y =  vf_subpop_out, col = "Visibility factor estimate without outliers (subpop)")) + 
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the visibility factor",
       x = "Visibility factor",
       y = "Mean Squared Error (MSE)")


################################################################################

###### Bias analysis ######


graph_data_bias = data.frame(data = simulation_data$data, 
                             vf_reach = vf_analysis_reach$bias,
                             vf_reach_out = vf_analysis_reach_out$bias,
                             vf_subpop = vf_analysis_subpop$bias,
                             vf_subpop_out = vf_analysis_subpop_out$bias)



ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  vf_reach, col = "Visibility factor estimate with remplacement (reach)")) + 
  geom_line(aes(x = data, y =  vf_reach_out, col = "Visibility factor estimate without outliers (reach)")) + 
  geom_line(aes(x = data, y =  vf_subpop, col = "Visibility factor estimate with remplacement (subpop)")) + 
  geom_line(aes(x = data, y =  vf_subpop_out, col = "Visibility factor estimate without outliers (subpop)")) + 
  geom_line(aes(x = data, y =  vf, col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the visibility factor",
       x = "Visibility factor",
       y = "Visibility factor value")


################################################################################

#### Standard deviation analysis ####

#Dataframe creation

graph_data_sd = data.frame(data = simulation_data$data, 
                           vf_reach = vf_analysis_reach$sd,
                           vf_reach_out = vf_analysis_reach_out$sd,
                           vf_subpop = vf_analysis_subpop$sd,
                           vf_subpop_out = vf_analysis_subpop_out$sd)


ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  vf_reach, col = "Visibility factor estimate with remplacement (reach)")) + 
  geom_line(aes(x = data, y =  vf_reach_out, col = "Visibility factor estimate without outliers (reach)")) + 
  geom_line(aes(x = data, y =  vf_subpop, col = "Visibility factor estimate with remplacement (subpop)")) + 
  geom_line(aes(x = data, y =  vf_subpop_out, col = "Visibility factor estimate without outliers (subpop)")) +
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the visibility factor",
       x = "Visibility factor",
       y = "Standard deviation")

#################################################################################

#### Absolute error ####


graph_data_abserror = data.frame(data = simulation_data$data, 
                                 vf_reach = vf_analysis_reach$abs_error,
                                 vf_reach_out = vf_analysis_reach_out$abs_error,
                                 vf_subpop = vf_analysis_subpop$abs_error,
                                 vf_subpop_out = vf_analysis_subpop_out$abs_error)


ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  vf_reach, col = "Visibility factor estimate with remplacement (reach)")) + 
  geom_line(aes(x = data, y =  vf_reach_out, col = "Visibility factor estimate without outliers (reach)")) + 
  geom_line(aes(x = data, y =  vf_subpop, col = "Visibility factor estimate with remplacement (subpop)")) + 
  geom_line(aes(x = data, y =  vf_subpop_out, col = "Visibility factor estimate without outliers(subpop)")) + 
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the visibility factor",
       x = "Visibility factor",
       y = "Absolute error")

