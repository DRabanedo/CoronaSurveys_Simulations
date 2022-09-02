library(dplyr)
library(matrixStats)
library(ggplot2)

simulation_data = read.csv("~/GitHub/CoronaSurveys_Simulations/R programs/Not disjoint subpopulations/Visibility factor/Analysis/Csv archives/seed_207/Simulation_visibilityfactorestimate_subpop_08_notdisjoint")

######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

simulation_data_v =  select(simulation_data, -starts_with("vf_subpop_") & -starts_with("X") & -starts_with("data"))
simulation_data_v_out = select(simulation_data, starts_with("vf_subpop_out"))

# We use -1 because we are comparing with the real value:  1
visibility_factor = 0.8
vf = visibility_factor

Visibility_factor_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(simulation_data_v-visibility_factor))),
                                        mse = rowMeans(as.matrix((simulation_data_v-visibility_factor)^2)),
                                        bias = rowMeans(as.matrix(simulation_data_v)),
                                        sd = rowSds(as.matrix(simulation_data_v)))

Visibility_factor_analysis_out = data.frame(abs_error = rowMeans(as.matrix(abs(simulation_data_v_out-visibility_factor))),
                                            mse = rowMeans(as.matrix((simulation_data_v_out-visibility_factor)^2)),
                                            bias = rowMeans(as.matrix(simulation_data_v_out)),
                                            sd = rowSds(as.matrix(simulation_data_v_out)))


################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### MSE #####

#Dataframe creation

graph_data_mse = data.frame(data = simulation_data$data, 
                            visibility_factor = Visibility_factor_analysis$mse, 
                            visibility_factor_out = Visibility_factor_analysis_out$mse)



ggplot(graph_data_mse) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate with remplacement")) + 
  geom_line(aes(x = data, y =  visibility_factor_out, col = "Visibility factor estimate without outliers")) + 
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the subpopulation memory factor",
       x = "Subpopulation memory factor",
       y = "Mean Squared Error (MSE)")


################################################################################

###### Bias analysis ######


graph_data_bias = data.frame(data = simulation_data$data, 
                             visibility_factor = Visibility_factor_analysis$bias,
                             visibility_factor_out = Visibility_factor_analysis_out$bias)



ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  rep(vf,length(simulation_data$data)), col = "Visibility factor")) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate with remplacement")) + 
  geom_line(aes(x = data, y =  visibility_factor_out, col = "Visibility factor estimate without outliers")) +
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the subpopulation memory factor",
       x = "Subpopulation memory factor",
       y = "Visibility factor value")



################################################################################

#### Standard deviation analysis ####

#Dataframe creation


graph_data_sd = data.frame(data = simulation_data$data, 
                           visibility_factor = Visibility_factor_analysis$sd,
                           visibility_factor_out = Visibility_factor_analysis_out$sd)



ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate with remplacement")) + 
  geom_line(aes(x = data, y =  visibility_factor_out, col = "Visibility factor estimate without outliers")) +
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the subpopulation memory factor",
       x = "Subpopulation memory factor",
       y = "Standard deviation")

#################################################################################

#### Absolute error ####


graph_data_abserror = data.frame(data = simulation_data$data, 
                                 visibility_factor = Visibility_factor_analysis$abs_error,
                                 visibility_factor_out = Visibility_factor_analysis_out$abs_error)



ggplot(graph_data_abserror) + 
  geom_line(aes(x = data, y =  visibility_factor, col = "Visibility factor estimate with remplacement")) + 
  geom_line(aes(x = data, y =  visibility_factor_out, col = "Visibility factor estimate without outliers")) + 
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the subpopulation memory factor",
       x = "Subpopulation memory factor",
       y = "Absolute error")

