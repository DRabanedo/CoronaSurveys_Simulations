library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)

simulation_data = read.csv("C:/Users/David Rabanedo/Downloads/Simulation_visibilityfactorestimate_reach_notdisjoint_2022.csv")

seed_number = 2022

######### Data analysis ##########
# This way of presenting the data allows us to carry out a more detailed analysis
# of each estimator.

#Binomial aproach
vf_bin_data = select(simulation_data, starts_with("vf_bin"))

# Double truncation
vf_dt_bin_data = select(simulation_data, starts_with("vf_dt_bin"))
vf_dt_fbin_data = select(simulation_data, starts_with("vf_dt_fbin"))
vf_dt_dt_data = select(simulation_data, starts_with("vf_dt_dt"))
vf_dt_ut_data = select(simulation_data, starts_with("vf_dt_ut"))


# Unilateral truncation 
vf_ut_bin_data = select(simulation_data, starts_with("vf_ut_bin"))
vf_ut_fbin_data = select(simulation_data, starts_with("vf_ut_fbin"))
vf_ut_dt_data = select(simulation_data, starts_with("vf_ut_dt"))
vf_ut_ut_data = select(simulation_data, starts_with("vf_ut_ut"))

##################################

visibility_factor = 0.8
vf = visibility_factor



#Binomial aproach
vf_bin_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_bin_data-visibility_factor))),
                             mse = rowMeans(as.matrix((vf_bin_data-visibility_factor)^2)),
                             bias = rowMeans(as.matrix(vf_bin_data)),
                             sd = rowSds(as.matrix(vf_bin_data)))

# Double truncation
vf_dt_bin_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_dt_bin_data-visibility_factor))),
                                mse = rowMeans(as.matrix((vf_dt_bin_data-visibility_factor)^2)),
                                bias = rowMeans(as.matrix(vf_dt_bin_data)),
                                sd = rowSds(as.matrix(vf_dt_bin_data)))

vf_dt_fbin_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_dt_fbin_data-visibility_factor))),
                                 mse = rowMeans(as.matrix((vf_dt_fbin_data-visibility_factor)^2)),
                                 bias = rowMeans(as.matrix(vf_dt_fbin_data)),
                                 sd = rowSds(as.matrix(vf_dt_fbin_data)))

vf_dt_dt_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_dt_dt_data-visibility_factor))),
                               mse = rowMeans(as.matrix((vf_dt_dt_data-visibility_factor)^2)),
                               bias = rowMeans(as.matrix(vf_dt_dt_data)),
                               sd = rowSds(as.matrix(vf_dt_dt_data)))

vf_dt_ut_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_dt_ut_data-visibility_factor))),
                               mse = rowMeans(as.matrix((vf_dt_ut_data-visibility_factor)^2)),
                               bias = rowMeans(as.matrix(vf_dt_ut_data)),
                               sd = rowSds(as.matrix(vf_dt_ut_data)))


# Unilateral truncation 
vf_ut_bin_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_ut_bin_data-visibility_factor))),
                                mse = rowMeans(as.matrix((vf_ut_bin_data-visibility_factor)^2)),
                                bias = rowMeans(as.matrix(vf_ut_bin_data)),
                                sd = rowSds(as.matrix(vf_ut_bin_data)))

vf_ut_fbin_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_ut_fbin_data-visibility_factor))),
                                 mse = rowMeans(as.matrix((vf_ut_fbin_data-visibility_factor)^2)),
                                 bias = rowMeans(as.matrix(vf_ut_fbin_data)),
                                 sd = rowSds(as.matrix(vf_ut_fbin_data)))

vf_ut_dt_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_ut_dt_data-visibility_factor))),
                               mse = rowMeans(as.matrix((vf_ut_dt_data-visibility_factor)^2)),
                               bias = rowMeans(as.matrix(vf_ut_dt_data)),
                               sd = rowSds(as.matrix(vf_ut_dt_data)))

vf_ut_ut_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_ut_ut_data-visibility_factor))),
                               mse = rowMeans(as.matrix((vf_ut_ut_data-visibility_factor)^2)),
                               bias = rowMeans(as.matrix(vf_ut_ut_data)),
                               sd = rowSds(as.matrix(vf_ut_ut_data)))


################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### MSE #####

#Dataframe creation

plot_name = str_c("Simulation_visibilityfactor_estimate_", seed_number, "_mf_ ", memory_factor, "_notdisjoint_mse.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_mse = data.frame(data = simulation_data$data, 
                            
                            vf_bin     = vf_bin_analysis$mse,
                            
                            vf_dt_bin  = vf_dt_bin_analysis$mse,
                            vf_dt_fbin = vf_dt_fbin_analysis$mse,
                            vf_dt_dt   = vf_dt_dt_analysis$mse,
                            vf_dt_ut   = vf_dt_ut_analysis$mse,
                            
                            vf_ut_bin  = vf_ut_bin_analysis$mse, 
                            vf_ut_fbin = vf_ut_fbin_analysis$mse,
                            vf_ut_dt   = vf_ut_dt_analysis$mse,
                            vf_ut_ut   = vf_ut_ut_analysis$mse)



ggplot(graph_data_mse) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_ut, col = "vf_dt_ut")) + 
  
  geom_line(aes(x = data, y =  vf_ut_bin, col = "Vf_ut_bin")) + 
  geom_line(aes(x = data, y =  vf_ut_fbin, col = "Vf_ut_fbin")) +
  geom_line(aes(x = data, y =  vf_ut_dt, col = "Vf_ut_dt")) + 
  geom_line(aes(x = data, y =  vf_ut_ut, col = "Vf_ut_ut")) +
  
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Memory factor ", memory_factor), 
       x = "Memory factor",
       y = "Mean Squared Error (MSE)")

dev.off()

################################################################################

###### Bias analysis ######

plot_name = str_c("Simulation_visibilityfactor_estimate_" , seed_number, "_mf_ ", memory_factor, "_notdisjoint_bias.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_bias = data.frame(data = simulation_data$data, 
                             
                             vf_bin     = vf_bin_analysis$bias,
                             
                             vf_dt_bin  = vf_dt_bin_analysis$bias,
                             vf_dt_fbin = vf_dt_fbin_analysis$bias,
                             vf_dt_dt   = vf_dt_dt_analysis$bias,
                             vf_dt_ut   = vf_dt_ut_analysis$bias,
                             
                             vf_ut_bin  = vf_ut_bin_analysis$bias, 
                             vf_ut_fbin = vf_ut_fbin_analysis$bias,
                             vf_ut_dt   = vf_ut_dt_analysis$bias,
                             vf_ut_ut   = vf_ut_ut_analysis$bias)



ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_ut, col = "vf_dt_ut")) + 
  
  geom_line(aes(x = data, y =  vf_ut_bin, col = "Vf_ut_bin")) + 
  geom_line(aes(x = data, y =  vf_ut_fbin, col = "Vf_ut_fbin")) +
  geom_line(aes(x = data, y =  vf_ut_dt, col = "Vf_ut_dt")) + 
  geom_line(aes(x = data, y =  vf_ut_ut, col = "Vf_ut_ut")) +
  
  geom_line(aes(x = data, y =  rep(vf,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Memory factor ", memory_factor),
       x = "Memory factor",
       y = "Visibility factor value")

dev.off()

################################################################################

#### Standard deviation analysis ####

#Dataframe creation

plot_name = str_c("Simulation_visibilityfactor_estimate_", seed_number,"_mf_ ", memory_factor, "_notdisjoint_sd.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_sd = data.frame(data = simulation_data$data, 
                           
                           vf_bin     = vf_bin_analysis$sd,
                           
                           vf_dt_bin  = vf_dt_bin_analysis$sd,
                           vf_dt_fbin = vf_dt_fbin_analysis$sd,
                           vf_dt_dt   = vf_dt_dt_analysis$sd,
                           vf_dt_ut   = vf_dt_ut_analysis$sd,
                           
                           vf_ut_bin  = vf_ut_bin_analysis$sd, 
                           vf_ut_fbin = vf_ut_fbin_analysis$sd,
                           vf_ut_dt   = vf_ut_dt_analysis$sd,
                           vf_ut_ut   = vf_ut_ut_analysis$sd)



ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_ut, col = "vf_dt_ut")) + 
  
  geom_line(aes(x = data, y =  vf_ut_bin, col = "Vf_ut_bin")) + 
  geom_line(aes(x = data, y =  vf_ut_fbin, col = "Vf_ut_fbin")) +
  geom_line(aes(x = data, y =  vf_ut_dt, col = "Vf_ut_dt")) + 
  geom_line(aes(x = data, y =  vf_ut_ut, col = "Vf_ut_ut")) +
  
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Memory factor ", memory_factor),
       x = "Memory factor",
       y = "Standard deviation")

dev.off()

#################################################################################

#### Absolute error ####

plot_name = str_c("Simulation_visibilityfactor_estimate_", seed_number, "_mf_ ", memory_factor, "_notdisjoint_abserror.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_abs_error = data.frame(data = simulation_data$data, 
                                  
                                  vf_bin     = vf_bin_analysis$abs_error,
                                  
                                  vf_dt_bin  = vf_dt_bin_analysis$abs_error,
                                  vf_dt_fbin = vf_dt_fbin_analysis$abs_error,
                                  vf_dt_dt   = vf_dt_dt_analysis$abs_error,
                                  vf_dt_ut   = vf_dt_ut_analysis$abs_error,
                                  
                                  vf_ut_bin  = vf_ut_bin_analysis$abs_error, 
                                  vf_ut_fbin = vf_ut_fbin_analysis$abs_error,
                                  vf_ut_dt   = vf_ut_dt_analysis$abs_error,
                                  vf_ut_ut   = vf_ut_ut_analysis$abs_error)



ggplot(graph_data_abs_error) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_ut, col = "vf_dt_ut")) + 
  
  geom_line(aes(x = data, y =  vf_ut_bin, col = "Vf_ut_bin")) + 
  geom_line(aes(x = data, y =  vf_ut_fbin, col = "Vf_ut_fbin")) +
  geom_line(aes(x = data, y =  vf_ut_dt, col = "Vf_ut_dt")) + 
  geom_line(aes(x = data, y =  vf_ut_ut, col = "Vf_ut_ut")) +
  
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Memory factor ", memory_factor),
       x = "Memory factor",
       y = "Absolute error")
dev.off()

