library(dplyr)
library(matrixStats)
library(ggplot2)
library(stringr)

simulation_data = read.csv("~/GitHub/CoronaSurveys_Simulations/R programs/Disjoint and no disjoint ensemble/Visibility factor/CSV/seed_2022/Simulation_visibilityfactorestimate_0.6_reach_notdisjoint_utc_2022.csv")
  
seed_number = 2022
visibility_factor = 0.6
vf = visibility_factor

######### Data analysis ##########
# This way of presenting the data allows us to carry outc a more detailed analysis
# of each estimator.

#Binomial aproach
vf_bin_data = select(simulation_data, starts_with("vf_bin"))

# Double truncation
vf_dt_bin_data = select(simulation_data, starts_with("vf_dt_bin"))
vf_dt_fbin_data = select(simulation_data, starts_with("vf_dt_fbin"))
vf_dt_dt_data = select(simulation_data, starts_with("vf_dt_dt"))
vf_dt_utc_data = select(simulation_data, starts_with("vf_dt_utc"))


# Unilateral truncation 
vf_utc_bin_data = select(simulation_data, starts_with("vf_utc_bin"))
vf_utc_fbin_data = select(simulation_data, starts_with("vf_utc_fbin"))
vf_utc_dt_data = select(simulation_data, starts_with("vf_utc_dt"))
vf_utc_utc_data = select(simulation_data, starts_with("vf_utc_utc"))

##################################





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

vf_dt_utc_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_dt_utc_data-visibility_factor))),
                               mse = rowMeans(as.matrix((vf_dt_utc_data-visibility_factor)^2)),
                               bias = rowMeans(as.matrix(vf_dt_utc_data)),
                               sd = rowSds(as.matrix(vf_dt_utc_data)))


# Unilateral truncation 
vf_utc_bin_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_utc_bin_data-visibility_factor))),
                                mse = rowMeans(as.matrix((vf_utc_bin_data-visibility_factor)^2)),
                                bias = rowMeans(as.matrix(vf_utc_bin_data)),
                                sd = rowSds(as.matrix(vf_utc_bin_data)))

vf_utc_fbin_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_utc_fbin_data-visibility_factor))),
                                 mse = rowMeans(as.matrix((vf_utc_fbin_data-visibility_factor)^2)),
                                 bias = rowMeans(as.matrix(vf_utc_fbin_data)),
                                 sd = rowSds(as.matrix(vf_utc_fbin_data)))

vf_utc_dt_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_utc_dt_data-visibility_factor))),
                               mse = rowMeans(as.matrix((vf_utc_dt_data-visibility_factor)^2)),
                               bias = rowMeans(as.matrix(vf_utc_dt_data)),
                               sd = rowSds(as.matrix(vf_utc_dt_data)))

vf_utc_utc_analysis = data.frame(abs_error = rowMeans(as.matrix(abs(vf_utc_utc_data-visibility_factor))),
                               mse = rowMeans(as.matrix((vf_utc_utc_data-visibility_factor)^2)),
                               bias = rowMeans(as.matrix(vf_utc_utc_data)),
                               sd = rowSds(as.matrix(vf_utc_utc_data)))


################################################################################

####### Graph representation #######
# The graphic representation allows us to compare the different results obtained
# by each one of the estimators

##### MSE #####

#Dataframe creation

plot_name = str_c("Simulation_visibilityfactor_estimate_", visibility_factor, " _", seed_number, "_notdisjoint_utc_mse.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_mse = data.frame(data = simulation_data$data, 
                            
                            vf_bin     = vf_bin_analysis$mse,
                            
                            vf_dt_bin  = vf_dt_bin_analysis$mse,
                            vf_dt_fbin = vf_dt_fbin_analysis$mse,
                            vf_dt_dt   = vf_dt_dt_analysis$mse,
                            vf_dt_utc   = vf_dt_utc_analysis$mse,
                            
                            vf_utc_bin  = vf_utc_bin_analysis$mse, 
                            vf_utc_fbin = vf_utc_fbin_analysis$mse,
                            vf_utc_dt   = vf_utc_dt_analysis$mse,
                            vf_utc_utc   = vf_utc_utc_analysis$mse)



ggplot(graph_data_mse) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_utc, col = "vf_dt_utc")) + 
  
  geom_line(aes(x = data, y =  vf_utc_bin, col = "Vf_utc_bin")) + 
  geom_line(aes(x = data, y =  vf_utc_fbin, col = "Vf_utc_fbin")) +
  geom_line(aes(x = data, y =  vf_utc_dt, col = "Vf_utc_dt")) + 
  geom_line(aes(x = data, y =  vf_utc_utc, col = "Vf_utc_utc")) +
  
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Visibility factor ", visibility_factor), 
       x = "Memory factor",
       y = "Mean Squared Error (MSE)")

dev.off()

################################################################################

###### Bias analysis ######

plot_name = str_c("Simulation_visibilityfactor_estimate_", visibility_factor, " _", seed_number, "_notdisjoint_utc_bias.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_bias = data.frame(data = simulation_data$data, 
                            
                            vf_bin     = vf_bin_analysis$bias,
                            
                            vf_dt_bin  = vf_dt_bin_analysis$bias,
                            vf_dt_fbin = vf_dt_fbin_analysis$bias,
                            vf_dt_dt   = vf_dt_dt_analysis$bias,
                            vf_dt_utc   = vf_dt_utc_analysis$bias,
                            
                            vf_utc_bin  = vf_utc_bin_analysis$bias, 
                            vf_utc_fbin = vf_utc_fbin_analysis$bias,
                            vf_utc_dt   = vf_utc_dt_analysis$bias,
                            vf_utc_utc   = vf_utc_utc_analysis$bias)



ggplot(graph_data_bias) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_utc, col = "vf_dt_utc")) + 
  
  geom_line(aes(x = data, y =  vf_utc_bin, col = "Vf_utc_bin")) + 
  geom_line(aes(x = data, y =  vf_utc_fbin, col = "Vf_utc_fbin")) +
  geom_line(aes(x = data, y =  vf_utc_dt, col = "Vf_utc_dt")) + 
  geom_line(aes(x = data, y =  vf_utc_utc, col = "Vf_utc_utc")) +
  
  geom_line(aes(x = data, y =  rep(vf,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Visibility factor ", visibility_factor),
       x = "Memory factor",
       y = "Visibility factor value")

dev.off()

################################################################################

#### Standard deviation analysis ####

#Dataframe creation

plot_name = str_c("Simulation_visibilityfactor_estimate_", visibility_factor, " _", seed_number, "_notdisjoint_utc_sd.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_sd = data.frame(data = simulation_data$data, 
                            
                            vf_bin     = vf_bin_analysis$sd,
                            
                            vf_dt_bin  = vf_dt_bin_analysis$sd,
                            vf_dt_fbin = vf_dt_fbin_analysis$sd,
                            vf_dt_dt   = vf_dt_dt_analysis$sd,
                            vf_dt_utc   = vf_dt_utc_analysis$sd,
                            
                            vf_utc_bin  = vf_utc_bin_analysis$sd, 
                            vf_utc_fbin = vf_utc_fbin_analysis$sd,
                            vf_utc_dt   = vf_utc_dt_analysis$sd,
                            vf_utc_utc   = vf_utc_utc_analysis$sd)



ggplot(graph_data_sd) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_utc, col = "vf_dt_utc")) + 
  
  geom_line(aes(x = data, y =  vf_utc_bin, col = "Vf_utc_bin")) + 
  geom_line(aes(x = data, y =  vf_utc_fbin, col = "Vf_utc_fbin")) +
  geom_line(aes(x = data, y =  vf_utc_dt, col = "Vf_utc_dt")) + 
  geom_line(aes(x = data, y =  vf_utc_utc, col = "Vf_utc_utc")) +
  
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Visibility factor ", visibility_factor),
       x = "Memory factor",
       y = "Standard deviation")

dev.off()

#################################################################################

#### Absolutce error ####

plot_name = str_c("Simulation_visibilityfactor_estimate_", visibility_factor, " _", seed_number, "_notdisjoint_utc_abserror.png")


png(filename = plot_name,
    width = 1000, height = 600)

graph_data_abs_error = data.frame(data = simulation_data$data, 
                            
                            vf_bin     = vf_bin_analysis$abs_error,
                            
                            vf_dt_bin  = vf_dt_bin_analysis$abs_error,
                            vf_dt_fbin = vf_dt_fbin_analysis$abs_error,
                            vf_dt_dt   = vf_dt_dt_analysis$abs_error,
                            vf_dt_utc   = vf_dt_utc_analysis$abs_error,
                            
                            vf_utc_bin  = vf_utc_bin_analysis$abs_error, 
                            vf_utc_fbin = vf_utc_fbin_analysis$abs_error,
                            vf_utc_dt   = vf_utc_dt_analysis$abs_error,
                            vf_utc_utc   = vf_utc_utc_analysis$abs_error)



ggplot(graph_data_abs_error) + 
  geom_line(aes(x = data, y =  vf_bin, col = "Vf_bin")) + 
  
  geom_line(aes(x = data, y =  vf_dt_bin, col = "Vf_dt_bin")) +
  geom_line(aes(x = data, y =  vf_dt_fbin, col = "Vf_dt_fbin")) + 
  geom_line(aes(x = data, y =  vf_dt_dt, col = "Vf_dt_dt")) +
  geom_line(aes(x = data, y =  vf_dt_utc, col = "vf_dt_utc")) + 
  
  geom_line(aes(x = data, y =  vf_utc_bin, col = "Vf_utc_bin")) + 
  geom_line(aes(x = data, y =  vf_utc_fbin, col = "Vf_utc_fbin")) +
  geom_line(aes(x = data, y =  vf_utc_dt, col = "Vf_utc_dt")) + 
  geom_line(aes(x = data, y =  vf_utc_utc, col = "Vf_utc_utc")) +
  
  #geom_line(aes(x = data, y =  rep(1,length(simulation_data$data)), col = "Visibility factor")) + 
  scale_color_discrete("Legend") + 
  labs(title = "Visibility factor estimates based on the memory factor",
       subtitle = str_c("Seed ", seed_number, ", Visibility factor ", visibility_factor),
       x = "Memory factor",
       y = "Absolutce error")
dev.off()

