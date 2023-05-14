##########################################################
## -- Basic SEIR model
## -- COVID-19
## -- Marita Zimmermann
## -- May 2023
##########################################################

rm(list = ls())

# -- Load packages -- #
library(tidyverse)   
library(deSolve)

# -- Input parameters -- #
infectious_period <- 8 # days                  
latent_period     <- 4 # days                     
Rnaught           <- 3

# -- Initial population -- #
X <- 9990 # Susceptible     
E <- 0    # Exposed       
Y <- 10   # Infectious        
Z <- 0    # Recovered

# -- Compute model parameters -- #
gamma_value <- 1 / infectious_period
beta_value  <- Rnaught * gamma_value
sigma_value <- 1 / latent_period
N <- X + E + Y + Z

# -- Function to compute SEIR model -- #
seir_model <- function(current_timepoint, state_values, parameters) { 
  # create state variables (local variables)
  S <- state_values[1]        # susceptibles
  E <- state_values[2]        # exposed 
  I <- state_values[3]        # infectious
  R <- state_values[4]        # recovered
  
  with (as.list(parameters), {
      # compute derivatives
      dS <- -beta * S * I /N 
      dE <- beta * S * I / N - sigma * E
      dI <- sigma * E - gamma * I
      dR <- gamma * I
      # combine results
      results = c(dS, dE, dI, dR)
      list(results)
    } ) }

# -- Solve SEIR model -- #
output <- lsoda(y     = c(S = X, E = E, I = Y, R = Z),
                times = seq(0, 150, by=1),
                func  = seir_model, 
                parms = c(beta = beta_value, gamma = gamma_value, sigma = sigma_value))

# -- Visualize results -- #
output %>% as.data.frame() %>%
  gather(var, val, -time) %>%
  ggplot() +
  geom_line(aes(x = time, y = val, group = var, color = var), linewidth = 2, alpha = 0.8) +
  theme_bw(base_size = 11) +
  labs(x = "Days", y = "Population", color = NULL) +
  scale_color_viridis_d(option = "magma") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

