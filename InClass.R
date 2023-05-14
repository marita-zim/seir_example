# In class question

# -- Add vaccination to the model if vaccination is occurring over time -- #
# assume people from any state can get vaccinated
# for this question, what would happen if people were vaccinated at a rate of 1% per day?
p = 0.01 # vaccination rate

# Function to compute SEIR model
seir_model <- function(current_timepoint, state_values, parameters) { 
  # create state variables (local variables)
  S <- state_values[1]        # susceptibles
  E <- state_values[2]        # exposed 
  I <- state_values[3]        # infectious
  R <- state_values[4]        # recovered
  
  with (as.list(parameters), {
    # compute derivatives
    dS <- -beta * S * I /N - p * S
    dE <- beta * S * I / N - sigma * E - p * E
    dI <- sigma * E - gamma * I - p * I
    dR <- gamma * I + p * (S+E+I)
    # combine results
    results = c(dS, dE, dI, dR)
    list(results)
  } ) }

# Solve model 
output_vax <- lsoda(y     = c(S = X, E = E, I = Y, R = Z),
                    times = seq(0, 500, by=1),
                    func  = seir_model, 
                    parms = c(beta = beta_value, gamma = gamma_value, sigma = sigma_value))

# Visualize results 
output_vax %>% as.data.frame() %>%
  gather(var, val, -time) %>%
  ggplot() +
  geom_line(aes(x = time, y = val, group = var, color = var), size = 2, alpha = 0.8) +
  theme_bw(base_size = 11) +
  labs(x = "Days", y = "Population", color = NULL) +
  scale_color_viridis_d(option = "magma") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
