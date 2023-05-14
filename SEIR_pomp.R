##########################################################
## -- SEIR model using pomp package
## -- COVID-19
## -- Marita Zimmermann
## -- May 2021
##########################################################

rm(list = ls())

# -- Load packages -- #
library(tidyverse)
library(pomp)
library(readxl)
library(scales)

# -- Load data -- #
# source: https://www.doh.wa.gov/Emergencies/COVID19/DataDashboard#downloads
tmp = tempfile(fileext = ".xlsx")
download.file(url = "https://www.doh.wa.gov/Portals/1/Documents/1600/coronavirus/data-tables/WA_COVID19_Cases_Hospitalizations_Deaths.xlsx", destfile = tmp, mode = "wb")
case_data <- read_excel(tmp, sheet = "Cases") %>%
  group_by(WeekStartDate) %>% summarise(reports = sum(ProbableCases, na.rm = T)) %>% # format case data. Only need date and number of cases.
  arrange(WeekStartDate) %>% mutate(week = as.numeric(row_number())) %>% select(-WeekStartDate) %>%
  filter(between(week,35,60)) %>% mutate(week = week - min(week)) # Select one epidemiologic peak (this is the the biggest peak over the winter) since it's hard to model the whole epidemic with a simple model (any ideas why?)



# -- Simulation -- #

# stochastic simulator for the unobserved state process
seir_step <- function(S, E, I, R, N, H, Beta, sigma, gamma, delta.t, ...) {
  dN_SE <- rbinom(n=1, size=S, prob = 1 - exp(-Beta * I/N * delta.t)) # S -> E = Beta * I /N
  dN_EI <- rbinom(n=1, size=E, prob = 1 - exp(-sigma * delta.t))      # E -> I = sigma
  dN_IR <- rbinom(n=1, size=I, prob = 1 - exp(-gamma * delta.t))      # I -> R = gamma
  S <- S - dN_SE
  E <- E + dN_SE - dN_EI
  I <- I + dN_EI - dN_IR
  R <- R + dN_IR;
  H <- H + dN_IR; # need to add another variable to track true (not reported) incidence, call it H, assume related to I -> R recovery
  c(S = S, E = E, I = I, R = R, H = H) }

# At day zero population, assume 1 infection, but don't know how many are susceptible, assume fraction eta
seir_init <- function(N, eta, ...) {
  c(S = round(N*eta), E = 0, I = 1, R = round(N*(1-eta)), H = 0) }

# Model the data as a binomial process
# need to add probability of reporting, rho
# Use dmeasure or an rmeasure component, or both:
dmeas <- function(reports, H, rho, log, ...) {
  dbinom(x=reports, size = H, prob=rho, log=log) }
rmeas <- function(H, rho, ...) {
  c(reports=rbinom(n=1, size = H, prob=rho)) }

# Fold these basic model components, with the data, into a pomp object
# Want H to tally only the incidence over the week, reset it to zero at the beginning of each week. Use accumvars argument
covidSEIR <- case_data %>%
  pomp(times      = "week",
       t0         = 0,
       rmeasure   = rmeas,
       dmeasure   = dmeas,
       accumvars  = "H",
       rprocess   = pomp::euler(seir_step, delta.t=1/7),
       rinit      = seir_init,
       paramnames = c("N","Beta","sigma","gamma","rho","eta"),
       statenames = c("S","E","I","R","H")) 

# make some assumptions about parameters to try the simulation and see if it works
sims <- covidSEIR %>%
  simulate(params = c(Beta  = 3.5,    
                      sigma = 13.1,   
                      gamma = 2.9,   
                      rho   = 0.25,     # probability of reporting
                      eta   = 1,        # initial susceptible proportion
                      N     = 50000),   # population 
           nsim = 20, format = "data.frame", include.data = TRUE) 

# Visualize results
sims %>%
  ggplot()+
  geom_line(aes(x=week, y=reports, group=.id, color=.id=="data", alpha=.id=="data")) +
  theme_bw(base_size = 11) + guides(alpha = "none") +
  scale_color_manual(NULL, labels = c("Model", "Data"), values = c("blue","black")) +
  scale_y_continuous("Cases reported", labels = comma) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





# -- Use C Snippets -- #
# A C snippet is a small piece of C code used to specify a basic model component.
# computational speed that we can obtain only by using compiled native code

# Need to install RTools if don't already have it
# https://cran.r-project.org/bin/windows/Rtools/
# If you're using a version between 3.3 and 3.6 (which most of us probably are), you should use this installer: https://cran.r-project.org/bin/windows/Rtools/Rtools35.exe

# re-run your pomp object above to replace four pieces of code above with the equivalent snippets below

seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-sigma*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")
  
  
seir_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")
  
dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
")
  
rmeas <- Csnippet("
  reports = rbinom(H,rho);
")
  


# -- parameter optimization by SSE -- #
# more on optimization https://kingaa.github.io/clim-dis/parest/parest.html#feature-based-parameter-estimation

# create a list of all the parameter combinations you want to try
grid <- expand.grid(Beta  = seq(from=1.5,to=4.5,length=7),
                    sigma = seq(from=0.1,to=15,length=17),
                    gamma = seq(from=0.1,to=3,length=15),
                    rho   = 0.25,
                    eta   = 1,
                    N     = 50000)

# simulate all of your results with the pomp object you made above (note that you'll want more than 5 sims each)
sim.results <- simulate(covidSEIR, params=t(grid), format="data.frame", nsim = 5, include.data = F)

# aggregate results and calculate SSE
sse.results <- sim.results %>%
  mutate(.id = as.character(.id)) %>% separate(".id", into = c("id", "sim")) %>%
  left_join(as.data.frame(covidSEIR) %>% rename(reports.data = reports), by="week") %>%
  group_by(id) %>% summarise(sse=sum((reports-reports.data)^2)) %>%
  bind_cols(grid)
View(sse.results) 
