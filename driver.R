################################################## 
## Stochastic SIR metapopulation model simulation
#################################################

## Clear workspace
rm(list=ls())

## Load required packages
source("packages.R")

## Load default model parameters
source("par.R")

## Load functions for running model
source("funcs.R")

## run a single 
out <- funcs$SEIR(states, par)

## Run the model with 100 repeats
out <- funcs$runSEIR(repeats=100, states=states, par=par)

## Plot the modeloutput
plot_susceptible <- ggplot(out$S, aes(x = time, y = mean, color = patch)) + 
  #geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=patch)) +
  geom_line() + 
  labs(y = "Susceptible")


plot_exposed <-  ggplot(out$I, aes(x = time, y = mean, color = patch)) + 
  # geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=patch)) +
  geom_line() + 
  labs(y = "Exposed")


plot_infected <-  ggplot(out$I, aes(x = time, y = mean, color = patch)) + 
  # geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=patch)) +
  geom_line() + 
  labs(y = "Infected")

plot_recovered <- ggplot(out$R, aes(x = time, y = mean, color = patch)) + 
  # geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=patch)) +
  geom_line() + 
  labs(y = "Recovered")

plot_pop <- ggplot(out$N, aes(x = time, y = mean, color = patch)) + 
  # geom_ribbon(aes(x=time, ymin=lwr, ymax=upr, fill=patch)) +
  geom_line() + 
  labs(y = "Population size") +
  scale_y_continuous(limits=c(0, max(states$N0)))

grid.arrange(plot_susceptible, plot_exposed, plot_infected,
             plot_recovered, plot_pop, ncol = 2)

