#####################################
## functions
####################################
funcs <- list(
  #######################
  # stochastic SEIR model simulation 
  #######################
  SEIR = function(states, par) {
    
    #######################
    # initialise and populate compartments
    #######################
    S <- E <- I<- R <- N <- vector("integer", length=par$n_patches)
    
    for (i in 1:par$n_patches) {
      S[i] <- states$N0[i] - states$I0[i]- states$E0[i]
      E[i] <- states$E0[i]
      I[i] <- states$I0[i]
      R[i] <- 0
      N[i] <- S[i]+E[i]+I[i]+R[i]
    }
    
    # recovery rate
    gamma <- 1/par$dur_infectious
    
    # proggression rate 
    sigma <- 1/par$dur_latent
    
    # mortality rate
    mu <- par$mort*gamma/(1-par$mort)
    
    #  contact matrix
    beta <- par$R0mat*(gamma+mu)
    
    # initialise storage matrices
    S_series <- E_series <- I_series <- R_series <- N_series <- matrix(0, nrow = par$T, ncol = par$n_patches)
    
    
    #######################
    # main time loop
    #######################
    for (t in 1:par$T) {
      S_series[t, ] <- S
      E_series[t, ] <- E
      I_series[t, ] <- I
      R_series[t, ] <- R
      N_series[t, ] <- S+I+R
      
      new_exposed <- new_endogenous_infection <-   new_exogenous_infection <- new_loss <- new_recovered <- new_dead <- rep(0, par$n_patches)
      
      for (patch in 1:par$n_patches) {
        
        
        #######################
        ## Event probabilities
        #######################
        
        if (N[patch]>0) {
          
          #exposure event
          p_exposure <- 1 - exp(-(
            (beta[patch, patch] * I[patch] / N[patch]) +
              sum(beta[patch, -patch] * I[-patch] / N[-patch], na.rm = T)))
          
          
          # any infection event
          p_infection <- 1 - exp(-(
            (beta[patch,patch] * I[patch] / N[patch]) +            # endogenous infection
              (sum(beta[patch,-patch]*I[-patch] / N[-patch], na.rm=T))      # exogenous infection
          ))
          
          # endogenous infection event
          p_endogenous_infection <- 1 - exp(-(
            (beta[patch,patch] * I[patch] / N[patch])
          ))
          
          # exogenous infection event
          p_exogenous_infection <- 1 - exp(-(
            (sum(beta[patch,-patch]*I[-patch] / N[-patch])) 
          ))
          
          # any loss event (recovery or death)
          
          p_loss <- 1 - exp(-(
            gamma +                                            # recovery 
              mu                                               # death
          ))
          
          # recovery event
          p_recovery <- 1 - exp(-(
            gamma
          ))
          
          # death event
          p_death <- 1 - exp(-(
            mu
          ))
          
          
        } else {
          
          p_exposure <- p_infection <- p_endogenous_infection <- p_exogenous_infection <-  p_loss <- p_recovery <-  p_death <- 0
          
        }
        
        #######################
        ## numbers of events
        #######################
        
        if (p_infection>0) {
          # new infection (endogenous or exogenous)
          new_infection[patch] <- rbinom(1, S[patch], p_infection)
          # endogenous infections
          new_endogenous_infection[patch] <- rbinom(1, new_infection[patch], p_endogenous_infection/p_infection)
          # exogenous infections
          new_exogenous_infection[patch] <- new_infection[patch] - new_endogenous_infection[patch]
        } else {
          new_infection[patch] <- new_endogenous_infection[patch] <- new_exogenous_infection[patch] <- 0
        }
        
        if (p_exposure > 0) {
          new_exposed[patch] <- rbinom(1, S[patch], p_exposure)
        } else {
          new_infection[patch] <- 0
        }
        
        if (p_loss>0) {
          # new loss (recovery or death)
          new_loss[patch] <- rbinom(1, I[patch], p_loss)
          # new recovery
          new_recovered[patch] <- rbinom(1, new_loss[patch], p_recovery/p_loss)
          # new death 
          new_dead[patch] <- new_loss[patch] - new_recovered[patch]
        } else {
          new_loss[patch] <- new_recovered[patch] <- new_dead[patch] <- 0
        }
        
      }
      
      #######################
      ## updates states
      #######################
      S <- S - new_exposed
      E <- E + new_exposed - new_infection 
      I <- I + new_endogenous_infection + new_exogenous_infection - new_recovered - new_dead
      R <- R + new_recovered
      N <- S + E + I + R
      
    }
    
    #######################
    ## outpur time series
    #######################
    
    list(S=S_series,E=E_series, I=I_series, R=R_series, N=N_series)
    
  }, 
  runSEIR = function(repeats, states, par) {
    
    multi <- replicate(repeats, funcs$SEIR(states=states, par=par))
    
    S <- as.data.frame(cbind(do.call(rbind, multi["S", 1:repeats]), time=seq(1, par$T)))
    E <- as.data.frame(cbind(do.call(rbind, multi["E", 1:repeats]), time=seq(1, par$T)))
    I <- as.data.frame(cbind(do.call(rbind, multi["I", 1:repeats]), time=seq(1, par$T)))
    R <- as.data.frame(cbind(do.call(rbind, multi["R", 1:repeats]), time=seq(1, par$T)))
    N <- as.data.frame(cbind(do.call(rbind, multi["N", 1:repeats]), time=seq(1, par$T)))
    
    S_df <- as.data.frame(S) %>% pivot_longer(-time, names_to="patch", values_to="S")
    E_df <- as.data.frame(E) %>% pivot_longer(-time, names_to="patch", values_to="E")
    I_df <- as.data.frame(I) %>% pivot_longer(-time, names_to="patch", values_to="I")
    R_df <- as.data.frame(R) %>% pivot_longer(-time, names_to="patch", values_to="R")
    N_df <- as.data.frame(N) %>% pivot_longer(-time, names_to="patch", values_to="N")
    
    S_out <- group_by(S_df, patch, time) %>% summarise(mean=mean(S), lwr=quantile(S, probs=0.025), upr=quantile(S, probs=0.975))
    E_out <- group_by(E_df, patch, time) %>% summarise(mean=mean(E), lwr=quantile(E, probs=0.025), upr=quantile(E, probs=0.975))
    I_out <- group_by(I_df, patch, time) %>% summarise(mean=mean(I), lwr=quantile(I, probs=0.025), upr=quantile(I, probs=0.975))
    R_out <- group_by(R_df, patch, time) %>% summarise(mean=mean(R), lwr=quantile(R, probs=0.025), upr=quantile(R, probs=0.975))
    N_out <- group_by(N_df, patch, time) %>% summarise(mean=mean(N), lwr=quantile(N, probs=0.025), upr=quantile(N, probs=0.975))
    
    return(list(S=S_out, E=E_out, I=I_out, R=R_out, N=N_out))
  }
)










