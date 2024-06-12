#####################################
## parameters
####################################
par <- list(n_patches=7,    
            T=100, 
            dur_infectious=13, ## duration of infectiousness (days)
            dur_latent=7,      ## duratino of latent period (days)
            mort=0.68,         ## proportion dying  
            R0intra=5)         ## intra-patch R0

par$R0mat <- matrix(c( par$R0intra, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                       0.5, par$R0intra, 0.5, 0.5, 0.5, 0.5, 0.5, 
                       0.5, 0.5, par$R0intra, 0.5, 0.5, 0.5, 0.5,
                       0.5, 0.5, 0.5, par$R0intra, 0.5, 0.5, 0.5, 
                       0.5, 0.5, 0.5, 0.5, par$R0intra, 0.5, 0.5,
                       0.5, 0.5, 0.5, 0.5, 0.5, par$R0intra, 0.5,
                       0.5, 0.5, 0.5, 0.5, 0.5, 0.5, par$R0intra),
                    nrow=par$n_patches, 
                    ncol=par$n_patches,byrow=T)

#####################################
## states
####################################
states = list(N0=  c(2, 6, 4, 1, 10, 5, 3), 
              I0 = c(0, 0, 0, 0, 1,  0, 0), 
              E0= c(0, 0, 0, 0, 1,  0, 0))