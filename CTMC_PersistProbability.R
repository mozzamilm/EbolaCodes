  
  
###############################################################################
#       Stochastic simulation of the Seasonal Ebola model 

#                     Mozzamil Mohammed
#           Drake Lab, CEID, University of Georgia, USA
#                     April 2024
###############################################################################

# This code is to estimate the disease persistence probability in the human 
# settlement region. 
###############################################################################




results_seasonality <- list() # to save the model simulation results

extinction_probability <- vector() # to save the extinction probability

seasonality_vector <- seq(0.1,5,by = (5-0.1)/9) # 10 values of seasonality amplitude



for(season_indx in 1:10){ # looping over seasonality amplitude
  


num_stoch_real = 1 # initialize the number of stochastic realizations

while(num_stoch_real <= 20){ # looping over the number of stochastic realization

# initial population sizes in patches A and B
  
  N_A = 100
  N_B = 1000
  
  #############################################################################
  # initial values for the state variables
  
  S_A = 95
  I_A = N_A - S_A
  R_A = 0   
  S_B = 995
  I_B = N_B - S_B
  R_B = 0   
  

  #############################################################################
  # Model parameters 
  
  params <- list(
    N_A = 100,
    N_B = 1000,
    mu_ = 0.02,
    gamma_A = 12,
    gamma_B = 12,
    beta_A_base = 60, 
    beta_B_base = 5, 
    phi_AB_base = 0.4, 
    phi_BA_base = 0.04,
    seasonality = seasonality_vector[season_indx]
  )
  
  
  #############################################################################
  # saving simulation outputs and time
  
  results_ <- list(susceptible_A = c(S_A), infected_A = c(I_A), recovered_A = c(R_A),
                   susceptible_B = c(S_B), infected_B = c(I_B), recovered_B = c(R_B),
                   time_ = c(simulation_time))
  
  
  #############################################################################
  # stochastic sampling using while loop
  

  simulation_time = 0 # simulation time
    
  t_max = 10  # maximal simulation time in years
  
  while(simulation_time <= t_max){
    
    
    ###########################################################################
    # transmission rate with seasonality
    # change seasonality to zero for the system without seasonality
    
    beta_A = params$beta_A_base*(1+(params$seasonality*cos((2*pi*simulation_time))))
    
    beta_B = params$beta_B_base*(1+(params$seasonality*cos((2*pi*simulation_time))))

    ###########################################################################
    # Dispersal rate with seasonality
    #change seasonality to zero for the system without seasonality

    
    phi_AB = params$phi_AB_base*(1 + (params$seasonality*cos((2*pi*simulation_time))))  
    
    phi_BA = params$phi_BA_base*(1 + (params$seasonality*(cos((2*pi*simulation_time)))))  
    
    ###########################################################################
    # Model rates
    
    suceptible_birth_in_patch_A <- params$mu_*N_A
    susceptible_death_in_patch_A <- params$mu_*S_A
    infected_death_in_patch_A <- params$mu_*I_A
    recovered_death_in_patch_A <- params$mu_*R_A

    transmission_in_patch_A <- params$beta_A*(I_A/N_A)*S_A + 15 # we added a constant to avoid fluctuation-driven extinction
    recovery_in_patch_A <- params$gamma_A*I_A
    

    ###########################################################################
    
    suceptible_birth_in_patch_B <- params$mu_*N_B
    susceptible_death_in_patch_B <- params$mu_*S_B
    infected_death_in_patch_B <- params$mu_*I_B
    recovered_death_in_patch_B <- params$mu_*R_B
    
    transmission_in_patch_B <- params$beta_B*(I_B/N_B)*S_B
    recovery_in_patch_B <- params$gamma_B*I_B
    
    ###########################################################################
    
    susceptible_movement_AB <- params$phi_AB*S_A
    susceptible_movement_BA <- params$phi_BA*S_B
    infected_movement_AB <- params$phi_AB*I_A
    infected_movement_BA <- params$phi_BA*I_B
    recovered_movement_AB <- params$phi_AB*R_A
    recovered_movement_BA <- params$phi_BA*R_B
    
    rates_ <- c(suceptible_birth_in_patch_A, susceptible_death_in_patch_A, infected_death_in_patch_A,
                recovered_death_in_patch_A, transmission_in_patch_A, recovery_in_patch_A,   
                suceptible_birth_in_patch_B, susceptible_death_in_patch_B, infected_death_in_patch_B,
                recovered_death_in_patch_B, transmission_in_patch_B, recovery_in_patch_B,
                susceptible_movement_AB, susceptible_movement_BA,
                infected_movement_AB, infected_movement_BA,
                recovered_movement_AB, recovered_movement_BA)
    
    
    ###########################################################################
    # Sum of rates
    
    rates_sum = sum(rates_)
    
    ###########################################################################
    # Draw two uniformally distributed random numbers 
    
    R_1 = runif(1)
    R_2 = runif(1)
    
    ###########################################################################
    # randomly choose time to next event: drawn from an exponential distribution (Gillespie 1977)
    
    tau_ = -log(R_1)/rates_sum
  
  
    ###########################################################################
    # update simulation time
  
    simulation_time = simulation_time + tau_
    
    results_$time_ = append(results_$time_, simulation_time)
  
  
    ###########################################################################
    # randomly choose an event to occur: Gillespie algorithm condition (Gillespie 1977)
  
    check_condition = cumsum(rates_) > R_2*sum(rates_)
    which_satisfy = which(check_condition)

    
    ###########################################################################
    # which event to occur first: Gillespie algorithm condition
    event_to_occur = which_satisfy[1]

    
    ###########################################################################
    # updating state variables
  
   
   if (event_to_occur == 1) { # susceptible birth in patch A
     S_A = S_A + 1
   } else if (event_to_occur == 2) { # susceptible death in patch A
     S_A = S_A - 1
   } else if (event_to_occur == 3) { # Infected death in patch A
     I_A = I_A - 1
   } else if (event_to_occur == 4) { # recovered death in patch A
     R_A = R_A - 1
   } else if (event_to_occur == 5) { # transmission in patch A
     S_A = S_A - 1
     I_A = I_A + 1
   } else if (event_to_occur == 6) { # recovery in patch A
     I_A = I_A - 1
     R_A = R_A + 1
   }

    #############################################################################################
   
    
    if (event_to_occur == 7) { # susceptible birth in patch B
      S_B = S_B + 1
    } else if (event_to_occur == 8) { # susceptible death in patch B
      S_B = S_B - 1
    } else if (event_to_occur == 9) { # Infected death in patch B
      I_B = I_B - 1
    } else if (event_to_occur == 10) { # recovered death in patch B
      R_B = R_B - 1
    } else if (event_to_occur == 11) { # transmission in patch B
      S_B = S_B - 1
      I_B = I_B + 1
    } else if (event_to_occur == 12) { # recovery in patch B
      I_B = I_B - 1
      R_B = R_B + 1
    }
    
    ###################################################################################
   
   
   else if (event_to_occur == 13) { # susceptible movement from patch A to B
     S_A = S_A - 1
     S_B = S_B + 1
   } else if (event_to_occur == 14) { # susceptible movement from patch B to A
     S_A = S_A + 1
     S_B = S_B - 1
   } else if (event_to_occur == 15) { # infected movement from patch A to B
     I_A = I_A - 1
     I_B = I_B + 1
   } else if (event_to_occur == 16) { # infected movement from patch B to A
     I_A = I_A + 1
     I_B = I_B - 1
   } else if (event_to_occur == 17) { # recovered movement from patch A to B
     R_A = R_A - 1
     R_B = R_B + 1
   } else if (event_to_occur == 18) { # recovered movement from patch B to A
     R_A = R_A + 1
     R_B = R_B - 1
   }

    
    # update state variables
    
    results_$susceptible_A = append(results_$susceptible_A, S_A)
  
    results_$infected_A = append(results_$infected_A, I_A)
  
    results_$recovered_A = append(results_$recovered_A, R_A)
    
    results_$susceptible_B = append(results_$susceptible_B, S_B)
    
    results_$infected_B = append(results_$infected_B, I_B)
    
    results_$recovered_B = append(results_$recovered_B, R_B)

      
  } # ending while loop


  
  ##############################################################################################################
  # Disease extinction probability per unit time
  
  
  indx = as.integer(results_$infected_B > 0) # counting how many times does the disease persist in human settlement region
  
  extinction_probability[num_stoch_real] <- sum(indx)/length(indx) # calculating the persistence probability
  
  num_stoch_real <- num_stoch_real + 1 # advancing the number of stochastic realizations
  
  }
 
results_seasonality[[season_indx]] <- extinction_probability  # saving persistence probability as a function of seasonality amplitude
}

results_seasonality


################################################################################
# Plotting the outputs
################################################################################

par(mfrow = c(1,1), mai = c(0.9,0.9,0.4,0.2))
  boxplot(results_seasonality,
          xlab = "Seasonality strength",
          ylab = "Disease persistence probability",
          names = c(0.1,0.6,1.2,1.7,2.3,2.8,3.4,3.9,4.5,5),
          cex.lab = 1.2, cex.axis = 1.2,
          ylim = c(0,0.6),
          col = "orange",
          border = "brown"
  )
  
  title("Low-risk settlement patch", cex.main = 1.6, font.main = 1)