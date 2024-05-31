

###############################################################################
# the following function is for the two-patch deterministic model presented 
# in the main text of the manuscript
# it will be called in the main script called Deterministic_Stochastic_Model.Rmd
###############################################################################


Two_Patch_Ebola_Model <- function(time, state_var, params) {
  with(as.list(c(state_var, params)), {
    
    
    ###########################################################################
    # transmission rate with seasonality
    # change seasonality to zero for the system without seasonality
    
    beta_A <- function(time){
      
      return(beta_A_base*(1+(seasonality*cos((2*pi*time)))))
    }
    
    beta_B <- function(time){
      
      return(beta_B_base*(1+(seasonality*cos((2*pi*time)))))
    }
    
    ###########################################################################
    # Dispersal rate with seasonality
    #change seasonality to zero for the system without seasonality
    
    
    phi_AB <- function(time){
      
      return(phi_AB_base*(1 + (seasonality*cos((2*pi*time)))))
    }
    
    phi_BA <- function(time){
      
      return(phi_BA_base*(1 + (seasonality*(cos((2*pi*time))))))
    }
    
    ###########################################################################
    
    # reservoir patch (patch A)
    
    dS_A <- (mu_ * N_A) - ((beta_A(time))*(I_A/N_A)*(S_A)) -
      ((mu_ + phi_AB(time))*S_A) + (phi_BA(time)*S_B) 
    
    dI_A <- (beta_A(time)*(I_A/N_A)*S_A) - ((mu_ + gamma_A + 
      phi_AB(time)) * I_A) + (phi_BA(time)*I_B) 
    
    dR_A <- (gamma_A*I_A)-((mu_ + phi_AB(time))*R_A) +
      (phi_BA(time)*R_B)  
    
    ###########################################################################
    
    # human settlement patch (patch B)
    
    dS_B <-(mu_ * N_B) - ((beta_B(time))*(I_B/N_B)*(S_B)) - 
      ((mu_ + phi_BA(time))*S_B) + (phi_AB(time)*S_A) 
    
    dI_B <-(beta_B(time)*(I_B/N_B)*S_B) - ((mu_ + gamma_B +
        phi_BA(time))*I_B) + ((time)*I_A) 
    
    dR_B <-(gamma_B*I_B)-((mu_ + phi_BA(time))*R_B) +
      (phi_AB(time)*R_A) 
    
    return(list(c(dS_A, dI_A, dR_A, dS_B, dI_B, dR_B)))
  })
}