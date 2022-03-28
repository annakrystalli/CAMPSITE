
cs_sim_next_trait_value  <- function(traits, active_lineages, trait_diff, ou, 
                                     pars, bounds, step_size) {
  
  trait_diff_sgn <- cs_trait_diff_sgn(trait_diff)
  n_lineages <- length(active_lineages)
  
  for (i in 1:n_lineages) {
    lin_id <- active_lineages[i]
    lin_i <- which(colnames(trait_diff) == lin_id)
    sign_t <- trait_diff_sgn[lin_i, -lin_i] 
    diff_t <- trait_diff[lin_i, -lin_i]
    
    trait_values <- traits[[lin_id]]$trait_values
    current_trait_value <- tail(trait_values, 1)
    
    next_trait_value <- evolve_trait(x = current_trait_value, alpha2 = pars$alpha2, 
                                     m = pars$m, s = pars$s, opt = ou$opt, sign_t, diff_t, bounds)
    
    if(!is.null(ou$opt)){
      next_trait_value <- add_ou_effect(x1 = next_trait_value, x0 = current_trait_value, opt= ou$opt, 
                                        alpha4 = ou$alpha4, step_size = step_size)
    }
    traits[[lin_id]]$trait_values <- c(trait_values, next_trait_value)
  }
  traits
}


cs_evolve_inc_lineage <- function(pars, lin_id, step_size, t, e) {
  trait <- e$traits[[lin_id]]
  stopifnot(trait$status == -1)
  
  #cur_trait_value <- current_trait_value(trait)
  
  parent_trait <- e$traits[[trait$relative_id]]
  trait$parent_trait_value <- current_trait_value(parent_trait) # update current trait of parental lineage
  
  lambda2 <-  calc_lamda2(trait, pars$tau0, pars$beta)
  
  mui <- calc_mui(trait, trait_diff = trait_diff_extr(e$trait_diff_m, lin_id), 
                  alpha1 = pars$alpha1, mui0 = pars$mui0, alpha3 = pars$alpha3, 
                  mui1 = pars$mui1, muibg = pars$muibg)
  
  probs_i <- calc_probs(mui, lambda2)
  
  if (cease_state(lambda2, mui, step_size)) { 
    event <- sample(1:2, size = 1, prob = probs_i)
    if (event == 1) { #speciation completed
      # update lineage status to good in lineage matrix, update speciation completion times
      e$lineages[lin_id, c("status", "spec_or_ext_ct", "spec_ct")] <- c(1, t, t) 
      # update lineage status to good in trait list
      trait$status <- 1
     
    }
    if (event == 2) { #extinction
      # update lineage status to extinct in lineage matrix, update ending time and extinction time
      e$lineages[lin_id, c("end_time", "status", "spec_or_ext_ct")] <- c(t, -2, t) 
      # update lineage status to extinct in trait list
      trait$status <- - 2
    }
  }
  e$traits[[lin_id]] <- trait
  e$dead_lin <- lineages_dead(e$lineages)
  return(e)
 
}



cs_evolve_good_lineage <- function(pars, lin_id, step_size, t, ou, e) {
  
  trait <- e$traits[[lin_id]]
  stopifnot(trait$status == 1)
  
  trait_diff_opt <- calc_diff_trait_opt(trait, ou$opt)
  trait_diff <- trait_diff_extr(e$trait_diff_m, lin_id)
  
  mu <- calc_mu(trait_diff, pars$alpha1, pars$alpha3, pars$mu0, pars$mu1, pars$mubg, 
                trait_diff_opt)
  
  probs <- calc_probs(mu, pars$lambda1)
  
  if (cease_state(pars$lambda1, mu, step_size)) { #probability that lineage does not remain good
    event <- sample(1:2, size = 1, prob = probs)
    if (event == 1) { #speciation initiated
      #the current lineage is ended, and is saved as the ancestral node - in its stead 
      # are created two new lineages, the new parental lineage (same lineage, new 
      # number) and the new incipient lineage:
      e$lineages[lin_id, "descendant_node"] <- lineage_next_parent(e$lineages)  # add descendant node number
      e$lineages[lin_id, "end_time"] <- t #update ending time
      e$lineages <- lineage_add(e$lineages, 
                                #  parent lineage
                                new_lineage(c(e$lineages[lin_id, "descendant_node"], 
                                              0, t, -1, 1, t, t)),
                                # incipient lineage
                                new_lineage(c(e$lineages[lin_id, "descendant_node"], 
                                              0, t, -1, -1, NA, NA)))
  
      #update trait list to include new parental lineages
      e$traits[[lineage_last_parent(e$lineages)]] <- new_trait(lin_id = lin_id_to_int(lineage_last_parent(e$lineages)), 
                                                               status = 1, 
                                                               # because this is parental, store number of daughter lineage 
                                                               relative_id = lin_id_to_int(lineage_last_incipient(e$lineages)), 
                                                               parent_trait_value = NA, #trait value for parent lineage
                                                               #trait values for all previous time steps (NA since it didn't exist)
                                                               #trait value for next time step (trait of parent lineage)
                                                               trait_values = c(rep(NA, times = e$step_count),  current_trait_value(trait)))
      
      # update trait list to include new incipient lineages
      e$traits[[lineage_last_incipient(e$lineages)]] <- new_trait(lin_id = lin_id_to_int(lineage_last_incipient(e$lineages)), 
                                                                  status = -1, 
                                                                  # because this is parental, store number of parent lineage 
                                                                  relative_id = lin_id_to_int(lineage_last_parent(e$lineages)), 
                                                                  parent_trait_value = current_trait_value(trait), #trait value for parent lineage
                                                                  #trait values for all previous time steps (NA since it didn't exist)
                                                                  #trait value for next time step (trait of parent lineage)
                                                                  trait_values = c(rep(NA, times = e$step_count),  
                                                                                   current_trait_value(trait)))
                                                    
      # extract lineages who have i as parent
      daughters <- traits_find_daughters(e$traits, lin_id, format = "lin_id") 
      daughters <- daughters[e$lineages[daughters, "descendant_node"] == 0] #keep only living (drop all non-tips)
      if (length(daughters) > 0) {
        for (daughter in daughters) {
          # update parental lineage number to its new number
          e$traits[[daughter]]$relative_id <- lin_id_to_int(lineage_last_parent(e$lineages))
        }
      }
    }
    if (event == 2) { #extinction
      # update lineage status to extinct in lineage matrix, update ending time and extinction time
      e$lineages[lin_id, c("end_time", "status", "spec_or_ext_ct")] <- c(e$t, -2, e$t) 
      e$traits[[lin_id]]$status <- -2 #update lineage status to extinct in trait list
      
      daughter <- e$traits[[lin_id]]$relative_id #find daughter lineage of extinct lineage
      if (e$lineages[daughter, "status"] == -1) { #if daughter is incipient & alive, becomes good:
        # update daughter lineage status to good in lineage matrix, update speciation completion times
        e$lineages[daughter, c("status", "spec_or_ext_ct", "spec_ct")] <- c(1, e$t, e$t) 
        e$traits[[daughter]]$status <- 1 #update daughter lineage status to good in trait list
      }
    }
  }
  return(e)
}



#calculate absolute value of difference from trait value of parental lineage
calc_parental_abs_trait_diff <- function(trait) {
  abs(current_trait_value(trait) - trait$parent_trait_value)
}

# calculate rate of completion of speciation for species - dependent on distance from parent
calc_lamda2 <- function(trait, tau0, beta) {
  tau0 * exp(beta * (calc_parental_abs_trait_diff(trait)) ^ 2)
}

calc_mui <- function(trait, trait_diff, alpha1,
                     mui0, alpha3, mui1, muibg) {
  # calculate rate of competitive-dependent extinction for species - depends on distance from lineages
  alpha1 * mui0 * sum(exp(-alpha1 * trait_diff ^ 2)) +
    # calculate rate of selective-dependent extinction for species - depends on distance from optima
    alpha3 * mui1 * (1 - sum(exp(-alpha3 * calc_parental_abs_trait_diff(trait) ^ 2))) + 
    # add background extinction rate
    muibg 
}

# calculate rate of competitive-dependent extinction for species - depends on distance from lineages
# calculate rate of selective-dependent extinction for species - depends on distance from optima
# add background extinction rate
calc_mu <- function(trait_diff, alpha1, alpha3, mu0, mu1, mubg, trait_diff_opt) {
 mu <- alpha1 * mu0 * sum(exp(-alpha1 * (trait_diff) ^ 2)) + 
    alpha3 * mu1 * (1 - sum(exp(-alpha3 * trait_diff_opt ^ 2))) + 
    mubg 
 
 if (mu0 != 0 & mu == 0) { #captures the case when there's only one lineage alive & extinction is activated
   mu <- 0.02 * pars$mu0 #when alone, a lineage has basal extinction rate (equal to having infinite distance with neighbors & no selection)
 }
 
 mu
}

calc_diff_trait_opt <- function(trait, opt) {
  current_trait_value(trait) - opt
}

# Calculate vector of speciation completion & extinction probabilities
calc_probs <- function(mu, lambda) {
  c(lambda/(lambda + mu), mu/(lambda + mu)) 
}

# randomly determine a probability that lineage does not remain in same state
cease_state <- function(lambda, mu, step_size) {
  runif(1) <= (lambda + mu) * step_size
}


# extract lineage trait difference vector from trait diff matrix
trait_diff_extr <- function(trait_diff_m, lin_id) {
  trait_diff_m[as.character(lin_id), colnames(trait_diff_m) != as.character(lin_id)]
}


evolve_trait <- function(x, alpha2, m, s, opt, sign_t, diff_t, bounds) {
  x + alpha2 * m * sum(sign_t * exp(-alpha2 * (diff_t)^2)) + #add competitive element where the trait value 
    # is pushed away from trait values of other lineages based on the differences - the closer two traits are, 
    # the stronger competition will drive them away
    rnorm(1, 0, s) + # add Brownian Motion element
    boundary_effect(x, bounds)
}

#calculate distance of current trait value from bounds & calculate the bound effect
# this is added to the simulated trait value to make sure the trait is kept inside the bounds  
boundary_effect <- function(x, bounds) {
  3 * sum(sign(x - bounds) * exp(-2 * (x - bounds)^2)) 
}

# add selection element where the trait value is pulled towards the optima based 
# on how far away it is - the farther from the optimum it is, the stronger the pull
add_ou_effect <- function(x1, x0, opt, alpha4, step_size) {
  for (j in seq_along(opt)){
    x1 <- x1 + alpha4[j] * (opt[j] - x0) * step_size
  }
  x1
}



cs_current_trait_diff <- function(traits, active_lineages) {
  n_lineages = length(active_lineages)
  
  # fill in matrix by calculating differences in trait values between each pairwise combination of lineages
  trait_diff <- matrix(0, n_lineages, n_lineages, 
                       dimnames = list(active_lineages, active_lineages)) 
  
  
  for (i in seq_along(active_lineages)) {
    for (j in seq_along(active_lineages)){
      
      lin_i <- active_lineages[i]
      lin_j <- active_lineages[j]
      
      trait_diff[i, j] <- tail(traits[[lin_i]]$trait_values, 1) - 
        tail(traits[[lin_j]]$trait_values, 1) 
    } 
  }
  
  diag(trait_diff) <- NA
  
  trait_diff
}


cs_trait_diff_sgn <- function(trait_diff) {
  # extract signs from matrix: which trait in pair is lower and which is higher
  trait_diff_sgn <- sign(trait_diff) 
  
  # if any pair of lineages have identical trait values: 
  # assigns a random sign to each member of a pair of identical lineages, to
  # make sure they are still driven further away from each other
  if (any(trait_diff_sgn == 0, na.rm = T)) { 
    trait_diff_sgn[lower.tri(trait_diff_sgn)] <- NA # turn lower triangle into NA to remove duplicates for rest of loop
    eq_ind <-  trait_diff_sgn == 0 & !is.na(trait_diff_sgn) # identify which lineages make up identical lineage pairs
    trait_diff_sgn[eq_ind] <- sign(rnorm(sum(eq_ind, na.rm = TRUE))) #assign random sign to lineage pair
    
    trait_diff_sgn[lower.tri(trait_diff_sgn)] <- -trait_diff_sgn[upper.tri(trait_diff_sgn)] #reassign lower triangle
  } 
  trait_diff_sgn
}

