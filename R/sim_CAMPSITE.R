sim_CAMPSITE <- function (pars, ou = list(opt = NULL, alpha4 = NULL), root.value = 0, age.max = 50, 
                          age.ext = NULL, step.size = 0.01, bounds = c(-Inf, Inf), 
                          plot = TRUE, ylims = NULL, full.sim = FALSE) 
{
  pars <- CAMPSITE::cs_pars()
  load("attic/dev-setup.RData")
  
  validate_bounds(bounds, root.value)
  ou <- validate_ou(ou)
  
  pars$s = sqrt(pars$sig2 * step.size) #variance of BM per step size
  pars$m = pars$m * step.size #relative contribution of competition with respect to BM per step size
  
  process_dead = FALSE #parameter signifying whether simulation is still running
  
  
  # e <- rlang::new_environment(list(traits = list(lin1 = new_trait(1, 1, 2, NA, root.value),
  #                                                lin2 = new_trait(2, -1, 1, root.value, root.value)), 
  #                                  lineages = new_lineages(
  #                                    new_lineage(c(1, 0, 0, -1, 1, 0, 0)), 
  #                                    new_lineage(c(1, 0, 0, -1, -1, 0, NA))
  #                                    )
  #                                  )
  #                             )
  
  e <- rlang::env(traits = list(lin1 = new_trait(1, 1, 2, NA, root.value),
                                lin2 = new_trait(2, -1, 1, root.value, root.value)), 
                  lineages = new_lineages(
                    new_lineage(c(1, 0, 0, -1, 1, 0, 0)), 
                    new_lineage(c(1, 0, 0, -1, -1, 0, NA))
                  ),
                  dead_lin = NULL, #reset vector of non-active lineages
                  born_lin = NULL, #reset vector of new lineages
                  t =  0 + step.size,
                  active_lineages = c("lin1", "lin2"),
                  n_good = 2L,
                  step_count = 1L
  )
  
  #t <- 0 + step.size #add step size to current time to simulate trait values and 
  # extinction/speciation for the next time step. 
  # The next bit of code will run until maximum age is reached:
  while ((age.max - e$t) > -step.size/2) { #t must go one step beyond age.max to ensure simulation of last step. /2 is to avoid numerical precision issues
    
    e$trait_diff_m <- cs_current_trait_diff(e$traits, e$active_lineages)
    #list of differences - each element will store the differences 
    #in trait value between each lineage and all other co-occurring lineages
    #diff_me = list()
    
    #simulate new trait value for each lineage:
    e$traits <- cs_sim_next_trait_value(traits = e$traits, e$active_lineages, e$trait_diff_m, ou, 
                                        params, bounds)
    
    for (i in e$active_lineages) {
      
      
      if (is_incipient(e$lineages, lin_id = i)) { #if lineage is incipient
        cs_evolve_inc_lineage(pars, lin_id = i, step.size, e$t)
      }
      if (is_good(e$lineages, lin_id = i)) { #if lineage is good
        cs_evolve_good_lineage(pars, lin_id = i, step.size, e$t)
      }
    }
    e$active_lineages <- lineages_active(e$lineages)
    e$step_count <- e$step_count + 1L #update step count to next step
    e$t <- e$t + step.size #update current time to next step
    cat("\r", "time:", t) #print out progress
    if (lenth(e$active_lineages) == 0) { #kill process if number of lineages has dropped to zero
      print("process died")
      e$process_dead = TRUE
      (break)()
    }
  }
  
  
  t = t - step.size #retract step size to return back to current time
  
  ##BUILD TREES & GET TIPS TRAIT VECTORS
  #complete process tree
  row.names(lineages) <- NULL
  colnames(lineages) <- NULL
  edges_mat <- lineages[, 1:2] #matrix of parental and descendant nodes for each lineage
  active_lineages <- sort(c(active_lineages, which(lineages[, 5] == -2))) #sort active lineages in increasing order, after including extinct lineages
  n_tips = length(active_lineages) #how many tips in phylogeny
  #update the node numbers to comply with phytools:
  edges_mat[, 1] <- edges_mat[, 1] + n_tips #update parental node numbers to number sequentially from 1 + highest tip number
  edges_mat[, 2][which(edges_mat[, 2] != 0)] <- edges_mat[, 2][which(edges_mat[, 2] != 0)] + n_tips #update ancestral node numbers for all ancestral lineages to number sequentially from 1 + highest tip number
  edges_mat[, 2][which(edges_mat[, 2] == 0)] <- 1:n_tips #change numbers of all tips to number from 1 to number of tips
  edg1 <- as.integer(edges_mat[, 1]) #vector of parental nodes per lineage
  edg2 <- as.integer(edges_mat[, 2]) #vector of node numbers per lineage (tip or internal)
  edges_mat <- cbind(edg1, edg2) #create new matrix of nodes
  dimnames(edges_mat) <- NULL
  lineages[, 4][which(lineages[, 4] == -1)] <- t #update ending time for all incipient lineages at end of process
  lin_length <- round(lineages[, 4] - lineages[, 3], 2) #calculate branch lengths per lineage (ending time - starting time)
  #generate the phylogeny:
  tree <- list(edge = edges_mat,
               edge.length = lin_length, 
               Nnode = (n_tips - 1),
               tip.label = paste("t", as.character(1:n_tips), sep = ""))
  class(tree) <- "phylo" #change class of tree
  tip_values = NULL
  tip_values <- sapply(traits[active_lineages], function(x) (x[length(x)])) #extract tip values as trait values at the last time step
  names(tip_values) <- tree$tip.label #name vector of tip values
  isp_todrop <- c(which(lineages[, 5] == -1),
                  which(lineages[, 5] == -2 & is.na(lineages[, 7]))) #create vector of all lineages that were incipient at the end of the process, or went extinct before completing speciation ('incomplete lineages')
  tip_ids <- edges_mat[isp_todrop, 2] #extract node numbers of incomplete lineages
  tree_gsp_fossil <- drop.tip(tree, tip = tip_ids) #prune tree to remove incomplete lineages
  tip_values_gsp_fossil <- tip_values[names(tip_values) %in% 
                                        tree_gsp_fossil$tip.label] #prune vector of tip values to remove incomplete lineages
  extinct_todrop <- c(which(lineages[, 5] == -1),
                      which(lineages[, 5] == -2)) #create vector of all lineages that are incipient at end of process or went extinct during the process ('extinct lineages')
  if (length(extinct_todrop) == n_tips) { #if all tips are 'extinct' (extinct or incipient)
    print("process died")
    process_dead = T
  }
  if (process_dead == T) {
    tree_gsp_extant <- "process died"
    tip_values_gsp_extant <- "process died"
  }
  else {
    tip_ext_ids <- edges_mat[extinct_todrop, 2] #extract node numbers of extinct lineages
    tree_gsp_extant <- drop.tip(tree, tip = tip_ext_ids) #prune tree to remove extinct lineages
    tip_values_gsp_extant <- tip_values[names(tip_values) %in% 
                                          tree_gsp_extant$tip.label] #prune vector of tip values to remove extinct lineages
    tree_gsp_fossil <- read.tree(text = write.tree(tree_gsp_fossil)) #save fossil tree (no incomplete lineages)
    tree_gsp_extant <- read.tree(text = write.tree(tree_gsp_extant)) #save extant tree (no extinct lineages)
  }
  tree <- read.tree(text = write.tree(tree)) #save full tree
  res <- list(all = list(tree = tree, tips = tip_values), 
              gsp_fossil = list(tree = tree_gsp_fossil, tips = tip_values_gsp_fossil), 
              gsp_extant = list(tree = tree_gsp_extant, tips = tip_values_gsp_extant)) #save simulation outputs
  if (full.sim == T) {
    res$all$trait_mat <- traits #save trait list to simulation output
    res$all$lin_mat <- lineages #save lineage matrix to simulation output
  }
  if (plot) {
    plotSimu <- function(traitmat, linmat, step_size, ylims = NULL) {
      #plots a simulation, with incipient lineages in red, good in black
      steps <- max(linmat[, 4])/step_size #calculate number of steps in simulation
      max_trait <- NULL
      min_trait <- NULL
      for (i in 1:nrow(linmat)) {
        max_trait <- c(max_trait,
                       max(traitmat[[i]][-(1:4)], 
                           na.rm = T)) #extract maximal trait values for each lineage
        min_trait <- c(min_trait,
                       min(traitmat[[i]][-(1:4)], 
                           na.rm = T)) #extract minimal trait values for each lineage
      }
      if (is.null(ylims)) {
        plot(1,
             type = "n",
             xlim = c(1, steps),
             ylim = c(min(min_trait), 
                      max(max_trait)),
             ylab = "trait value",
             xlab = "Time") #generate plot area
      }
      else {
        plot(1, type = "n",
             xlim = c(1, steps),
             ylim = ylims, 
             ylab = "trait value",
             xlab = "Time") #generate plot area with pre-determined limits for trait space
      }
      completion <- linmat[, 6] #extract vector of speciation or extinction times for each lineage
      #handle and plot each lineage based on its end point:
      for (i in 1:length(completion)) {
        if (is.na(completion[i])) { #extinct incipient lineages
          lines(x = (5:length(traitmat[[i]])), #branch length
                y = traitmat[[i]][5:length(traitmat[[i]])], #trait values
                col = "red", #incipient lineages in red
                cex = 0.5)
        }
        else if (linmat[i, 5] == -2 &&
                 !is.na(linmat[i, 7]) &&
                 completion[i] != linmat[i, 7]) { #extinct good lineages
          lines(x = 5:((linmat[i, 7]/step_size) + 5), #branch length of incipient stage
                y = traitmat[[i]][5:((linmat[i, 7]/step_size) + 5)], #trait values of incipient stage
                col = "red", #incipient stage in red
                cex = 0.5)
          lines(x = ((linmat[i, 7]/step_size) + 5):length(traitmat[[i]]), #branch length of good stage
                y = traitmat[[i]][((linmat[i, 7]/step_size) + 5):length(traitmat[[i]])], #trait value of good stage
                col = "black", #good stage in black
                cex = 0.5)
        }
        else { #extant good and incipient lineages
          lines(x = 5:((completion[i]/step_size) + 5), #branch length of incipient stage
                y = traitmat[[i]][5:((completion[i]/step_size) + 5)], #trait values of incipient stage
                col = "red", #incipient stage in red
                cex = 0.5)
          if (length(traitmat[[i]]) > (round(completion[i]/step_size) +  5)) { #if completed speciation before end of process
            lines(x = ((completion[i]/step_size) + 5):length(traitmat[[i]]), #branch length of good stage
                  y = traitmat[[i]][((completion[i]/step_size) + 5):length(traitmat[[i]])], #trait values of good stage
                  col = "black", #good stage in black
                  cex = 0.5)
          }
        }
      }
    }
    if (length(ylims) < 2) { #make sure limits on y axis aren't too small
      ylims = NULL
    }
    plot <- plotSimu(traits, lineages, step.size, ylims) #plot the simulation
  }
  return(res) #return simulation output
}





runCAMPSITE <- function(comp, selec){
  
  
  
  cores = detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  set.seed(42)
  sim <- foreach(i=1:100, .packages = c("ape")) %dopar% {
    repeat {
      model <- sim_CAMPSITE(pars = c(0.25,
                                     0.01,
                                     0.6,
                                     0.2,
                                     0.2,
                                     0.01,
                                     0.4,
                                     0.4,
                                     0.02,
                                     comp,
                                     comp,
                                     selec,
                                     0.5,
                                     20),
                            ou = list(c(0),
                                      c(selec)),
                            bounds = c(-Inf,
                                       Inf),
                            root.value = 0,
                            age.max = 50,
                            step.size = 0.01,
                            full.sim = T,
                            plot = F
      )
      if (!isTRUE(model$gsp_extant$tree == c("process died")) && model$gsp_extant$tree$Nnode > 4) {break}
    }
    model
  }
  stopCluster(cl)
  saveRDS(sim, paste("output/sims/comp", comp, "_selec", selec, ".rds", sep = ""))
  rm(sim)
}
