#' Run CAMPSITE simulation
#' 
#' Run simulation using the Competition And Multiple-Peak Selection Integrated Trait Evolution (CAMPSITE) model
#' @param pars an object of class `cs_pars`.
#' @param ou a named list of length 2 containing a vector of `optima` and `alpha4` values.
#' @param root.value numeric. The initial trait value of the root node. Must be within the specified bounds.
#' @param age.max numeric. maximum age of simulation
#' @param age.ext age of extinction
#' @param step_size size of each simulation step
#' @param bounds numeric vector of length 2. Trait value bounds
#' @param plot logical. Whether to plot resulting tree
#' @param ylims ylim
#' @param full_results logical. Whether to return full results, including lineage and
#'
#' @return an object of class `cs_sim_results` containing the following:
#' - `trees` trees resulting from the simulated lineage evolution. Each tree element contains the tree, an object of class `phylo` (`tree`),
#'   and a numeric vector of tip values (`tips`).
#'     - `all` tree of all lineages
#'     - `gsp_fossil` tree that excludes incomplete lineages, i.e. lineages that were 
#'     incipient at the end of the simulation or went extinct before completing speciation.
#'     - `gsp_extant` tree including only extant species. Both `tree` and `tips` will 
#'       be `NULL` if the process dies and there are no extant species at the end of the simulation.
#'     
#' @export
#'
#' @examples
#' \dontrun{
#' pars <- cs_pars()
#' 
#' cs_simulate(pars, ou = list(opt = NULL, alpha4 = NULL), root.value = 0, age.max = 10, 
#' age.ext = NULL, step_size = 0.01, bounds = c(-Inf, Inf), 
#' plot = TRUE, ylims = NULL, full_results = TRUE) 
#' }
#' 
cs_simulate <- function (pars, ou = list(opt = NULL, alpha4 = NULL), root.value = 0, age.max = 50, 
                          age.ext = NULL, step_size = 0.01, bounds = c(-Inf, Inf), 
                          plot = TRUE, ylims = NULL, full_results = FALSE) 
{
  
  validate_bounds(bounds, root.value)
  ou <- validate_ou(ou)
  
  pars$s = sqrt(pars$sig2 * step_size) #variance of BM per step size
  pars$m = pars$m * step_size #relative contribution of competition with respect to BM per step size
  
  process_dead <- FALSE #parameter signifying whether simulation is still running
  
  e <- list(traits = list(lin1 = new_trait(1, 1, 2, NA, root.value),
                                lin2 = new_trait(2, -1, 1, root.value, root.value)), 
                  lineages = new_lineages(
                    new_lineage(c(1, 0, 0, -1, 1, 0, 0)), 
                    new_lineage(c(1, 0, 0, -1, -1, 0, NA))
                  ),
                  t =  0 + step_size,
                  active_lineages = c("lin1", "lin2"),
                  step_count = 1L
  )
  
  # Initialise progress bar
  usethis::ui_info("Simulation initiated \n")
  pb <- progress::progress_bar$new(total = ceiling(age.max / step_size), show_after = 0,
                                   clear = FALSE)
  #pb$tick(0)
  #t must go one step beyond age.max to ensure simulation of last step. /2 is to avoid numerical precision issues
  while ((age.max - e$t) > -step_size / 2) { 
    
    # list of differences - each element will store the differences 
    # in trait value between each lineage and all other co-occurring lineages
    e$trait_diff_m <- cs_current_trait_diff(e$traits, e$active_lineages)
    
    
    
    #simulate new trait value for each lineage:
    e$traits <- cs_sim_next_trait_value(traits = e$traits, e$active_lineages, e$trait_diff_m, ou = ou, 
                                        pars = pars, bounds = bounds, step_size = step_size)
    
    for (lineage in e$active_lineages) {
      
      if (is_incipient(e$lineages, lin_id = lineage)) { #if lineage is incipient
        e <- cs_evolve_inc_lineage(pars, lin_id = lineage, step_size, e$t, e)
      }
      if (is_good(e$lineages, lin_id = lineage)) { #if lineage is good
        e <- cs_evolve_good_lineage(pars, lin_id = lineage, step_size, e$t, ou, e)
      }
    }
    e$active_lineages <- lineages_active(e$lineages)
    e$step_count <- e$step_count + 1L 
    e$t <- e$t + step_size 
    
    # kill process if number of lineages has dropped to zero
    if (length(e$active_lineages) == 0) {
      usethis::ui_info("Process dead. No active lineages left \n")
      process_dead <- TRUE
      (break)()
    }
    
    pb$tick()
  }
  usethis::ui_done("Simulation complete \n")
  # ---- BUILD TREES ----
  # retract step size to return back to current time
  t_end <- e$t - step_size 
  trees <- cs_build_trees(e$lineages, e$traits, t_end)

  usethis::ui_done("Building tree complete \n")
  
  if (plot) {
    
    if (length(ylims) < 2) { 
      # make sure limits on y axis aren't too small
      ylims <- NULL
    }
    # plot the simulation
    cs_plot_trait_evolution(traits = e$traits, lineages = e$lineages, 
                            t = t_end, step_size = step_size, ylims = ylims) 
  }
  
  if (full_results == TRUE) {
    
    return(
      new_sim_result(trees = trees, t_end = t_end, step_size = step_size, pars = pars,
                     process_dead = process_dead, lineages = e$lineages, traits = e$traits, 
                     full_results = TRUE, end_active_lineages = TRUE) 
    ) 
  } else {
    return(
      new_sim_result(trees = trees, t_end = t_end, step_size = step_size, pars = pars,
                     process_dead = process_dead, lineages = NULL, traits = NULL, 
                     full_results = FALSE, end_active_lineages = TRUE) 
    ) 
  }
}
