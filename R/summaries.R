#' Summarise results of CAMPSITE simulation
#'
#' @param x results object of class `"cs_sim_results"` output from `cs_simulate()` 
#' or list of `"cs_sim_results"` objects.
#' @param max_trait_len length each trait vector should be expanded by. 
#'
#' @return an object (if `x` is `"cs_sim_results"` object) or list objects 
#' (if `x` is list of `"cs_sim_results"` objects) of class `cs_result_summaries` containing:
#' - `lineages`: object of class `cs_lineages`.
#' - `VAR`: variance of trait values across all lineages at each time step
#' - `MNND`: Mean Nearest Neighbour Distance of trait values across all lineages at each time step.
#' - `VNND`: Variance of Nearest Neighbour Distance of trait values across all lineages at each time step.
#' - `tip_traits`: tip trait values of extant lineages
#' - `Nnode_extant`: extant tree number of nodes
#' - `tree_extant`: object of class `phylo` containing extant lineages only
#' - `tree_fossil`: object of class `phylo` containing tree that excludes incomplete lineages
#' @export
cs_summarise_results <- function(x, max_trait_len = NULL){
  
  if(inherits(x, "cs_sim_results")){
    return(cs_summarise_result(x, max_trait_len))
  }

  res_class <- sapply(x, FUN = function(i){inherits(i, "cs_sim_results")})
  res_class_n <- sum(res_class)
  
  
  if(sum(!res_class) == length(x)){
    usethis::ui_stop("No objects of class cs_sim_results in {usethis::ui_field('x')}")
  }
  
  if(res_class_n != length(x)){
    usethis::ui_warn("{res_class_n} of {length(x)} elements of {usethis::ui_field('x')} not objects of class cs_sim_results. Ignored")
  }
  
  return(
    lapply(x[res_class], function(i){
      cs_summarise_result(i, max_trait_len)})
  )
}


cs_summarise_result <- function(x,  max_trait_len = NULL) {
  checkmate::assert_class(x, "cs_sim_results")
  
  traits_mat <- traits_compile_values_matrix(x$traits,
                                             max_trait_len = max_trait_len)
  
  summaries <- list(
    lineages = x$lineages,
    VAR = apply(traits_mat, 2, stats::var, na.rm = TRUE),
    MNND = apply(traits_mat, 2, MNND),
    VNND = apply(traits_mat, 2, VNND),
    tip_traits = x$trees$gsp_extant$tips,
    Nnode_extant = x$trees$gsp_extant$tree$Nnode,
    tree_extant = x$trees$gsp_extant$tree,
    tree_fossil = x$trees$gsp_fossil$tree,
    competition = x$pars$alpha1,
    selection = x$pars$alpha3,
    t_end = x$t_end,
    step_size = x$step_size,
    replicate = x$replicate
    
  )
  ## Set the name for the class
  class(summaries) <- append(class(summaries),"cs_result_summaries")
  return(summaries)
}


MNND <- function(x){
  a <- x[!is.na(x)]
  if(length(a) > 1){
    a <- as.matrix(stats::dist(a))
    diag(a) <- NA
    return(mean(apply(a, 2, min, na.rm = TRUE)))
  }
  else(return(0))
}

VNND <- function(x){
  a <- x[!is.na(x)]
  if(length(a) > 1){
    a <- as.matrix(stats::dist(a))
    diag(a) <- NA
    return(stats::var(apply(a, 2, min, na.rm = TRUE)))
  }
  else(return(0))
}

sumBL <- function(x){
  tree <- x$tree_fossil
  max.bl <- max(adephylo::distRoot(tree, method = "patristic"))
  
  bl <- numeric()
  for (i in 0:(max.bl - 1)){
    if (i < round(max.bl - 1)) {
      bl[[i + 1]] <- sum(adephylo::distRoot(
        paleotree::timeSliceTree(tree, i, plot = F), method = "patristic")) -
        sum(adephylo::distRoot(suppressWarnings(paleotree::timeSliceTree(tree, i + 1, plot = F)), 
                               method = "patristic"))
    } else {
      bl[[i+1]] <- sum(adephylo::distRoot(suppressWarnings(paleotree::timeSliceTree(tree, i, plot = F)),
                                          method = "patristic")) + 1*(ceiling(max.bl)-max.bl)
    }
  }
  
  bl[is.na(bl)] <- 1
  
  return(data.frame(time_bin = seq_along(bl),
                    branch_len = rev(bl)))
}


