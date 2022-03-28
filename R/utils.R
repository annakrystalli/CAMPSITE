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
    VAR = apply(traits_mat, 2, var, na.rm = TRUE),
    MNND = apply(traits_mat, 2, MNND),
    VNND = apply(traits_mat, 2, VNND),
    tip_traits = x$trees$gsp_extant$tips,
    Nnode_extant = x$trees$gsp_extant$tree$Nnode,
    tree_extant = x$trees$gsp_extant$tree,
    tree_fossil = x$trees$gsp_fossil$tree,
    competition = x$pars$alpha1,
    selection = x$pars$alpha3,
    t_end = x$t_end,
    step_size = x$step_size
    
  )
  ## Set the name for the class
  class(summaries) <- append(class(summaries),"cs_result_summaries")
  return(summaries)
}


MNND <- function(x){
  a <- x[!is.na(x)]
  if(length(a) > 1){
    a <- as.matrix(dist(a))
    diag(a) <- NA
    return(mean(apply(a, 2, min, na.rm = TRUE)))
  }
  else(return(0))
}

VNND <- function(x){
  a <- x[!is.na(x)]
  if(length(a) > 1){
    a <- as.matrix(dist(a))
    diag(a) <- NA
    return(var(apply(a, 2, min, na.rm = TRUE)))
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
        sum(adephylo::distRoot(paleotree::timeSliceTree(tree, i + 1, plot = F), 
                               method = "patristic"))
    } else {
      bl[[i+1]] <- sum(adephylo::distRoot(paleotree::timeSliceTree(tree, i, plot = F),
                                          method = "patristic")) + 1*(ceiling(max.bl)-max.bl)
    }
  }
  
  bl[is.na(bl)] <- 1
  
  return(data.frame(time_bin = seq_along(bl),
                    branch_len = rev(bl)))
}


extractSignal <- function(model, method = c("K", "lambda")){
  tree <- model$tree_extant
  tips <- model$tip_traits
  tips <- tips[tree$tip.label]
  if (method == "K")
  {K <- phylosig(tree, tips, "K") 
  return(as.numeric(K))} else
    if (method == "lambda")
    {lambda <- phylosig(tree, tips, "lambda")
    return(lambda$lambda)}
}

fitTraitModels <- function(model){
  models <- numeric()
  tree <- model$tree_extant
  tree$edge.length <- sapply(tree$edge.length, function(x) x + rnorm(1, sd = 0.001))
  tips <- as.matrix(model$tip_traits)
  models[[1]] <- transformPhylo.ML(tips, tree, model = "BM")$AICc
  models[[2]] <- transformPhylo.ML(tips, tree, model = "OU")$AICc
  models[[3]] <- transformPhylo.ML(tips, tree, model = "ACDC")$AICc
  models[[4]] <- fit_t_comp(tree, tips[,1], model = "MC")$aicc
  models[[5]] <- fit_t_comp(tree, tips[,1], model = "DDlin")$aicc
  models[[6]] <- fit_t_comp(tree, tips[,1], model = "DDexp")$aicc
  names(models) <- c("BM", "OU", "ACDC", "MC", "DDlin", "DDexp")
  if (names(sort(models))[[1]] == "ACDC"){
    if (transformPhylo.ML(tips, tree, model = "ACDC")$ACDC[,3] < 0) {names(models)[[3]] <- "ACDCdec"} 
    else {names(models)[[3]] <- "ACDCinc"}
  }
  return(names(sort(models))[[1]])
}

fitDivModels <- function(model){
  models <- numeric()
  tree <- model$tree_extant
  f.null <- function(t,y){0}
  f.cst <- function(t,y){y[1]}
  f.lin <- function(t,y){abs(y[1] + y[2] * t)}
  f.exp <- function(t,y){y[1] * exp(y[2] * t)}
  models[[1]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.cst, f.mu = f.null, lamb_par = c(0.09), mu_par = c(), f = 1, cst.lamb = T, fix.mu = T, cond = "stem")$aicc
  models[[2]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.lin, f.mu = f.null, lamb_par = c(0.09, 0.001), mu_par = c(), f = 1, fix.mu = T, cond = "stem")$aicc
  models[[3]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.exp, f.mu = f.null, lamb_par = c(0.05, 0.01), mu_par = c(), f = 1, expo.lamb = T, fix.mu = T, cond = "stem")$aicc
  models[[4]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.cst, f.mu = f.cst, lamb_par = c(0.09), mu_par = c(0.005), f = 1, cst.lamb = T, cst.mu = T, cond = "stem")$aicc
  models[[5]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.lin, f.mu = f.cst, lamb_par = c(0.09, 0.001), mu_par = c(0.005), f = 1, cst.mu = T, cond = "stem")$aicc
  models[[6]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.exp, f.mu = f.cst, lamb_par = c(0.05, 0.01), mu_par = c(0.005), f = 1, expo.lamb = T, cst.mu = T, cond = "stem")$aicc
  models[[7]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.cst, f.mu = f.lin, lamb_par = c(0.09), mu_par = c(0.005, 0.0001), f = 1, cst.lamb = T, cond = "stem")$aicc
  models[[8]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.lin, f.mu = f.lin, lamb_par = c(0.09, 0.001), mu_par = c(0.005, 0.0001), f = 1, cond = "stem")$aicc
  models[[9]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.exp, f.mu = f.lin, lamb_par = c(0.05, 0.01), mu_par = c(0.005, 0.0001), f = 1, expo.lamb = T, cond = "stem")$aicc
  models[[10]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.cst, f.mu = f.exp, lamb_par = c(0.09), mu_par = c(0.0035, 0.001), f = 1, cst.lamb = T, expo.mu = T, cond = "stem")$aicc
  models[[11]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.lin, f.mu = f.exp, lamb_par = c(0.09, 0.001), mu_par = c(0.0035, 0.001), f = 1, expo.mu = T, cond = "stem")$aicc
  models[[12]] <- fit_bd(phylo = tree, tot_time = 50, f.lamb = f.exp, f.mu = f.exp, lamb_par = c(0.05, 0.01), mu_par = c(0.0035, 0.001), f = 1, expo.lamb = T, expo.mu = T, cond = "stem")$aicc
  names(models) <- c("Bcst", "Blin", "Bexp", "BcstDcst", "BlinDcst", "BexpDcst", "BcstDlin", "BlinDlin", "BexpDlin", "BcstDexp", "BlinDexp", "BexpDexp")
  return(names(sort(models))[[1]])
}
