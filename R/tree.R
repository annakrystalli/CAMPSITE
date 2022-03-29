# BUILD TREES & GET TIPS TRAIT VECTORS
#' @export
cs_build_trees <- function(lineages, traits, t) {
  # sort active lineages in increasing order, after including extinct lineages
  edge_lineages <- lineages_edges(lineages)
  edges_mat <- cs_edge_matrix(lineages)
  n_tips <- lineages_tips_n(lineages)
  
  # update ending time for all incipient lineages at end of process
  lineages[lineages[, "end_time"] < 0, "end_time"] <- t 
  
  # calculate branch lengths per lineage (ending time - starting time)
  lin_length <- round(lineages[, "end_time"] - lineages[, "start_time"], 2) 
  
  # generate the phylogeny:
  tree <- list(edge = edges_mat,
               edge.length = lin_length, 
               Nnode = (n_tips - 1),
               tip.label = paste("t", as.character(1:n_tips), sep = ""))
  
  
  class(tree) <- "phylo" #change class of tree
  
  tip_values <- NULL
  # extract tip values as trait values at the last time step
  tip_values <- sapply(edge_lineages,  FUN = function(x){current_trait_value(traits[[x]])}) 
  names(tip_values) <- tree$tip.label #name vector of tip values
  
  # create vector of all lineages that were incipient at the end of the process, 
  # or went extinct before completing speciation ('incomplete lineages')
  incomplete_todrop <- rownames(lineages)[is_incipient(lineages) | is_extinct(lineages) & 
                                     is.na(lineages[, "spec_ct"])]

  # extract node numbers of incomplete lineages                
  tip_ids <- edges_mat[incomplete_todrop, "descendant_node"] 
  # prune tree to remove incomplete lineages
  tree_gsp_fossil <- ape::drop.tip(tree, tip = tip_ids) 
  # prune vector of tip values to remove incomplete lineages
  tip_values_gsp_fossil <- tip_values[names(tip_values) %in% 
                                        tree_gsp_fossil$tip.label] 
  
  # create vector of all lineages that are incipient at end of process or went 
  # extinct during the process ('extinct lineages')
  extinct_todrop <- rownames(lineages)[is_incipient(lineages) | is_extinct(lineages)]
  
  # if all tips are 'extinct' (extinct or incipient)
  if (length(extinct_todrop) == n_tips) { 
    usethis::ui_info("Process dead. All tips extinct or incipient at end of process. No extant tree available \n")
    tree_gsp_extant <- NULL
    tip_values_gsp_extant <- NULL
    process_dead <- TRUE
  }
  else {
    tip_ext_ids <- edges_mat[extinct_todrop, "descendant_node"]  
    tree_gsp_extant <- ape::drop.tip(tree, tip = tip_ext_ids) 
    tip_values_gsp_extant <- tip_values[names(tip_values) %in% 
                                          tree_gsp_extant$tip.label] 
    process_dead <- FALSE
    # tree_gsp_fossil <- ape::read.tree(text = ape::write.tree(tree_gsp_fossil)) #save fossil tree (no incomplete lineages)
    # tree_gsp_extant <- ape::read.tree(text = ape::write.tree(tree_gsp_extant)) #save extant tree (no extinct lineages)
  }
  #tree <- read.tree(text = write.tree(tree)) #save full tree
  
  return(new_tree_results(tree, tip_values, 
                          tree_gsp_fossil, tip_values_gsp_fossil,
                          tree_gsp_extant, tip_values_gsp_extant))
}

#' @export
cs_edge_matrix <- function(lineages) {

  is_tip <- function(x) {
    x[, "descendant_node"] == 0
  }
  
  is_not_tip <- function(x) {
    x[, "descendant_node"] != 0
  }
  
  # matrix of parental and descendant nodes for each lineage
  edges_mat <- lineages[, c("parent_node", "descendant_node")]
  
  n_tips <- lineages_tips_n(lineages)
  stopifnot(sum(is_tip(edges_mat)) == n_tips)
  
  # update the node numbers to comply with phytools:
  # update parental node numbers to number sequentially from 1 + highest tip number
  edges_mat[, "parent_node"] <- as.integer(edges_mat[, "parent_node"] + n_tips)
  # update ancestral node numbers for all ancestral lineages to number sequentially from 1 + highest tip number
  edges_mat[is_not_tip(edges_mat), "descendant_node"] <- as.integer(edges_mat[is_not_tip(edges_mat), "descendant_node"] + n_tips)
  # change numbers of all tips to number from 1 to number of tips
  edges_mat[is_tip(edges_mat), "descendant_node"] <- as.integer(1:n_tips) 

  return(edges_mat)
  
}