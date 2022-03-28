#' Create a trait object
#'
#' @param x a vector of four or more values.
#'
#' @return
#' An object of class "trait" with the following values
#' @param lin_id lineage identifying number
#' @param status status (incipient=-1/good=1/extinct=-2)
#' @param relative_id relative node (lineage number). If trait is associated with incipient 
#' lineage, `relative_id` indicates the parent lineage. if trait is associated with 
#' good lineage, `relative_id` indicates the daughter lineage.
#' @param trait_values current trait value of parental lineage, and all consecutive 
#' elements contain trait value at each time step (always `root.value` at the beginning)
#' @examples
#' new_trait(1, 1, 2, c(NA, root.value))
new_trait <- function(lin_id = integer(1), status = integer(1), relative_id = integer(1), 
                      parent_trait_value = numeric(1), trait_values = numeric()) {
  checkmate::assert_numeric(lin_id, len = 1)
  checkmate::assert_numeric(status, len = 1)
  checkmate::assert_numeric(relative_id, len = 1)
  checkmate::assert_numeric(parent_trait_value, len = 1)
  checkmate::assert_numeric(trait_values)
  
  x <- list(as.integer(lin_id), as.integer(status), as.integer(relative_id), 
            parent_trait_value, trait_values)
  names(x) <- c("lin_id", "status", "relative_id", "parent_trait_value", "trait_values")
  
  class(x) <- append(class(x), "cs_trait")
  return(x)
  
}


current_trait_value <- function(x) {
  UseMethod("current_trait_value", x)
}
#' @export
current_trait_value.cs_trait <- function(x){
  return(tail(x$trait_values, 1))
}




min_trait_value <- function(x) {
  UseMethod("min_trait_value", x)
}
#' @export
min_trait_value.cs_trait <- function(x){
  return(min(x$trait_values, na.rm = TRUE))
}



max_trait_value <- function(x) {
  UseMethod("max_trait_value", x)
}
#' @export
max_trait_value.cs_trait <- function(x){
  return(max(x$trait_values, na.rm = TRUE))
}



extract_trait_values <- function(x) {
  UseMethod("extract_trait_values", x)
}
#' @export
extract_trait_values.cs_trait <- function(x){
  return(x$trait_values)
}


traits_compile_values <- function(x) {
  stopifnot(all(sapply(x, inherits, what = "cs_trait")))
  
  sapply(x, extract_trait_values)
}

traits_compile_values_matrix <- function(x, max_trait_len = NULL) {
  checkmate::assert_number(max_trait_len, null.ok = TRUE)
  
  out <- traits_compile_values(x)
  if(is.null(max_trait_len)){
    max_trait_len <- max(sapply(out, length))
  }
  
  out <- lapply(out, function(v, max_trait_len){
    length(v) <- max_trait_len
    return(v)
  }, max_trait_len) %>%
    do.call(rbind, .)
  
}

#' Create a lineage object
#'
#' @param x numeric vector of length 7.
#' @return
#' An object of class "lineage" with the following values
#' - `parent_node` parental node (lineage number)`
#' - `descendant_node` descendent node (0 if tip)
#' - `start_time` starting time
#' - `end_time` ending time (-1 if still active)
#' - `status` status (incipient=-1/good=1/extinct=-2)
#' - `spec_or_ext_ct` speciation completion or extinction time
#' - `spec_ct` speciation completion time (NA if still incipient)

#' @examples
#' new_lineage(c(1, 0, 0, -1, 1, 0, 0))
new_lineage <- function(x = numeric(length = 7L)) {
  checkmate::assert_numeric(x, len = 7)
  x[c(1, 2, 5)] <- as.integer(x[c(1, 2, 5)])
  names(x) <- c("parent_node", "descendant_node", "start_time", "end_time", 
                "status", "spec_or_ext_ct", "spec_ct")
  
  structure(x, class = "cs_lineage")
}



new_lineages <- function(...) {
  
  lapply(list(...), 
         FUN = function(obj){checkmate::assert_class(obj, classes = "cs_lineage")})
  checkmate::assert_true(...length() > 0)
  
  x <- rbind(...)
  rownames(x) <- int_to_lin_id(1:...length())
  
  class(x) <- append(class(x), "cs_lineages")
  
  return(x)
}



lineage_next_parent <- function(x) {
  UseMethod("lineage_next_parent", x)
}

lineage_next_parent.cs_lineages <- function(x){
  return(max(x[, "parent_node"]) + 1)
}



lineage_next <- function(x) {
  UseMethod("lineage_next", x)
}

lineage_next.cs_lineages <- function(x){
  return(max(lin_id_to_int(rownames(x))) + 1)
}

#' @export
lineage_add <- function(x, ...) {
  UseMethod("lineage_add", x)
}

lineage_add.cs_lineages <- function(x, ...){
  checkmate::assert_class(x, classes = "cs_lineages")
  
  lapply(list(...), 
         FUN = function(obj){checkmate::assert_class(obj, classes = "cs_lineage")})
  
  lin_names <- rownames(x)
  lin_max <- max(lin_id_to_int(lin_names))
  lin_names <- c(lin_names, int_to_lin_id(1:...length() + lin_max))
  
  
  y <- rbind(x, ...)
  rownames(y) <- lin_names
  
  structure(y, class = c('matrix', 'array',"cs_lineages"))
  
}



lineage_last_speciated <- function(x) {
  UseMethod("lineage_last_speciated", x)
}

lineage_last_speciated.cs_lineages <- function(x){
  return(rownames(tail(x, 2)))
}


lineage_last_incipient <- function(x) {
  UseMethod("lineage_last_incipient", x)
}

lineage_last_incipient.cs_lineages <- function(x){
  x <- x[lineage_last_speciated(x), ]
  return(rownames(x)[is_incipient(x)])
}


lineage_last_parent <- function(x) {
  UseMethod("lineage_last_parent", x)
}

lineage_last_parent.cs_lineages <- function(x){
  x <- x[lineage_last_speciated(x), ]
  return(rownames(x)[!is_incipient(x)])
}


lineages_dead <- function(x) {
  UseMethod("lineages_dead", x)
}

lineages_dead.cs_lineages <- function(x){
  return(rownames(x)[x[ , "end_time"] > 0])
}


lineages_extinct <- function(x) {
  UseMethod("lineages_extinct", x)
}

lineages_extinct.cs_lineages <- function(x){
  return(rownames(x)[x[ , "status"] < -1])
}


lineages_born <- function(x) {
  UseMethod("lineages_born", x)
}

lineages_born.cs_lineages <- function(x){
  return(rownames(x)[x[ , "start_time"] > 0])
}



lineages_active <- function(x) {
  UseMethod("lineages_active", x)
}

lineages_active.cs_lineages <- function(x){
  return(rownames(x)[x[ , "status"] > -2 & x[ , "end_time"] < 0])
}


lineages_edges <- function(x) {
  UseMethod("lineages_edges", x)
}

lineages_edges.cs_lineages <- function(x){
  x <- lin_id_to_int(c(lineages_active(x), lineages_extinct(x)))
  return(int_to_lin_id(sort(x)))
}



lineages_tips_n <- function(x) {
  UseMethod("lineages_tips_n", x)
}

lineages_tips_n.cs_lineages <- function(x){
  return(length(lineages_edges(x)))
}



is_incipient <- function(lineages, lin_id = NULL){
  if(is.null(lin_id)) {
    return(lineages[ , "status"] == -1)
  }
  stopifnot(length(lin_id) < nrow(lineages))
  
  return(lineages[lin_id, "status"] == -1)
}


is_good <- function(lineages, lin_id = NULL){
  if(is.null(lin_id)) {
    return(lineages[ , "status"] == 1)
  }
  stopifnot(length(lin_id) < nrow(lineages))
  
  return(lineages[lin_id, "status"] == 1)
}



is_extinct <- function(lineages, lin_id = NULL){
  if(is.null(lin_id)) {
    return(lineages[ , "status"] == -2)
  }
  stopifnot(length(lin_id) < nrow(lineages))
  
  return(lineages[lin_id, "status"] == -2)
}


lin_id_to_int <- function(lin_id) {
  as.integer(gsub("lin", "", lin_id))
}


int_to_lin_id <- function(i){
  paste0("lin", i)
}


traits_find_daughters <- function(traits, lin_id, format = c("int", "lin_id")) {
  format <- match.arg(format)
  
  daughters <- names(traits)[sapply(traits, 
                                    function(x) (x[["relative_id"]]) == lin_id_to_int(lin_id))]
  if (format == "int") {
    return(lin_id_to_int(daughters))}else{
      return(daughters)
    }
}


new_sim_result <- function(trees, t_end, step_size, pars, process_dead = FALSE, lineages = NULL, 
                           traits = NULL, full_results = FALSE, end_active_lineages = TRUE) {
  
  checkmate::assert_number(t_end)
  checkmate::assert_class(trees, "cs_sim_trees")
  checkmate::assert_class(pars, "cs_pars")
  checkmate::assert_class(lineages, "cs_lineages", null.ok = TRUE)
  checkmate::assert_true(is.null(traits) == is.null(lineages))
  checkmate::assert_logical(process_dead, len = 1)
  checkmate::assert_logical(full_results, len = 1)
  checkmate::assert_logical(end_active_lineages, len = 1)
  
  if(!is.null(traits)){
    checkmate::assert_true(length(traits) == nrow(lineages))
  }
  stopifnot(full_results & !is.null(lineages) & !is.null(traits))
  
  if(end_active_lineages & full_results){
    lineages[lineages[, "end_time"] < 0, "end_time"] <- t_end
  }
  
  x <- list(trees = trees, pars = pars, process_dead = process_dead, t_end = t_end, 
            step_size = step_size, lineages = lineages, traits = traits, 
            full_results = full_results)
  
  class(x) <- append(class(x),"cs_sim_results")
  
  return(x)
}

new_tree_results <- function(tree, tip_values, tree_gsp_fossil, tip_values_gsp_fossil,
                             tree_gsp_extant, tip_values_gsp_extant){
  
  checkmate::assert_class(tree, "phylo")
  checkmate::assert_class(tree_gsp_fossil, "phylo")
  checkmate::assert_class(tree_gsp_extant, "phylo", null.ok = TRUE)
  checkmate::assert_numeric(tip_values)
  checkmate::assert_numeric(tip_values_gsp_fossil)
  checkmate::assert_numeric(tip_values_gsp_extant, null.ok = TRUE)
  
  x <- list(all = list(tree = tree, tips = tip_values), 
            gsp_fossil = list(tree = tree_gsp_fossil, tips = tip_values_gsp_fossil), 
            gsp_extant = list(tree = tree_gsp_extant, tips = tip_values_gsp_extant))
  
  class(x) <- append(class(x), "cs_sim_trees")
  
  return(x)
}


