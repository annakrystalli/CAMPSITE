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
#' @export
#' @examples
#' new_trait(1, 1, 2, c(NA, root.value))
new_trait <- function(lin_id = integer(1), status = integer(1), relative_id = integer(1), 
                      parent_trait_value = numeric(1), trait_values = numeric()) {
  checkmate::assert_numeric(lin_id, len = 1)
  checkmate::assert_numeric(status, len = 1)
  checkmate::assert_numeric(relative_id, len = 1)
  checkmate::assert_numeric(parent_trait_value, len = 1)
  checkmate::assert_numeric(trait_values)
  
  x <- list(lin_id, status, relative_id, parent_trait_value, trait_values)
  names(x) <- c("lin_id", "status", "relative_id", "parent_trait_value", "trait_values")
  structure(x, class = "cs_trait")

  
  
}

#' @export
current_trait_value <- function(x) {
  UseMethod("current_trait_value", x)
}
#' @export
current_trait_value.cs_trait <- function(x){
  return(tail(x$trait_values, 1))
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
#' @export
#' @examples
#' new_lineage(c(1, 0, 0, -1, 1, 0, 0))
new_lineage <- function(x = numeric(length = 7L)) {
  checkmate::assert_numeric(x, len = 7)
  names(x) <- c("parent_node", "descendant_node", "start_time", "end_time", 
                      "status", "spec_or_ext_ct", "spec_ct")
  structure(x, class = "cs_lineage")
}


#' @export
new_lineages <- function(...) {
  
  lapply(list(...), 
         FUN = function(obj){checkmate::assert_class(obj, classes = "cs_lineage")})
  checkmate::assert_true(...length() > 0)

  x <- rbind(...)
  rownames(x) <- int_to_lin_id(1:...length())
  
  structure(x, class = "cs_lineages")
}


#' @export
lineage_next_parent <- function(x) {
  UseMethod("lineage_next_parent", x)
}
#' @export
lineage_next_parent.cs_lineages <- function(x){
  return(max(x[, "parent_node"]) + 1)
}


#' @export
lineage_next <- function(x) {
  UseMethod("lineage_next", x)
}
#' @export
lineage_next.cs_lineages <- function(x){
  return(max(lin_id_to_int(rownames(x))) + 1)
}

#' @export
lineage_add <- function(x, ...) {
  UseMethod("lineage_add", x)
}
#' @export
lineage_add.cs_lineages <- function(x, ...){
  checkmate::assert_class(x, classes = "cs_lineages")
  
  lapply(list(...), 
         FUN = function(obj){checkmate::assert_class(obj, classes = "cs_lineage")})
  
  lin_names <- rownames(x)
  lin_max <- max(lin_id_to_int(lin_names))
  lin_names <- c(lin_names, int_to_lin_id(1:...length() + lin_max))
  
  
  y <- rbind(x, ...)
  rownames(y) <- lin_names
  
  structure(y, class = "cs_lineages")
  
}


#' @export
lineage_last_speciated <- function(x) {
  UseMethod("lineage_last_speciated", x)
}
#' @export
lineage_last_speciated.cs_lineages <- function(x){
  return(tail(rownames(x), 2))
}

#' @export
lineage_last_incipient <- function(x) {
  UseMethod("lineage_last_incipient", x)
}
#' @export
lineage_last_incipient.cs_lineages <- function(x){
  x <- x[lineage_last_speciated(x), ]
  return(rownames(x)[is_incipient(x)])
}

#' @export
lineage_last_parent <- function(x) {
  UseMethod("lineage_last_parent", x)
}
#' @export
lineage_last_parent.cs_lineages <- function(x){
  x <- x[lineage_last_speciated(x), ]
  return(rownames(x)[!is_incipient(x)])
}

#' @export
lineages_dead <- function(x) {
  UseMethod("lineages_dead", x)
}
#' @export
lineages_dead.cs_lineages <- function(x){
  return(rownames(x)[x[ , "end_time"] > 0])
}


#' @export
lineages_born <- function(x) {
  UseMethod("lineages_born", x)
}
#' @export
lineages_born.cs_lineages <- function(x){
  return(rownames(x)[x[ , "start_time"] > 0])
}


#' @export
lineages_active <- function(x) {
  UseMethod("lineages_active", x)
}
#' @export
lineages_active.cs_lineages <- function(x){
  return(rownames(x)[x[ , "status"] > -2])
}


#' @export
is_incipient <- function(lineages, lin_id = NULL){
  if(is.null(lin_id)) {
    return(lineages[ , "status"] == -1)
  }
  stopifnot(length(lin_id) < nrow(lineages))
  stopifnot(is.numeric(lin_id))
  
  return(lineages[lin_id, "status"] == -1)
}

#' @export
is_good <- function(lineages, lin_id = NULL){
  if(is.null(lin_id)) {
    return(lineages[ , "status"] == 1)
  }
  stopifnot(length(lin_id) < nrow(lineages))
  stopifnot(is.numeric(lin_id))
  
  return(lineages[lin_id, "status"] == 1)
}

#' @export
lin_id_to_int <- function(lin_id) {
  as.integer(gsub("lin", "", lin_id))
}

#' @export
int_to_lin_id <- function(i){
  paste0("lin", i)
}

#' @export
traits_find_daughters <- function(traits, lin_id, format = c("int", "lin_id")) {
 format <- match.arg(format)
  
  daughters <- names(traits)[sapply(traits, 
                       function(x) (x[["relative_id"]]) == lin_id_to_int(lin_id))]
  if (format == "int") {
    return(lin_id_to_int(daughters))}else{
      return(daughters)
    }
}

