
validate_bounds <- function(bounds, root.value) {
  
  checkmate::assert_numeric(bounds, len = 2)

  if (root.value < bounds[1] || root.value > bounds[2]) {
    usethis::ui_warn("Root value outside boundaries; continuing simulation")
  }
  
}

validate_ou <- function(ou) {
  
  checkmate::assert_list(ou)
  checkmate::assert_subset(names(ou), choices = c("opt", "alpha4"))
  
  if(is.null(ou$opt)){
    usethis::ui_warn("no OU parameters supplied; continuing without OU process")
    ou$alpha4 <- 0
    return(ou)
  }
  
    if (is.null(ou$alpha4)){
      usethis::ui_stop("Alpha values for OU process must be provided")
    } 
  
  checkmate::assert_numeric(ou$opt)
  checkmate::assert_numeric(ou$alpha4)
  
  if (length(ou$alpha4) == length(ou$opt)) {
    return(ou)
  }
  
  if (length(ou$alpha4) > length(ou$opt)){
      usethis::ui_stop("Number of optima and alpha values must match")
    } else if (length(ou$alpha4) == 1) {
      ou$alpha4 = rep(ou$alpha4, length(ou$opt))
    }
  return(ou)
}


