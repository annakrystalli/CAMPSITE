#' CAMSITE simulation parameters
#'
#' @param lambda1 speciation initiation rate
#' @param tau0 basal speciation completion rate
#' @param beta effect of trait differences on the speciation completion rate
#' @param mu0 competitive extinction parameter for good species
#' @param mu1 selective extinction parameter for good species
#' @param mubg background good species extinction rate
#' @param mui0 competitive extinction parameter for incipient species
#' @param mui1 selective extinction parameter for incipient species
#' @param muibg background incipient species extinction rate
#' @param alpha1 competition effect on extinction (competition strength)
#' @param alpha2 competition effect on trait evolution (competition strength)
#' @param alpha3 selection effect on extinction (selection strength)
#' @param sig2 variance (rate) of Brownian motion (BM)
#' @param m relative contribution of character displacement (competition) with respect to stochastic (brownian) evolution
#' @param s variance of BM per step size
#' @param step_size size of step of BM
#' @param age.max maximum age
#'
#' @return object of class `pars` of parameter values
#' @export
#'
#' @examples
#' cs_pars()
cs_pars <- function(lambda1 = 0.25, tau0 = 0.01, beta = 0.6, mu0 = 0.2, mu1 = 0.2,
                    mubg = 0.01, mui0 = 0.4, mui1 = 0.4, muibg = 0.02, 
                    alpha1 = 0, alpha2 = 0, alpha3 = 0, sig2 = 0.5, m = 20) {
  
  
  validate_pars(new_pars(lambda1 = lambda1, tau0, beta, mu0, mu1,
           mubg, mui0, mui1, muibg, 
           alpha1, alpha2, alpha3, sig2, m))
  
}

new_pars <- function(lambda1 = numeric(), tau0 = numeric(), beta = numeric(), mu0 = numeric(), 
                     mu1 = numeric(), mubg = numeric(), mui0 = numeric(), mui1 = numeric(),
                     muibg = numeric(), alpha1 = numeric(), alpha2 = numeric(), 
                     alpha3 = numeric(), sig2 = numeric(), m = numeric()) {
  
  structure(as.list(environment()), class = "cs_pars")
  
  }

validate_pars <- function(pars) {
  
  if(any(!sapply(pars, is.numeric))){
    usethis::ui_stop("non-numeric values supplied for parameters {usethis::ui_field(names(pars[!sapply(pars, is.numeric))]))}")
  }
  
  if(any(is.null(pars))){
    usethis::ui_stop("{usethis::ui_value('NULL')} values supplied for parameters {usethis::ui_field(names(pars[is.null(pars)]))}")
  }
  
  if(any(is.na(pars))){
    usethis::ui_stop("{usethis::ui_value(NA)} supplied for parameters {usethis::ui_field(names(pars[is.na(pars)]))}")
  }
  
  pars
}


