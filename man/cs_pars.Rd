% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pars.R
\name{cs_pars}
\alias{cs_pars}
\title{CAMSITE simulation parameters}
\usage{
cs_pars(
  lambda1 = 0.25,
  tau0 = 0.01,
  beta = 0.6,
  mu0 = 0.2,
  mu1 = 0.2,
  mubg = 0.01,
  mui0 = 0.4,
  mui1 = 0.4,
  muibg = 0.02,
  alpha1 = 0,
  alpha2 = 0,
  alpha3 = 0,
  sig2 = 0.5,
  m = 20
)
}
\arguments{
\item{lambda1}{speciation initiation rate}

\item{tau0}{basal speciation completion rate}

\item{beta}{effect of trait differences on the speciation completion rate}

\item{mu0}{competitive extinction parameter for good species}

\item{mu1}{selective extinction parameter for good species}

\item{mubg}{background good species extinction rate}

\item{mui0}{competitive extinction parameter for incipient species}

\item{mui1}{selective extinction parameter for incipient species}

\item{muibg}{background incipient species extinction rate}

\item{alpha1}{competition effect on extinction (competition strength)}

\item{alpha2}{competition effect on trait evolution (competition strength)}

\item{alpha3}{selection effect on extinction (selection strength)}

\item{sig2}{variance (rate) of Brownian motion (BM)}

\item{m}{relative contribution of character displacement (competition) with respect to stochastic (brownian) evolution}
}
\value{
object of class \code{pars} of parameter values
}
\description{
CAMSITE simulation parameters
}
\examples{
cs_pars()
}
