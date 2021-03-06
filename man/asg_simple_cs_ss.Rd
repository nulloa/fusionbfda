% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asg_simple_cs_ss.R
\name{asg_simple_cs_ss}
\alias{asg_simple_cs_ss}
\title{asg_simple_cs_ss}
\usage{
asg_simple_cs_ss(y, x, group, season, priors, niter = 4000,
  nwarmup = niter/2, nchains = 3, thin = 1, inits = NULL)
}
\arguments{
\item{y}{response variable which follows binomial dist}

\item{x}{explanatory variable}

\item{group}{group of response}

\item{season}{season of response}

\item{priors}{list of priors}

\item{niter}{number of interations to be run (default=4000)}

\item{nchains}{number of chains to be run (default=3)}

\item{thin}{when you want to thin (default=1)}

\item{inits}{Add specific initial values}

\item{count}{n in binomial dist}

\item{nwarmiup}{number of iterations to be used as warmup (see link below)}
}
\value{
A MCMC object
}
\description{
Runs Asymmetric Gaussian MCMC with a region mean structure accross the seasons
}
\examples{
priors <- list()
priors$m0
priors$C0

}
\seealso{
\url{http://mc-stan.org/users/documentation/}
}
