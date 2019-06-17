#' asg_simple_cs_ss
#' 
#' Runs Asymmetric Gaussian MCMC with a region mean structure accross the seasons
#'
#' @param y response variable which follows binomial dist
#' @param x explanatory variable
#' @param count n in binomial dist
#' @param group group of response
#' @param season season of response
#' @param priors list of priors
#' @param niter number of interations to be run (default=4000)
#' @param nchains number of chains to be run (default=3)
#' @param nwarmiup number of iterations to be used as warmup (see link below)
#' @param thin when you want to thin (default=1)
#' @param inits Add specific initial values
#' 
#' @seealso \url{http://mc-stan.org/users/documentation/}
#'
#' @return A MCMC object
#'
#' @examples
#' priors <- list()
#' priors$m0
#' priors$C0
#'
#' @export


asg_simple_cs_ss <- function(y, x, group, season, priors, niter=4000, nwarmup=niter/2, nchains=3, thin=1, inits=NULL){
  # Load Library
  require(rstan)
  
  if(nchains>1){
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
  }
  
  # Setup data for model
  dat = list(y = y, 
             x = x,
             group = group, 
             seas = as.numeric(season),
             k=6,
             nG = length(unique(group)), 
             n = length(y),
             nS = length(unique(season))
             )
  # Set priors
  dat <- c(dat, priors)
  
  # Set up the model in stan
  m <- stan(file = system.file("model", "asgSimpleCS_ss.stan", package = "fusionbfda"), 
            data = dat, iter = niter, warmup=nwarmup, thin=thin, chains = nchains, init = inits, 
            control = list(adapt_delta = 0.99, max_treedepth = 15))
  return(m)
}
