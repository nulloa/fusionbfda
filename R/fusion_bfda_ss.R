#' fusion_bfda_ss
#' 
#' Runs BFDA with horseshoe prior on betas for forecasting
#'
#' @param W matrix of response
#' @param E Eigen functions that correspond to W
#' @param model Which hierarchy structure to be used: 'Indep' (default), 'Group Indep', or 'Group Hier'
#' @param prior Which prior to use to specify model: 'HS'(Horseshoe), 'HSt'(Horseshoe with t), 'LASSO', 'FHS'(Finnish HS)
#' @param group Grouping variable to be included if model is Region or Season 
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
#'
#' @export


fusion_bfda_ss <- function(W, E, model="Indep", prior="HSt", group=NULL, niter=6000, nwarmup=niter/2, nchains=3, thin=1, inits=NULL){
  # Load Library
  require(rstan)
  
  if(nchains>1){
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
  }
  
  # Setup data for model
  W = reshape2::melt(W)
  W$subj = W$Var1
  W$week = W$Var2
  W$Var1 = W$Var2 = NULL
  W = W[complete.cases(W),]
  
  dat = list(W = W$value,
             subj = W$subj,
             week = W$week,
             E = E,
             dim_space = dim(E)[2],
             N_weeks = length(unique(W$week)),
             N_obs = nrow(W),
             N_subj = length(unique(W$subj))
  )
  
  
  # Include Group variable if necessary
  if(!(model %in% c("Indep", "Group Indep"))){
    
    group_list = list(group=group,
                      N_group=length(unique(group)))
    dat <- c(dat, group_list)
  }
  
  
  # Set up the model in stan
  model.name <- NULL
  if(model=="Indep"){
    model.name <- "IndepFC_ss.stan"
  }else if(model=="Group Indep"){
    if(prior=="HSt"){
      model.name <- "GroupIndepHStFC_ss.stan"
    }else{
      cat("Misspecified prior.")
      stop()
    }
  }else if(model=="Group Hier"){
    if(prior=="HSt"){
      model.name <- "GroupHierHStFC_ss.stan"
    }else{
      cat("Misspecified prior.")
      stop()
    }
  }else{
    cat("Misspecified model. Model needs to be 'Indep', 'Group Indep', or 'Group Hier'")
    stop()
  }
  
  m <- stan(file = system.file("model", model.name, package = "fusionbfda"), 
            data = dat, iter = niter, warmup=nwarmup, thin=thin, chains = nchains,# init = inits, 
            control = list(adapt_delta = 0.99, max_treedepth = 15))
  return(m)
}