#' horseshoe_bfda_fc
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


fusion_bfda <- function(W, Z, E, model="Indep", prior="HSt", group=NULL, niter=6000, nwarmup=niter/2, nchains=3, thin=1, inits=NULL){
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
  
  Z = reshape2::melt(Z)
  Z$subj = Z$Var1
  Z$week = Z$Var2
  Z$Var1 = Z$Var2 = NULL
  Z = Z[complete.cases(Z),]
  
  dat = list(W = W$value,
             subj = W$subj,
             week = W$week,
             Z = Z$value,
             subj2 = Z$subj,
             week2 = Z$week,
             E = E,
             dim_space = dim(E)[2],
             N_weeks = length(unique(W$week)),
             N_obs = nrow(W),
             N_obs2 = nrow(Z),
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
    model.name <- "IndepFC.stan"
  }else if(model=="Group Indep"){
    if(prior=="HS"){
      model.name <- "GroupIndepHSFC.stan"
    }else if(prior=="HSt"){
      model.name <- "GroupIndepHStFC.stan"
    }else if(prior=="FHS"){
      model.name <- "GroupIndepFHSFC.stan"
    }else if(prior=="LASSO"){
      model.name <- "GroupIndepLASSOFC.stan"
    }
  }else if(model=="Group Hier"){
    if(prior=="HS"){
      model.name <- "GroupHierHSFC.stan"
    }else if(prior=="HSt"){
      model.name <- "GroupHierHStFC.stan"
    }else if(prior=="FHS"){
      model.name <- "GroupHierFHSFC.stan"
    }else if(prior=="LASSO"){
      model.name <- "GroupHierLASSOFC.stan"
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