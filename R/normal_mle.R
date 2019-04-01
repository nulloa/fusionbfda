#' normal_mle
#' 
#' A function which finds and formats the mle of ASG ff using normal data model
#' 
#' @param y response variable which follows binomial dist
#' @param x explanatory variable
#' @param group groups of response
#' @param season seasons in response
#' @param inits initial values for optim (defualt=NULL)
#' @param ll likelihood to be optimized (defualt=ASG)
#'
#' @return A function
#'
#' @examples
#' asg_mle(y, x, group, seas)
#'
#' @export


normal_mle <- function(y, x, group, season, inits=NULL, ll=NULL, correction=TRUE){
  
  log.lik <- function(y, x, par){
    theta <- boot::inv.logit(asg(x, par[1], par[2], par[3], par[4], par[5], par[6]))
    # ll <- sum(lchoose(num, y)) + sum(y*log(theta)) + sum((num-y)*log(1-theta))
    n <- length(y)
    ll <- -(n/2)*log(2*pi) - (n/2)*log(par[7]^2) - (1/(2*par[7]^2))*sum((y - theta)^2)
    if(par[3] <= 0 | par[5] <= 0 | par[6] <= 0 | par[7] <= 0){ll <- -Inf}
    return(-ll)
  }
  
  if(is.null(ll)){ll <- log.lik}
  if(is.null(inits)){inits <- c(-4.5, -4.5, 15, -2.5, 4, 5, 0.1)}
  
  df <- data.frame(y=y, x=x, group=group, season=season)
  opt <- array(NA, dim=c(length(unique(df$group)), length(unique(df$season)), 6))
  
  for(s in 1:length(unique(df$season))){
    for(g in 1:length(unique(df$group))){
      dat <- subset(df, season==paste0(unique(df$season)[s]) & group==paste0(unique(df$group)[g]))
      opts <- optim(par=inits, ll, y=dat$y, x=dat$x)
      newopts <-  c(opts$par[1], opts$par[2], opts$par[4], opts$par[3], opts$par[5], opts$par[6])
      opt[unique(df$group)[g], unique(df$season)[s], ] <- newopts
    }
  }
  
  
  if(correction==TRUE){
    s <- length(unique(df$season))
    for(g in 1:length(unique(df$group))){
      whichoff <- abs(opt[g, s, ] - opt[g, s-1, ]) > 3
      opt[g, s, ][whichoff] <- opt[g, s-1, ][whichoff]
    }
  }
  
  return(opt)
}





#' asg
#' 
#' A function for the ff of the ASG Distribution
#' 
#' @param x explanatory variable
#' @param beta1 intercept in first half
#' @param beta2 intercept in second half
#' @param mu peak week
#' @param h peak
#' @param sigma1 variance in first half
#' @param sigma2 variance in second half
#'
#' @return A function
#'
#' @examples
#' asg(x=c(1:33), beta1=-4, beta2=-1, mu=15, h=10, sigma1=10, sigma2=15)
#'
#' @export



asg <- Vectorize(function(x, beta1, beta2, mu, h, sigma1, sigma2){
  top <- beta1 + (h - beta1)*exp(-((x - mu)^2)/(2*sigma1^2))
  bot <- beta2 + (h - beta2)*exp(-((x - mu)^2)/(2*sigma2^2))
  ifelse(x < mu,return(top),return(bot))
})
