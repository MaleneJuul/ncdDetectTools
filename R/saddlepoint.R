#' Saddlepoint Approximation
#'
#' Let \eqn{X_1,X_2,\ldots X_N} be independent Bernouilli trials
#' with \eqn{X_i \sim p_i}.
#' We evaluate the probability
#' \deqn{P(s_1X_1+s_2X_2\ldots s_NX_N > t)}.
#' 
#' @param t   Point where tail probability is evaluated (note this function is not vectorized over t)
#' @param p   Probabilities of each Bernouilli trial
#' @param s   Score associated with each trial
#' @param lattice Lattice size in minimal lattice
#' @param log Return log probability
#' 
#' @examples
#' saddlepoint(300, rep(0.01, 10000))
#' pbinom(299, 10000, 0.01, lower.tail = F, log.p = T)
#' @importFrom Rcpp evalCpp
#' @export
saddlepoint <- function(t, dat, lattice = 1L, log = T){
  # Sort data frame after x
  dat <- dat[order(x),]
  
  x <- dat$x
  p <- dat$probability
  s <- dat$y
  
  stopifnot(all(s >= 0))
  
  # Return 0 and 1 if outside range
  
  # Giver saddelpunkts approksimationen paa log-skala
  # Use a combination of bisection and Newton raphson
  
  theta_old <- -Inf
  theta <- 0
  
  # 
  theta_largest_negative <- -Inf # Largest theta that gives a negative value so far
  theta_smallest_positive <- Inf # Smallest theta that gives a positive value so far
  
  iter <- 0
  while(abs(theta - theta_old) > 1e-8){
    cumDer <- cumulantDerivatives(theta, x, p, s)
    # For convergence properties
    iter <- iter + 1
    theta_old <- theta
    
    # Newton raphson
    val <- cumDer[2] - t
    if(val < 0 && theta > theta_largest_negative)
      theta_largest_negative <- theta
    if(val > 0 && theta < theta_smallest_positive)
      theta_smallest_positive <- theta
    deriv <- cumDer[3]
    step <- - val / deriv 
    
    # Next theta if newton raphson
    theta <- theta + step
    if(theta < theta_largest_negative || theta > theta_smallest_positive){
      # Fall back to bisection
      theta <- (theta_largest_negative+theta_smallest_positive)/2
    }
  }
  
  cumDer <- cumulantDerivatives(theta, x, p, s)
  
  # El approximazione
  v <- cumDer[3]
  
  # Uden lattice-korrektion
  if(lattice == 0){
    ret <- cumDer[1]-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = F, log.p = T)
  } else{
    ret <- cumDer[1]-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = F, log.p = T) +
            log(abs(theta*lattice))- log((1-exp(-lattice*abs(theta))))
  }
  
  ifelse(log, ret, exp(ret))
}
