#' Logit function
logit <- function(p){
  log(p) - log1p(-p)
}

#' Logistic function
logistic <- function(x){
  1 / (1 + exp(-x))
}

#' Add discrete stochastic random variables of arbitrary dimensions with overdispersion
#' 
#' @param dat A data.table with columns 'x', 'y', and 'probabilities'. Column 'x' is a numeric indicator unique for each random variable to be convoluted together. Column 'y' is the outcome value, and colum 'probabilities' contains the probability corresponding to each value of the 'y' column.
#' @param threshold An optional number to set a threshold. If set, the probabilities of all possible outcomes of the sum of the random variables are only caluclated below this value.
#' @param verbose A logical indicator which defaults to FALSE. Set to TRUE to print information.
#' @param overdispersion The standard deviation of a random effect in logit space. The random effect is not added to the first column.
#' @return the convoluted distribution.
#' @examples
#' convolution(dat = data.table(x = c(rep(1, 2), rep(2, 3), rep(3, 4)), y = c(1:2, 1:3, 1:4), 
#'     probability = c(3/4, 1/4, rep(1/3, 3), (1:4)/10)))
#' convolution(dat = data.table(x = c(rep(1, 2), rep(2, 3), rep(3, 4)), y = c(1:2, 1:3, 1:4), 
#'     probability = c(3/4, 1/4, rep(1/3, 3), (1:4)/10)), threshold = 6)
#' @export
ncdDetectOverdispersionSaddlepoint <- function(predictions, scores, overdispersion = 0, N = 10, method = c("naive", "numeric", "trapez"), observed_score) {
  # Make initial checks on input --------------------------------------------


# - data dimenstions
if (!identical(nrow(predictions), nrow(scores))) {
  stop("the row numbers of predictions and scores must be identical")
}

# - negative scores (cannot be handled)
if (any(scores < 0)) {
  stop("ncdDetect can currently not handle negative scores")
}

# row sums of predictions matrix
if (unique(round(rowSums(predictions), 10)) != 1) {
  stop("the row sums of the prediction matrix must sum to one")
}



# Convert data into data.tables -------------------------------------------

predictions <- as.data.table(predictions)
scores <- data.table(scores)



# Add x-values to data ----------------------------------------------------  
# (needed in the convolution step in order to know which positions go together)

scores[, x := 1:.N]



# Get scores in long format -------------------------------------------------

scores_long <- melt(scores, id.vars = c("x"), measure.vars = setdiff(names(scores), "x"),
                    variable.name = "mutation_type", value.name = "y", variable.factor = F)


# - set key columns for merging
setkeyv(scores_long, c("x", "mutation_type"))



# Integrate out random effect ---------------------------------------------


###
# Integrate out random effect
###

# Pre-computation for methods
if (method == "numeric"){
  val <- seq(-3*overdispersion, 3*overdispersion, length.out = N)
  m <- diff(pnorm(c(-Inf,val,Inf), sd = overdispersion)) # N+1
  weight <- (m[1:N]+m[2:(N+1)])/2
  weight[1] <- weight[1] + m[1]/2 # Make constant approximation for remainder of integrand
  weight[N] <- weight[N] + m[N+1]/2 # Same
}
if (method == "trapez"){
  # Transform density to bounded interval (-pi/2; pi/2)
  y <- seq(-pi/2,pi/2,length.out = N+2)[2:(N+1)]
  val <- tan(y) / overdispersion / 15 # x
  weight <- pi/(N+1)/2*rep(2,N) # Trapezoidal weight
  weight <- weight * (1+tan(y)^2)/overdispersion/15 * dnorm(val, 0, overdispersion)
}
if (method == "naive"){
  val <- rnorm(N, 0, sd = overdispersion)
  weight <- rep(1/N, N)
}

p_values <- vector("list", N)

for(i in 1:N){
  # - print progress
  #cat("calculating p-values using overdispersion. Iteration", i, "of", N, ".\n")
  
  random_effect <- val[i] # Specific method determines which random_effects to evaluate and their weights
  
  # - add random effect to predictions in logit space
  pred_cols <- names(predictions)[-1]
  pred_od <- data.table(logistic(logit(as.matrix(predictions[, pred_cols, with = F])) + random_effect))
  
  # - add column with probability of no mutation    
  no_mut_prob <- data.table(NO = 1-rowSums(pred_od))
  setnames(no_mut_prob, names(predictions)[1])
  
  pred_od <- cbind(no_mut_prob, pred_od)
  
  
  # - add x-values to predictions
  pred_od[, x := 1:.N]
  
  # - get predictions in long format
  predictions_long <- melt(pred_od, id.vars = c("x"), measure.vars = setdiff(names(pred_od), "x"),
                           variable.name = "mutation_type", value.name = "probability", variable.factor = F)
  
  
  # - set key columns for merging
  setkeyv(predictions_long, c("x", "mutation_type"))
  
  # - merge predictions and scores; set key columns
  dat <- predictions_long[scores_long]
  setkeyv(dat, c("x", "mutation_type"))
  
  # - perform significance evaluation using saddlepoint approximation
  p_values[[i]] <- saddlepoint(t = observed_score, dat, log = F) 

  
  # - clean up
  rm(random_effect, pred_cols, pred_od, no_mut_prob, predictions_long, dat)
  
}

out <- data.table(p = unlist(p_values), weight = weight)
final_p <- out[, sum(p * weight)]



# Return results ----------------------------------------------------------

return(p = final_p)


}



