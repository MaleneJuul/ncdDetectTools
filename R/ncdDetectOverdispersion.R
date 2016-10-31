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
ncdDetectOverdispersion <- function(predictions, scores, overdispersion = 0, N = 10, method = c("naive", "is", "numeric"), observations = NA, thres = NA) {
  # Make initial checks on input --------------------------------------------
  
  observation_available <- F
  if (!is.na(observations[1])) {
    observation_available <- T
  }
  
  
  # - data dimenstions
  if (!identical(nrow(predictions), nrow(scores))) {
    stop("the row numbers of predictions, scores and observations must be identical")
  }
  
  if (observation_available) {
    if (!identical(nrow(predictions), nrow(observations))) {
      stop("the row numbers of predictions, scores and observations must be identical") 
    }
  }
  
  # - negative scores (cannot be handled)
  if (any(scores < 0)) {
    stop("ncdDetect can currently not handle negative scores")
  }
  
  # row sums of predictions matrix
  if (unique(rowSums(predictions)) != 1) {
    stop("the row sums of the prediction matrix must sum to one")
  }
  
  # row sums of observations must sum to one (must contain only one 1; the rest must be zero)
  if (observation_available) {
    if (unique(rowSums(observations)) != 1) {
      stop("each row in the observation matrix must contain one 1, with the rest being zeros")
    }
  }
  
  
  
  
  # Convert data into data.tables -------------------------------------------
  
  predictions <- as.data.table(predictions)
  scores <- data.table(scores)
  if (observation_available) {
    observations <- data.table(observations)
  }
  
  # Add x-values to data ----------------------------------------------------  
  # (needed in the convolution step in order to know which positions go together)
  
  predictions[, x := 1:.N]
  scores[, x := 1:.N]
  if (observation_available) {
    observations[, x := 1:.N]
  }
  
  
  # Get data in long format -------------------------------------------------
  

  scores_long <- melt(scores, id.vars = c("x"), measure.vars = setdiff(names(scores), "x"),
                      variable.name = "mutation_type", value.name = "y", variable.factor = F)
  
  if (observation_available) {
    observations_long <- melt(observations, id.vars = c("x"), measure.vars = setdiff(names(observations), "x"),
                              variable.name = "mutation_type", value.name = "observation", variable.factor = F)
  }
  
  
  # Combine predictions, scores and observations ----------------------------
  
  # - set key columns for merging
  setkeyv(scores_long, c("x", "mutation_type"))
  
  
  
  
  ###
  # Integrate out random effect
  ###
  conv_results <- list()
  for(i in 1:N){
    random_effect <- rnorm(1, 0, sd = overdispersion)
    
    pred_cols <- setdiff(names(predictions), "x")[-1]
    pred_cols <- setdiff(names(predictions), "x")
    
    pred_od <- data.table(logistic(logit(as.matrix(predictions[, pred_cols, with = F])) + random_effect))
    
    ### malene fikser
    pred_od[, V1 := 1 - (V2 + V3)]
    
    pred_od[, x := 1:.N]
    
    
    predictions_long <- melt(pred_od, id.vars = c("x"), measure.vars = setdiff(names(pred_od), "x"),
                             variable.name = "mutation_type", value.name = "probability", variable.factor = F)
    
    
    # - set key columns for merging
    setkeyv(predictions_long, c("x", "mutation_type"))
    
    # - merge predictions and scores; set key columns
    dat <- predictions_long[scores_long]
    setkeyv(dat, c("x", "mutation_type"))
    
    # - finally, merge observations onto the data
    if (observation_available) {
      setkeyv(observations_long, c("x", "mutation_type"))
      dat <- dat[observations_long]
    }
    
   
    # Perform convolution ----------------------------------------------------- 
    conv_results[[i]] <- convolution(dat, threshold = thres)
    
  }
  
  # Average over results
  conv_results_df <- suppressWarnings(Reduce(function(x,y){merge(as.data.frame(x),as.data.frame(y), all = T, by = "y")}, conv_results))
  conv_results_df[is.na(conv_results_df)] <- 0
  score_dist <- rowMeans(conv_results_df[,-1])
  
  # Get observed score and p-value ------------------------------------------
  
  if (observation_available) {
    
    obs_score <- dat[observation == 1, sum(y)]
    
    p_value <- mean(sapply(conv_results, FUN = function(score_dist_df){score_dist_df[y >= 10, sum(probability)]}))
    
  }
  
  
  
  # Return results ----------------------------------------------------------
  
  if (observation_available) {
    return(list("score_dist" = score_dist, "obs_score" = obs_score, "p_value" = p_value))  
  } else {
    return(list("score_dist" = score_dist))  
  }
}