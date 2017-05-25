#' Add discrete stochastic random variables of same dimension
#' 
#' @param predictions A matrix in which each row corresponds to a discrete random variable and each column corresponds to a discrete outcome. The matrix contains the probabilities of each outcome, and each row must sum to one.
#' @param scores A matrix in which each row corresponds to a discrete random variable and each column corresponds to a discrete outcome. The matrix contains the values for each outcome for each random variable.
#' @param observations A matrix (optional) in which each row corresponds to a discrete random variable and each column corresponds to a discrete outcome. The matrix contains the observed outcome for each variable. Each row must contain exactly one 1, while the rest of the entries must be 0.
#' @param thres An optional number to set a threshold. If set, the probabilities of all possible outcomes of the sum of the random variables are only caluclated below this value.
#' @return score_dist The convoluted distribution of the sum of the discrete random variables defined by input matrices predictions and scores.
#' @return obs_score The observed score. Only returned if matrix observations is provided as input.
#' @return p_value The p-value resulting from evaluating the observed score in the convoluted distribution. Only returned if matrix observations is provided as input.
#' @examples
#' ncdDetect(predictions = matrix(rep(1/6, 12), nrow = 2), scores = matrix(rep(1:6, 2), nrow = 2, byrow = TRUE), 
#'     observations = matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0), nrow = 2, byrow = TRUE))
#' ncdDetect(predictions = matrix(rep(1/6, 12), nrow = 2), scores = matrix(rep(1:6, 2), nrow = 2, byrow = TRUE), thres = 6)
#' ncdDetect(predictions = matrix(rep(1/6, 12), nrow = 2), scores = matrix(rep(1:6, 2), nrow = 2, byrow = TRUE))
#' @useDynLib ncdDetectTools
#' @import data.table
#' @export
ncdDetect <- function(predictions, scores, observations = NA, thres = NA) {
  

  # Make initial checks on input --------------------------------------------

  observation_available <- F
  if ( ! any(is.na(observations)) ) {
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
  if (unique(round(rowSums(predictions),10)) != 1) {
    stop("the row sums of the prediction matrix must sum to one")
  }

  # row sums of observations must sum to one (must contain only one 1; the rest must be zero)
  if (observation_available) {
    if (unique(round(rowSums(observations), 10)) != 1) {
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

  predictions_long <- melt(predictions, id.vars = c("x"), measure.vars = setdiff(names(predictions), "x"),
                           variable.name = "mutation_type", value.name = "probability", variable.factor = F)
  
  scores_long <- melt(scores, id.vars = c("x"), measure.vars = setdiff(names(scores), "x"),
                      variable.name = "mutation_type", value.name = "y", variable.factor = F)
  
  if (observation_available) {
    observations_long <- melt(observations, id.vars = c("x"), measure.vars = setdiff(names(observations), "x"),
                              variable.name = "mutation_type", value.name = "observation", variable.factor = F)
  }
  
  
  # Combine predictions, scores and observations ----------------------------
  
  # - set key columns for merging
  setkeyv(predictions_long, c("x", "mutation_type"))
  setkeyv(scores_long, c("x", "mutation_type"))
 
  # - merge predictions and scores; set key columns
  dat <- predictions_long[scores_long]
  setkeyv(dat, c("x", "mutation_type"))
  
  # - finally, merge observations onto the data
  if (observation_available) {
    setkeyv(observations_long, c("x", "mutation_type"))
    dat <- dat[observations_long]
  } 
  

  # Perform convolution ----------------------------------------------------- 

  score_dist <- convolution(dat, threshold = thres)
    

  # Get observed score and p-value ------------------------------------------
  
  if (observation_available) {
    
    obs_score <- dat[observation == 1, sum(y)]
    
    p_value <- score_dist[y >= obs_score, sum(probability)]
    
  }
  

  
  # Return results ----------------------------------------------------------
  
  if (observation_available) {
    return(list("score_dist" = score_dist, "obs_score" = obs_score, "p_value" = p_value))  
  } else {
    return(list("score_dist" = score_dist))  
  }
  
}












