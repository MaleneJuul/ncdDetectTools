
# Info --------------------------------------------------------------------

# Function to define scores
# Input is data set (dat) with columns
# - position ("end")
# - mutation probabilities [multiple columns] ("probability_TYPE_mut")
# - probabilities for no mutation ("probability_nomut")
# - sample identification ("sample_id")
# and scoring.scheme defined as either "nMut", "phyloP", or "logLikelihoods"

define_scores <- function(dat, scoring_scheme) {
  
  # initial checks on input
  if (!("probability_nomut" %in% names(dat))) stop("there needs to be a column named 'probability_nomut")
  
  score_cols <- grep(pattern = "probability", x = names(dat))
  if (unique(round(rowSums(dat[, score_cols, with = F]), 10)) != 1) stop("columns with predictions must sum to one")
  
  prob_cols <- setdiff(names(dat)[score_cols], "probability_nomut")
  prob_cols_namecheks <- c(sapply(strsplit(prob_cols, "_"), FUN = function(x) x[1] == "probability"),
                           sapply(strsplit(prob_cols, "_"), FUN = function(x) x[3] == "mut"))
  if (sum(prob_cols_namecheks * 1) != length(prob_cols_namecheks)) stop("name your columns with predictions 'probability_TYPE_mut'")
  
  rm(score_cols, prob_cols, prob_cols_namecheks)
  
  
  
  
  discretize_score <- function(values, granularity) {
    posneg <- sign(values)
    z = abs(values)*10^granularity
    z = trunc(z)
    z*posneg
  }
  
  if (scoring_scheme == "log_likelihoods") {
    
    # copy data set
    dat_copy <- copy(dat)
    
    # get columns for which to calculate a score
    score_cols <- grep(pattern = "probability", x = names(dat_copy))
    
    for (i in 1:length(score_cols)) {
      
      # get column with positions
      pos_col <- which(names(dat_copy) == "end")
      
      # get column with sample_id
      sample_col <- which(names(dat_copy) == "sample_id")
      
      # take out position and sample columns as well as the i'th column with probabilities
      tmp <- dat_copy[, c(pos_col, sample_col, score_cols[i]), with = F]
      
      # change names of new tmp data set so it's easy to work with  
      setnames(tmp, c("end", "sample_id", "probability"))     
      
      # calculate score
      tmp[, score := discretize_score(-log(probability), 1)]
      tmp <- tmp[, .(end, sample_id, score)]
      
      # set new column id
      new_id <- gsub(pattern = "probability", replacement = "score", x = names(dat_copy)[score_cols[i]])
      setnames(tmp, c("end", "sample_id", new_id))
      
      # add information to original data set
      setkeyv(tmp, c("end", "sample_id"))
      setkeyv(dat, c("end", "sample_id"))
      
      dat <- tmp[dat]
    }  
    
    
  } else if (scoring_scheme == "n_mutations") {
    
    # copy data set
    dat_copy <- copy(dat)
    
    # get columns for which to calculate a score
    score_cols <- grep(pattern = "probability", x = names(dat_copy))
    
    for (i in 1:length(score_cols)) {
      
      # get column with positions
      pos_col <- which(names(dat_copy) == "end")
      
      # get column with sample_id
      sample_col <- which(names(dat_copy) == "sample_id")
      
      # take out position and sample columns as well as the i'th column with probabilities
      tmp <- dat_copy[, c(pos_col, sample_col, score_cols[i]), with = F]
      
      # change names of new tmp data set so it's easy to work with  
      setnames(tmp, c("end", "sample_id", "probability"))     
      
      # calculate score
      if (names(dat_copy)[score_cols[i]] == "probability_nomut") {
        tmp[, score := 0]
        tmp <- tmp[, .(end, sample_id, score)]
      } else {
        tmp[, score := 1]
        tmp <- tmp[, .(end, sample_id, score)]
      }
      
      # set new column id
      new_id <- gsub(pattern = "probability", replacement = "score", x = names(dat_copy)[score_cols[i]])
      setnames(tmp, c("end", "sample_id", new_id))
      
      # add information to original data set
      setkeyv(tmp, c("end", "sample_id"))
      setkeyv(dat, c("end", "sample_id"))
      
      dat <- tmp[dat]
    }  
    
    
  } else if (scoring_scheme == "phyloP") {
    
    # shift scores, so they live on the positive integer grid
    dat[, phyloP_shifted := phyloP + 20]
    
    # copy data set
    dat_copy <- copy(dat)
    
    # get columns for which to calculate a score
    score_cols <- grep(pattern = "probability", x = names(dat_copy))
    
    for (i in 1:length(score_cols)) {
      
      # get column with positions
      pos_col <- which(names(dat_copy) == "end")
      
      # get column with sample_id
      sample_col <- which(names(dat_copy) == "sample_id")
      
      # get column with phyloP score
      phyloP_col <- which(names(dat_copy) == "phyloP_shifted")
      
      # take out position and sample columns as well as the i'th column with probabilities
      tmp <- dat_copy[, c(pos_col, sample_col, phyloP_col), with = F]
      
      # change names of new tmp data set so it's easy to work with  
      setnames(tmp, c("end", "sample_id", "phyloP_shifted"))     
      
      # calculate score
      if (names(dat_copy)[score_cols[i]] == "probability_nomut") {
        tmp[, score := 0]
        tmp <- tmp[, .(end, sample_id, score)]
      } else {
        tmp[, score := discretize_score(phyloP_shifted, 1)]
        tmp <- tmp[, .(end, sample_id, score)]
      }
      
      # set new column id
      new_id <- gsub(pattern = "probability", replacement = "score", x = names(dat_copy)[score_cols[i]])
      setnames(tmp, c("end", "sample_id", new_id))
      
      # add information to original data set
      setkeyv(tmp, c("end", "sample_id"))
      setkeyv(dat, c("end", "sample_id"))
      
      dat <- tmp[dat]
    }  
    
    dat[, phyloP_shifted := NULL]
    
  } else {
    
    stop("Only scoring schemes 'phyloP', 'log_likelihoods', are 'n_mutations' is currently implemented.")  
    
  }
  
  setkey(dat, end)
  
  return(dat)
  
}


