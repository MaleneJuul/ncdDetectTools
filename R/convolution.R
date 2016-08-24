#' Add discrete stochastic random variables of arbitrary dimensions
#' 
#' @param dat A data.table with columns 'x', 'y', and 'probabilities'. Column 'x' is a numeric indicator unique for each random variable to be convoluted together. Column 'y' is the outcome value, and colum 'probabilities' contains the probability corresponding to each value of the 'y' column.
#' @param threshold An optional number to set a threshold. If set, the probabilities of all possible outcomes of the sum of the random variables are only caluclated below this value.
#' @param verbose A logical indicator which defaults to FALSE. Set to TRUE to print information.
#' @return the convoluted distribution.
#' @examples
#' convolution(dat = data.table(x = c(rep(1, 2), rep(2, 3), rep(3, 4)), y = c(1:2, 1:3, 1:4), 
#'     probability = c(3/4, 1/4, rep(1/3, 3), (1:4)/10)))
#' convolution(dat = data.table(x = c(rep(1, 2), rep(2, 3), rep(3, 4)), y = c(1:2, 1:3, 1:4), 
#'     probability = c(3/4, 1/4, rep(1/3, 3), (1:4)/10)), threshold = 6)

convolution <- function(dat, threshold=NA, verbose=F) {
  
  require(data.table)
  
  # make sure that the x's has no "jumps" (1:dat[. .N])
  # -------------
  tmp_x <- data.table(x = dat[, unique(x)])
  setkey(tmp_x,x)
  tmp_x[, x_new := 1:.N]
  setkey(dat, x)
  dat <- tmp_x[dat][, .(x_new, probability, y)]
  setnames(dat, c("x", "probability", "y"))
  setkey(dat, x)

  
  # set threshold
  # -------------
  max_possible_score <- dat[, max(y), by = x][, sum(V1)]
  if (is.na(threshold)) {
    thres <- max_possible_score
  } else {
    thres <- threshold
  }
  
  
  # set dummy vector ints
  # -------------
  ints <- 0:thres
  
  
  # do calculations
  # -------------

  out <- as.data.table(convolution_body(as.matrix(dat[, list(y, probability, x)]), threshold=thres, integers=ints))
  setnames(out, c("y", "probability"))
  
  
  # if out contains no rows - i.e. if the threshold is set too low - 
  # recalculate with threshold set to max possible score
  # -------------
  if (out[, .N] == 0) {
    cat("the original convoluted distribution contains no rows! Calculating again without a threshold\n\n")
    thres <- max_possible_score
    ints <- 0:thres
    out <- as.data.table(convolution_body(as.matrix(dat[, list(y, probability, x)]), threshold=thres, integers=ints))
    setnames(out, c("y", "probability"))
  }
  
  
  # only add line in the end if a threshold is set
  if (!is.na(threshold)) {
    
    # add line at the bottom to capture the remaining probability mass (only if needed)
    # -------------
    if (out[, sum(probability)] < 1) {
      
      if (verbose) cat("an extra line is added to the distribution to catch the remaining probability mass")
      
      lastLine <- data.table(y = out[, max(y) + 1], probability = out[, 1-sum(probability)])
      out <- rbindlist(list(out, lastLine))
    }
  }

  # return results
  # -------------
  return(out)
  
}
