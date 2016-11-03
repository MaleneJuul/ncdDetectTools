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
#' @export
convolution <- function(dat, threshold=NA, verbose=F) {
  
  
  # make sure that the x’s have no ‘jumps’ (1:N) ----------------------------
  tmp_x <- data.table(x = dat[, unique(x)])
  setkey(tmp_x,x)
  tmp_x[, x_new := 1:.N]
  setkey(dat, x)
  dat <- tmp_x[dat][, .(x_new, probability, y)]
  setnames(dat, c("x", "probability", "y"))
  setkey(dat, x)
  
  
  # set threshold -----------------------------------------------------------
  max_possible_score <- dat[, max(y), by = x][, sum(V1)]
  if (is.na(threshold)) {
    thres <- max_possible_score
  } else {
    thres <- threshold
  }
  
  
  # set dummy vector of integers --------------------------------------------
  ints <- 0:(thres+1)
  
  
  # perform convolution -----------------------------------------------------
  out <- as.data.table(convolution_body(as.matrix(dat[, list(y, probability, x)]), threshold=(thres+1), integers=ints))
  setnames(out, c("y", "probability"))
  
  
  
  # return results ----------------------------------------------------------
  return(out)
  
}
