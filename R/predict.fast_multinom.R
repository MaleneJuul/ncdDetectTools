predict.fast_multinom <- function (object, newdata = NULL) {
  if (is.null(newdata)) {
    y <- object$predictions
    if (is.null(y)) 
      warning("The predictions on the original data couldn't be retreived, because they weren't saved along with the fit (predictions=F in the function call fast_multinom).")
  }
  else {
    if (is.data.table(newdata)) {
      warning("The class of the new data is changed from data.table to data.frame.")
      data = as.data.frame(data)
    }
    M = length(object$formulae)
    n = nrow(newdata)
    lin.pred = matrix(NA, ncol = M, nrow = n)
    left.side = as.character(attr(terms(object$formulae[[1]]), 
                                  "variables")[2])
    refname.regexp = regexpr("(?<=, )\\w+(?=\\)$)", left.side, 
                             perl = T)
    refname = regmatches(left.side, refname.regexp)
    if (length(regmatches) == 0) {
      stop("Can't find out the name of the reference category")
    }
    outputnames = c(refname, names(object$coefficients))
    for (name in outputnames) {
      if (!name %in% colnames(newdata)) {
        currnames = names(newdata)
        newdata = cbind(newdata, numeric(n))
        names(newdata) = c(currnames, name)
      }
    }
    for (m in 1:M) {
      mm = model.Matrix(object$formulae[[m]], data = newdata, 
                        xlev = object$xlevels[[m]], contrasts.arg = object$contrasts[[m]], 
                        sparse = T)
      lin.pred[, m] = as.numeric(mm %*% object$coefficients[[m]])
    }
    logits = exp(lin.pred)
    phi0 = 1/(1 + rowSums(logits))
    y = cbind(phi0, logits * phi0)
    colnames(y) = c(all.vars(object$formulae[[1]])[2], names(object$coefficients))
  }
  return(y)
}
