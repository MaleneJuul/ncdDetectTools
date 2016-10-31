library(ncdDetectTools)
context("ncdDetectOverdispersion")

test_that("Betabinomial (overdispersion) distribution 2 columns",{
  predictions <- matrix(c(0.999, 0.001), 1000, 2, byrow = T)
  scores <- matrix(c(0,1), 1000, 2, byrow = T)
  observed_score <- 10
  #observations <- matrix(c(1,0), 1000, 2, byrow = T)
  #observations[1:10,2] <- 1; observations[1:10, 1] <- 0
  overdispersion <- 0.3
  
  
  res <- ncdDetectOverdispersion(predictions, scores, overdispersion, N = 100, method = "naive", observed_score = observed_score)[type == "naive", p_value]
  bbinom_p_value <- VGAM::dbetabinom(10, 1000, 0.001, overdispersion^2*0.001/(1-0.001), log = T)

  expect_lt(abs(log(res)-bbinom_p_value), log(3)) # Within a factor 2 using the naive approach
})


test_that("Betabinomial (overdispersion) distribution 3 columns",{
  predictions <- matrix(c(0.999, 0.0005, 0.0005), 1000, 3, byrow = T)
  scores <- matrix(c(0,1,1), 1000, 3, byrow = T)
#   observations <- matrix(c(1,0,0), 1000, 3, byrow = T)
#   observations[1:10,3] <- 1; observations[1:10, 1] <- 0
  observed_score <- 10
  overdispersion <- 0.3
  
  res <- ncdDetectOverdispersion(predictions, scores, overdispersion, N = 100, method = "naive", observed_score = observed_score)[type == "naive", p_value]
  bbinom_p_value <- VGAM::dbetabinom(10, 1000, 0.001, overdispersion^2*0.001/(1-0.001), log = T)
  
  expect_lt(abs(log(res)-bbinom_p_value), log(3)) # Within a factor 2 using the naive approach
})
