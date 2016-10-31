library(ncdDetectTools)
context("ncdDetectOverdispersion")

test_that("Betabinomial (overdispersion) distribution 2 columns",{
  predictions <- matrix(c(0.999, 0.001), 1000, 2, byrow = T)
  scores <- matrix(c(0,1), 1000, 2, byrow = T)
  observations <- matrix(c(1,0), 1000, 2, byrow = T)
  observations[1:10,2] <- 1; observations[1:10, 1] <- 0
  overdispersion <- 0.3
  
  res <- ncdDetectOverdispersion(predictions, scores, overdispersion, N = 100, method = "naive", observations = observations)
  bbinom_p_value <- VGAM::dbetabinom(10, 1000, 0.001, overdispersion^2*0.001/(1-0.001), log = T)
  
  expect_lt(abs(log(res$p_value)-bbinom_p_value), log(2)) # Within a factor 2 using the naive approach
})


test_that("Betabinomial (overdispersion) distribution 3 columns",{
  predictions <- matrix(c(0.999, 0.0005, 0.0005), 1000, 3, byrow = T)
  scores <- matrix(c(0,1,1), 1000, 3, byrow = T)
  observations <- matrix(c(1,0,0), 1000, 3, byrow = T)
  observations[1:10,3] <- 1; observations[1:10, 1] <- 0
  overdispersion <- 0.3
  
  # Naive sampling
  res <- ncdDetectOverdispersion(predictions, scores, overdispersion, N = 100, method = "naive", observations = observations)
  bbinom_p_value <- VGAM::dbetabinom(10, 1000, 0.001, overdispersion^2*0.001/(1-0.001), log = T)
  
  expect_lt(abs(log(res$p_value)-bbinom_p_value), log(2)) # Within a factor 2 using the naive approach
  
  # Numeric Integration
  res <- ncdDetectOverdispersion(predictions, scores, overdispersion, N = 10, method = "numeric", observations = observations)
  bbinom_p_value <- VGAM::dbetabinom(10, 1000, 0.001, overdispersion^2*0.001/(1-0.001), log = T)
  
  expect_lt(abs(log(res$p_value)-bbinom_p_value), log(2)) # Within a factor 2 using the naive approach
})
