library(ncdDetectTools)
context("ncdDetectOverdispersion")

test_that("Betabinomial (overdispersion) distribution 2 columns",{
  predictions <- matrix(c(0.999, 0.001), 1000, 2, byrow = T)
  scores <- matrix(c(0,1), 1000, 2, byrow = T)
  observed_score <- 10
  overdispersion <- 0.3
  
  res <- ncdDetectOverdispersion(predictions, scores, overdispersion, N = 100, method = "naive", observed_score = observed_score)[type == "naive", p_value]
  bbinom_p_value <- VGAM::dbetabinom(10, 1000, 0.001, overdispersion^2*0.001/(1-0.001), log = T)

  expect_lt(abs(log(res)-bbinom_p_value), log(3)) # Within a factor 2 using the naive approach
})


test_that("Betabinomial (overdispersion) distribution 3 columns",{
  predictions <- matrix(c(0.999, 0.0005, 0.0005), 1000, 3, byrow = T)
  scores <- matrix(c(0,1,1), 1000, 3, byrow = T)
  observed_score <- c(5, 10)
  overdispersion <- 0.3
  
  # Naive method with N = 10,000
  # score      p_value              type
  # 1:     5 7.346926e-03             naive
  # 2:    10 2.702902e-06             naive
  # 3:     5 3.636878e-03 no_overdispersion
  # 4:    10 1.074283e-07 no_overdispersion
  
  # Numeric Integration
  res_numer <- ncdDetectOverdispersion(predictions, scores, overdispersion, N = 100, method = "numeric", observed_score = observed_score)

  # Expect overdispersion correction gives less significant p-values
  res_numer_5  <- res_numer[score == 5 & type == "numeric", p_value]
  res_numer_10 <- res_numer[score == 10 & type == "numeric", p_value]
  expect_gt( res_numer_5, res_numer[score == 5 & type == "no_overdispersion", p_value] )
  expect_gt( res_numer_10, res_numer[score == 10 & type == "no_overdispersion", p_value] )
  
  # Expect reasonably close to value obtained by naive sampling
  expect_lt( abs(log(res_numer_5)  - log(7.346926e-03)), log(1.50)) # With-in 50% of the correct value
  expect_lt( abs(log(res_numer_10) - log(2.702902e-06)), log(1.50))
})
