library(ncdDetectTools)
context("saddlepoint")

test_that("Binomial distribution",{
  sp <- saddlepoint(51, dat = data.table(y = rep(1,100), probability = 0.4))
  pb <- pbinom(q = 50, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  
  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space
  
  sp <- saddlepoint(61, dat = data.table(y = rep(1,100), probability = 0.4))
  pb <- pbinom(q = 60, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  
  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space
})