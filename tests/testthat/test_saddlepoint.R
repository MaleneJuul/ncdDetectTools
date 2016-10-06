library(ncdDetectTools)
context("saddlepoint")

test_that("Binomial distribution",{
  sp <- saddlepoint(51, dat = data.table(x = rep(1:100, each = 2), 
                                         y = rep(c(1,0),100), 
                                         probability = rep(c(0.4,0.6), 100)))
  pb <- pbinom(q = 50, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  
  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space
  
  sp <- saddlepoint(61, dat = data.table(x = rep(1:100, each = 2), 
                                         y = rep(c(1,0),100), 
                                         probability = rep(c(0.4,0.6), 100)))
  pb <- pbinom(q = 60, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  
  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space
})

test_that("Unsorted x",{
  set.seed(1)
  
  dat <- data.table(x = rep(1:100, each = 2), 
                    y = rep(c(1,0),100), 
                    probability = rep(c(0.4,0.6), 100))
  dat <- dat[sample(200),]
  
  sp <- saddlepoint(61, dat = dat)
  pb <- pbinom(q = 60, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  
  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space
})