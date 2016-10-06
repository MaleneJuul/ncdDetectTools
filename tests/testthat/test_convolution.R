library(ncdDetectTools)
context("convolution")

test_that("Binomial distribution",{
  cv <- convolution(dat = data.table(x = rep(1:100, each = 2), 
                                     y = rep(c(1,0),100), 
                                     probability = rep(c(0.4,0.6), 100)))
  cv_50 <- log(sum(cv[52:101,]$probability))
  cv_60 <- log(sum(cv[62:101,]$probability))
  pb_50 <- pbinom(q = 50, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  pb_60 <- pbinom(q = 60, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  
  expect_lt(abs(cv_50-pb_50), log(1.05) ) # At most 5 % difference in probability space
  expect_lt(abs(cv_60-pb_60), log(1.05) )
})

test_that("Unsorted x",{
  dat <- data.table(x = rep(1:100, each = 2), 
                    y = rep(c(1,0),100), 
                    probability = rep(c(0.4,0.6), 100))
  dat <- dat[sample(200),]
  cv <- convolution(dat = dat)
  cv_50 <- log(sum(cv[52:101,]$probability))
  cv_60 <- log(sum(cv[62:101,]$probability))
  pb_50 <- pbinom(q = 50, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  pb_60 <- pbinom(q = 60, size = 100, prob = 0.4, lower.tail = F, log.p = T)
  
  expect_lt(abs(cv_50-pb_50), log(1.05) ) # At most 5 % difference in probability space
  expect_lt(abs(cv_60-pb_60), log(1.05) )
})