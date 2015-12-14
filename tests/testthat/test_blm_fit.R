context("Testing blm fit correctness")

seed <- as.integer(1000 * rnorm(1))
test_that(paste("Deviance of a blm fit can max be equal to lm fit using seed", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w1 * x + w0, 1/b)
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  
  expect_gte(deviance(mod), deviance(lm(y~x, data=data.frame(x=x, y=y))))
  
})


seed <- as.integer(1000 * rnorm(1))
test_that(paste("1/sd of residuals close to precision of data using seed", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 0.1
  x <- rnorm(1000)
  y <- rnorm(1000, w1 * x + w0, 1/b)
  mod1000 <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  
  expect_equal(1/sd(resid(mod1000)), b, tolerance=.1)
  
})


seed <- as.integer(1000 * rnorm(1))
test_that(paste("Distribution of a blm fit reflects the precision of the model using seed", seed), {
  
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(100)
  y <- rnorm(100, w1 * x + w0, 1/b)
  mod13 <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  covar.blm(mod13)
  
  b33 <- 3.3
  y33 <- rnorm(100, w1 * x + w0, 1/b33)
  mod33 <- blm(y33~x, beta=b33, data=data.frame(x=x, y=y33))
  covar.blm(mod33)
  
  #model with more precision should have smaller variance on the coefficients
  expect_equal(all(as.vector(covar.blm(mod13)) > as.vector(covar.blm(mod33))), 
               TRUE)
  
  #model provided with higher precision should have less variance on the 
  #coefficients than equal model provided with higher precision
  mod13b <- blm(y~x, beta=b33, data=data.frame(x=x, y=y))

  expect_equal(all(as.vector(covar.blm(mod13)) > as.vector(covar.blm(mod13b))), 
               TRUE)
  
})



seed <- as.integer(1000 * rnorm(1))
test_that(paste("Precision of predicted values drops with distance to data points using seed", seed), {
  
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(100)
  y <- rnorm(100, w1 * x + w0, 1/b)
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  
  x2 <- 0:100
  y2 <- predict(mod, data.frame(x=x2), report.var=T)
  expect_equal(cor(x2, y2$var, method='spearman'),1)

})
