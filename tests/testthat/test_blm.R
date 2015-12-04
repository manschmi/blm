context("Testing constructor, update and predict")

seed <- as.integer(1000 * rnorm(1))
test_that(paste("Evaluates construction of simple blm with intercept using", seed), {
  
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w1 * x + w0, 1/b)
  
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  expect_is(mod, 'blm')
  
  weights <- coef(mod)
  weights_covar <- covar.blm(mod)
  
  expect_equal(colnames(weights_covar), rownames(weights_covar))
  expect_equal(colnames(weights_covar), names(weights))
  
  x_new <- rnorm(50)
  y_new <- rnorm(50, w1 * x_new + w0, 1/b)
  
  mod_new <- blm(y~x, prior=mod, beta=b, 
                 data=data.frame(x=x_new, y=y_new))
  
  mod_new2 <- update.blm(mod, prior=mod, beta=b, 
                     data=data.frame(x=x_new, y=y_new))
  
  expect_equal(coef(mod_new), coef(mod_new2))
  expect_equal(mod_new$frame, mod_new2$frame)
  
  mod_new3 <- update.blm(mod, formula=y~x+0, prior=mod, beta=b, 
                         data=data.frame(x=x_new, y=y_new))
  
  mod_new4 <- blm(y~x+0, prior=mod, beta=b, 
                 data=data.frame(x=x_new, y=y_new))
  
})

