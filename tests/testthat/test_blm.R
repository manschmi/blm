context("Testing constructing blm, update etc on a range of formulas")

#checks expect_equal for coefficients and covariance of the coefficients of 2 blm models
compare_blm <- function(model1, model2) {
  expect_equal(coef(model1), coef(model2))
  expect_equal(covar.blm(model1), covar.blm(model2))
}


#checks return values (mostly names and dimensions) of generic 'getter' functions running on blm models: coef, covar.blm, resid, fitted, confint
check_blm_getter_funs <- function(model) {
  responseless_formula <- delete.response(terms(model$formula))
  var_names <- attr(responseless_formula,"term.labels")
  if ( attr(responseless_formula,"intercept") ) {
    var_names <- c( '(Intercept)', var_names)
  }
  
  cf <- coef(model)
  cov <- covar.blm(model)
  cfi <- confint(model)
  res <- resid(model)
  dev <- deviance(model)
  pred <- predict(model)
  fit <- fitted(model)
  
  expect_equal(names(cf), var_names)
  expect_equal(colnames(cov), var_names)
  expect_equal(rownames(cov), var_names)
  expect_equal(rownames(cfi), var_names)
  
  expect_equal(length(pred), nrow(model$frame))
  expect_equal(length(fit), nrow(model$frame))
  expect_equal(length(res), nrow(model$frame))
  
  expect_equal(fit, pred)
  
  expect_equal( (res + fit), model.response(model$frame) )
  
  expect_gte(dev, 0)
  expect_equal(dev, sum(res^2))
  
}


seed <- as.integer(1000 * rnorm(1))
test_that(paste("Evaluates construction of simple blm with intercept using seed", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w1 * x + w0, 1/b)
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  expect_is(mod, 'blm')
  
  weights <- coef(mod)
  weights_covar <- covar.blm(mod)
  
  expect_equal(colnames(weights_covar), rownames(weights_covar))
  expect_equal(colnames(weights_covar), names(weights))
  
  check_blm_getter_funs(mod)
})


seed <- as.integer(1000 * rnorm(1))
test_that(paste("Evaluates construction and updating of simple blm with intercept using seed", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w1 * x + w0, 1/b)
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  
  ##updating vs constructing a new model
  x_new <- rnorm(50)
  y_new <- rnorm(50, w1 * x_new + w0, 1/b)
  
  mod_new <- blm(y~x, prior=mod, beta=b, 
                 data=data.frame(x=x_new, y=y_new))
  
  mod_new2 <- update(mod, data=data.frame(x=x_new, y=y_new))
  
  compare_blm(mod_new, mod_new2)
  
  check_blm_getter_funs(mod_new)
  check_blm_getter_funs(mod_new2)
  
})

seed <- as.integer(1000 * rnorm(1))
test_that(paste("Evaluates construction and updating of simple blm with intercept to model without intercept using seed", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w1 * x + w0, 1/b)
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  
  ##updating vs constructing to a new model wo intercept
  mod2 <- update(mod, formula=y~x+0)
  
  mod3 <- blm(y~x+0, prior=mod, beta=b, 
                 data=data.frame(x=x, y=y))
  
  expect_equal(length(coef(mod2)), 1)
  expect_equal(length(coef(mod3)), 1)
  expect_equal(dim(covar.blm(mod2)), NULL)
  expect_equal(dim(covar.blm(mod3)), NULL)
  
  compare_blm(mod2, mod3)
  
  #check_blm_getter_funs(mod2) ##unnamed covar matrix for single weight models
  #check_blm_getter_funs(mod3) ##unnamed covar matrix for single weight models
  
})


seed <- as.integer(1000 * rnorm(1))
test_that(paste("Evaluates construction of polynomial models using seed", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; w2 <- 3.3 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w2 * x^2 + w1 * x + w0, 1/b)
  
  
  mod1 <- blm(y~poly(x,2,raw=T), beta=b, data=data.frame(x=x, y=y))
  mod2 <- blm(y~x+I(x^2), beta=b, data=data.frame(x=x, y=y))
  
  expect_equal(as.numeric(coef(mod1)), as.numeric(coef(mod2)))
  expect_equal(as.numeric(covar.blm(mod1)), as.numeric(covar.blm(mod2)))
  
  #check_blm_getter_funs(mod2) -> does not work for single weight models 
})


seed <- as.integer(1000 * rnorm(1))
test_that(paste("Evaluates models with 2 explanatory variables using seed", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; w2 <- 3.3 ; b <- 1.3
  x <- rnorm(50)
  z <- rnorm(50)
  y <- rnorm(50, w2 * z + w1 * x + w0, 1/b)
  
  
  mod <- blm(y~x+z, beta=b, data=data.frame(x=x, y=y, z=z))
  
  check_blm_getter_funs(mod)

  #simplify the model
  mod2 <- blm(y~x, prior=mod, beta=b, data=data.frame(x=x, y=y))
  
  #we actually do not need to specify 'data', it is fetched from environment, but this is probably a bad idea
  mod3 <- blm(y~x, prior=mod, beta=b) 
  
  compare_blm(mod2,mod3)
  
  #we can also use update
  mod4 <- update(mod, formula=y~x)
  
  compare_blm(mod2,mod4)
})


seed <- as.integer(1000 * rnorm(1))
test_that(paste("Evaluates updating when removing intercept from formula using", seed), {

  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w1 * x + w0, 1/b)
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  
  expect_message(update(mod, y~x+0, prior=mod$prior), "prior contains more variables than the model")
  
  new_mod <- blm(y~x+0, prior=mod$prior, beta=b)
  new_mod2 <- update(mod, y~x+0, prior=mod$prior) 
  new_mod3 <- update(mod, ~.+0, prior=mod$prior) 
  
  compare_blm(new_mod, new_mod2)
  compare_blm(new_mod, new_mod3)
})