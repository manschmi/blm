context("Testing plots")

seed <- as.integer(1000 * rnorm(1))
test_that(paste("Plots execute silently but diagnostic fails for lm class input", seed), {
  w0 <- 0.3 ; w1 <- 1.1 ; b <- 1.3
  x <- rnorm(50)
  y <- rnorm(50, w1 * x + w0, 1/b)
  mod <- blm(y~x, beta=b, data=data.frame(x=x, y=y))
  
  expect_silent(plot(mod))
  
  expect_silent(diagnostic_plots(mod))
  
  expect_error(diagnostic_plots(lm(y~x)))
  
})

