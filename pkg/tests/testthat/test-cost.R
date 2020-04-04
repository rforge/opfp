library(testthat)
context("cost")

x <- c(-1, 1)
test_that("penalty above 2 gets 1 segment", {
  penalty <- 3
  penalty <- 2.001
  fit <- fpop::Fpop(x, penalty)
  cost.vec <- c(-1, 0)+penalty
  expect_equal(length(fit$t.est), 1)
  expect_equal(fit$cost, cost.vec)
})
test_that("penalty below 2 gets 2 segments", {
  penalty <- 1.999
  fit <- fpop::Fpop(x, penalty)
  cost.vec <- c(-1, -2)+penalty*(1:2)
  expect_equal(length(fit$t.est), 2)
  expect_equal(fit$cost, cost.vec)
})
