library(testthat)
context("cost")

x <- c(-1, 1)
test_that("penalty above 2 gets 1 segment", {
  fit <- fpop::Fpop(x, 2.0001)
  expect_equal(length(fit$t.est), 1)
})
test_that("penalty below 2 gets 2 segments", {
  fit <- fpop::Fpop(x, 1.9999)
  expect_equal(length(fit$t.est), 2)
})
