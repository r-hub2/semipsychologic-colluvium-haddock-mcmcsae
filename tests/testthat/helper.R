
expect_between <- function(object, lower = -Inf, upper = Inf) {
  expect_true(all(object >= lower & object <= upper))
}
