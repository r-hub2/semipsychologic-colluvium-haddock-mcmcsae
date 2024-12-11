
context("Kronecker product templates")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("kronecker product closure template works", {
  M1 <- Cdiag(c(1, 1, 5))
  M2 <- c(2, 3)
  kron <- build_kron(M1, M2, q2=2)
  M1b <- Cdiag(c(1.5, 1, 6))
  M2b <- c(3, 2.1)
  expect_equal(kron(M1b, M2b), cross(Cdiag(M2b), M1b))
  expect_equal(kron(M1b, M2b)@x, as.numeric(outer(M2b, M1b@x)))
  M1 <- Q_RW1(3L)
  M2 <- 12.1
  kron <- build_kron(M1, M2, q2=1)
  M1b <- M1 + CdiagU(3)
  M2b <- 5.8
  expect_equal(kron(M1b, M2b), M2b * M1b)
  M2 <- c(1.1, 2.1)
  kron <- build_kron(M1, M2, q2=2)
  M2b <- Cdiag(c(3.1, 3.9))
  expect_equal(as(as(kronecker(M1b, M2b), "CsparseMatrix"), "symmetricMatrix"), cross(M2b, M1b))
  expect_equal(kron(M1b, M2b@x), cross(M2b, M1b))
  M2 <- rWishart(1L, 4, diag(4))[,,1L]
  kron <- build_kron(M1, M2, q2=4)
  M2b <- rWishart(1L, 4, diag(4))[,,1L]
  expect_equal(kron(M1b, M2b), as(as(kronecker(M1b, M2b), "CsparseMatrix"), "symmetricMatrix"))
  M1 <- Diagonal(x=runif(4L))
  M2 <- Q_RW2(7L)
  kron <- build_kron(M1, M2, q2=7L)
  M1b <- Diagonal(x=runif(4L))
  M2b <- 2*M2
  expect_equal(kron(M1b, M2b@x), as(as(kronecker(M1b, M2b), "CsparseMatrix"), "symmetricMatrix"))
  M1 <- rWishart(1L, 4, diag(4))[,,1L]
  kron <- build_kron(M1, M2, q2=7L)
  M1b <- rWishart(1L, 4, diag(4))[,,1L]
  expect_equal(kron(M1b, M2b@x), as(as(kronecker(M1b, M2b), "CsparseMatrix"), "symmetricMatrix"))
  M1 <- as.matrix(Q_RW1(5))
  M2 <- 1 + runif(1)
  kron <- build_kron(M1, M2, q2=7)
  M1b <- M1 + CdiagU(5)
  M2b <- M2 + runif(1)
  expect_equal(kron(M1b, M2b), as(as(kronecker(M1b, Cdiag(rep.int(M2b, 7))), "CsparseMatrix"), "symmetricMatrix"))
})

test_that("kronecker product closure template works, fixed M1", {
  M1 <- CdiagU(3)
  M2 <- c(2, 3)
  kron <- build_kron(M1, M2, q2=2, M1.fixed=TRUE)
  M2b <- c(3, 2.1)
  expect_equal(kron(M1, M2b)@x, as.numeric(outer(M2b, rep.int(1, 3))))
  M1 <- as(as(kronecker(Q_RW1(4), Q_iid(3)), "CsparseMatrix"), "symmetricMatrix")
  M2 <- 12.1
  kron <- build_kron(M1, M2, q2=1, M1.fixed=TRUE)
  M2b <- 5.8
  expect_equal(kron(M1, M2b), M2b * M1)
  M2 <- c(0.8, 0.7)
  kron <- build_kron(M1, M2, q2=2, M1.fixed=TRUE)
  M2b <- c(0.9, 5.7)
  expect_equal(kron(M1, M2b), cross(Cdiag(M2b), M1))
  M2 <- rWishart(1L, 4, diag(4))[,,1L]
  kron <- build_kron(M1, M2, q2=4, M1.fixed=TRUE)
  M2b <- rWishart(1L, 4, diag(4))[,,1L]
  expect_equal(kron(M1, M2b), as(as(kronecker(M1, M2b), "CsparseMatrix"), "symmetricMatrix"))
  M1 <- Diagonal(10)
  kron <- build_kron(M1, M2, q2=4L, M1.fixed=TRUE)
  expect_equal(kron(M1, M2b), as(as(kronecker(M1, M2b), "CsparseMatrix"), "symmetricMatrix"))
  M1 <- as.matrix(Q_RW2(7))
  M2 <- 1 + runif(1)
  kron <- build_kron(M1, M2, q2=12, M1.fixed=TRUE)
  M2b <- M2 + runif(1)
  expect_equal(kron(M1, M2b), as(as(kronecker(M1, Cdiag(rep.int(M2b, 12))), "CsparseMatrix"), "symmetricMatrix"))
})

test_that("kronecker product template does not fail with explicit zeros", {
  M1 <- Diagonal(x=c(1,0,3))
  M2 <- matrix(c(2,1,1,2), 2, 2)
  kron <- build_kron(M1, M2, q2=2)
  M2b <- rWishart(1L, 2, diag(2))[,,1L]
  expect_equal(sum(abs(kron(M1, M2b) - kronecker(M1, M2b))), 0)
})
