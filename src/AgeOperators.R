
require(Matrix)

shiftO <- function(N) {
  bandSparse(N, k = 1) * 1
}

stageO <- function(N, p) {
  bandSparse(N, k = 0) * (1 - p) + bandSparse(N, k = 1) * p
}

twoBlockO <- function(A, B, p) {
  nA <- dim(A)[1]
  calO <- bdiag(A, B)
  calO[nA, nA + 1] <- p
  calO
}

ages_1066 <- function(unit = 365) {
  x1 <- c(1:18) * 10
  x2 <- max(x1) + c(1:18) * (30 + 1 / 1.8)
  x3 <- max(x2) + c(1:18) * 365
  x4 <- max(x3) + c(1:12) * 365 * 5
  c(x1, x2, x3, x4) / unit
}

calO_1066 <- function() {
  A <- shiftO(18)
  B <- stageO(18, 1 / 3)
  C <- stageO(18, 1 / 36.5)
  D <- stageO(12, 1 / 36.5 / 5)
  calO <- twoBlockO(A, B, 1 / 3)
  calO <- twoBlockO(calO, C, 1 / 36.5)
  calO <- twoBlockO(calO, D, 1 / 36.5 / 5)
  return(calO)
}

ages_u10 <- function(unit = 365) {
  x1 <- c(1:3) * 10
  x2 <- max(x1) + c(1:11) * (30 + 5 / 11)
  x3 <- max(x2) + c(1:9) * 365
  x4 <- max(x3) + c(1:2) * 365 * 5
  c(x1, x2, x3, x4) / unit
}

calO_u10 <- function() {
  A <- shiftO(3)
  B <- stageO(11, 1 / 3)
  C <- stageO(9, 1 / 36.5)
  D <- stageO(2, 1 / 36.5 / 5)
  calO <- twoBlockO(A, B, 1 / 3)
  calO <- twoBlockO(calO, C, 1 / 36.5)
  calO <- twoBlockO(calO, D, 1 / 36.5 / 5)
  return(calO)
}

constructBirthO <- function(N) {
  nn <- matrix(0, N, N)
  nn[1, 1] <- 1
  return(nn)
}
