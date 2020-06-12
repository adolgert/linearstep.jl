
#' Takes hazard rate per day and time interval and returns probability.
#'
#' @param h hazard per day
#' @param tau time step duration
#' @return probability of event
#' @export
foi2ar <- function(h, tau = 10) {
  1 - exp(-h * tau)
}


#' Takes in probability and time interval and returns rate per day
#'
#' @param alpha a probability (likely the attack rate)
#' @param tau time step duration
#' @return hazard rate per day
#' @export
ar2foi <- function(alpha, tau = 10) {
  -log(1 - alpha) / tau
}


eir2foi <- function(eir, k = 3.8, b = 0.55, tau = 10) {
  log(1 + b * eir * k * tau) / k / tau
}

foi2eir <- function(h, k = 3.8, b = 0.55, tau = 10) {
  (exp(k * h * tau) - 1) / b / k / tau
}

ar2eir <- function(alpha, b = 0.55, k = 3.8, tau = 10) {
  ((1 - alpha)^(-k) - 1) / b / k / tau
}

eir2ar <- function(eir, b = 0.55, k = 3.8, tau = 10) {
  1 - (1 + b * k * eir * tau)^(-1 / k)
}

ek2vc_vec <- function(eir, kappa, pTau) {
  L <- length(eir)
  Vc <- c()
  for (i in 2:L) {
    vc <- (eir[i] * (1 - pTau) - pTau * (1 - pTau) * eir[i - 1]) / kappa[i - 1]
    Vc <- c(Vc, vc)
  }
}

ek2vc <- function(eir, kappa, pTau) {
  diff(tail(eir, 2) * c(1 - pTau, pTau * (1 - pTau))) / tail(kappa, 2)[1]
}

Vc2Vo <- function(Vc, phi) {
}
