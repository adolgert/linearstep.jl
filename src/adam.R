#' The Adam model is the least complicated model that can support operational needs.
#'
#' For documentation on how or why this model is built this way, look at the
#' RAMP Model Library. This description is a guide to names and functions in the code.
#' This is a discrete-time model which defaults to a ten-day time step.
#' We build up probability matrices and apply them to the state space.
#'
#' The \emph{compartmental state space} for a cohort has 22 states
#' \enumerate{
#'   \item U - uninfected
#'   \item X1-X9 - infected
#'   \item E0 - Exposed, who were previously uninfected.
#'   \item E1-E9 - exposed, meaning active infection plus superinfection.
#'   \item P1, P2 - people treated one or two steps ago.
#' }
#' The system state is larger. The function \code{\link{emptyX_Adam}} creates
#' a system state for a cohort, while \code{\link{emptyXX_Adam}} creates system
#' state for multiple age groups at once. That's \code{X} for cohort versus \code{XX}
#' for age-segmented. System state includes
#' \itemize{
#'   \item \code{X} or \code{XX}, compartmental state, above
#'   \item \code{pit} which is a vector to track immune states.
#'   \item \code{ageInDays} Well, the age in days.
#'   \item \code{V}
#'   \item \code{alpha2} Attack rate from two steps ago.
#'   \item \code{alpha1} Attack rate from one step ago.
#' }
#' Compartmental state is updated with a matrix multiplication and all others are
#' updated separately.
#'
#' When we apply the Adam model forward in time, it takes an attack rate, \code{ar},
#' and a system state, \code{XX}, and produces a new system state. This forward
#' function is \code{\link{PtimesX_Adam}}. It applies a forward-time Markov matrix
#' that comes from \code{\link{calP_Adam}}. That matrix has a block form constructed
#' of two matrices, \eqn{B_1} which derives from \code{\link{SuSblock_Adam}} and
#' \eqn{B_2} which derives from \code{\link{calP_SoI}}. Here, \eqn{\alpha} is attack
#' rate and \eqn{\chi} is the cumulant form \eqn{\chi=P[X_i|U]}. The vectors
#' \code{(T_1, T_2)} relate to treatment rates and bring people into the treated class.
#' \deqn{
#' \mathcal{P} = \left[\begin{array}{ccc}
#' (1-\alpha\chi)B_1 & (1-\alpha)B_2 & 0 & 1 \\
#' \alpha\chi B_1 & \alpha B_2 & 0 & 0 \\
#' T_1 & T_2 & 0 & 0 \\
#' 0 & 0 & 1 & 0
#' \end{array}\right]
#' }
#' The net result is a length 22 square matrix.
#'
#' If you don't know the previous state of the model, then there is a function
#' to find steady-state conditions that would have resulted in the given
#' attack rate. The function \code{\link{ar2stableWave_Adam}} iterates until
#' it converges on a system state.
#'
#' The forward step allows these transitions:
#' \itemize{
#'   \item U -> E0
#'   \item E0 -> X1
#'   \item Xi -> Xi+1
#'   \item Xi -> U
#'   \item Xi -> Ei
#'   \item Ei -> P1
#'   \item Xi -> P1
#'   \item P1 -> P2
#'   \item P2 -> U
#' }
#' Most of the work of constructing the model is controlling how the rates of these
#' transitions change with age and intervention.
#'
#' We usually report PfPR, not system state. This is denoted as \code{pr},
#' so the function \code{\link{XX2pr29}} turns system state into PfPR for
#' 2 to 9 year-olds.
#'
#' We also apply the Adam model in reverse, using it to infer AR from PfPR.
#' The function \code{\link{pr2ar_Adam}} does the reverse step.
#'
#' This code implements detection with projection operators. These
#' are vectors, where you take the dot product of the vector and the
#' compartmental state in order to find a scalar value. These are
#' denoted as \code{calD...} or \code{calT...} for detection or treatment.
#'
NULL


require(Matrix)
source("MetricConversions.R")
source("AgeWeights.R")
source("AgeOperators.R")


#' Stage-up function for superinfection.
#' The fraction of new infections entering each stage. So it's
#' \eqn{P[X_i|U]}.
#' @param x A single immunity level, where 0 is no immunity and 1 is high immunity.
#' @param k1 Shape parameter for distribution.
#' @param k2 Shape parameter for distribution.
#' @param N The number of stages of infection that could be entered.
#' @return A conditional probability to enter a stage given infection for N - 1 stages.
#' @export
xiF.x <- function(x, k1 = 3, k2 = 0.4, N = 9) {
  diff(pbeta(c(0:(N - 1)) / (N - 1), 1 - exp(-k1 * x) + k2 * x * (1 - x), exp(-k1 * x) + k2 * x^2))
}


#' Stage-up function for superinfection.
#' @param x a Single immunity level
#' @param k1 a shape parameter
#' @return Conditional probability for N - 1 infection stages.
#' @seealso \code{\link{xiF.x}}
xiF_0 <- function(x, k1 = .008, N = 9) {
  diff(pbeta(c(0:(N - 1)) / (N - 1), 1 - exp(-k1 * x), exp(-k1 * x)))
}


#' Stage-up function for superinfection.
#' @param x a single immunity level.
#' @param p Unused
#' @param N The number of infection stages.
#' @return Probability to enter a stage off infection for N stages.
#' @seealso \code{\link{xiF.x}}
xiF_1 <- function(x, p = 1, N = 9) {
  diff(pbeta(c(0:N) / N, x, 1 - x))
}


#' Stage-up function for superinfection.
#' @seealso \code{\link{xiF.x}}
xiF.2 <- function(V, k1 = .2, k2 = .2) {
  xi1 <- exp(-k1 * V[1])
  xi2 <- xiF_0(k1 * V[2], k2, 8)
  xi <- c(xi1, (1 - xi1) * xi2)
  return(xi)
}


#' Stage-up function for superinfection.
#' @seealso \code{\link{xiF.x}}
xiF.3 <- function(V, k1 = 5, p = 3, N = 9) {
  xi1 <- exp(-k1 * V)
  xi2 <- xiF_1(V, p, N - 1)
  xi <- c(xi1, (1 - xi1) * xi2)
  return(xi)
}


#' Stage-up function for superinfection.
#' smitdave, do you want this to be k1=.3, k2=.2? or .1 and .1?
#'
#' @param pit Two varaibles, of which the first is used.
#' @param k1 A shape parameter for the distribution.
#' @param k2 A shape parameter for the distribution.
#' @param p Unused.
#' @param N The number of stages of infection.
#' @return
#' @seealso \code{\link{xiF.x}}
xiF_h <- function(pit, k1 = .3, k2 = .2, p = 2, N = 9) {
  # smitdave - Should this be pit[1] or pit[2]? It's different in the two files.
  xi1 <- exp(-k1 * pit[1])
  xi2 <- xiF_0(pit[1], k2, N - 1)
  xi <- c(xi1, (1 - xi1) * xi2)
  return(xi)
}




#' SuS is stage-up superinfection.
#'
#' @parameter par Parameters including N.
#' @return An (N+1) square matrix.
#'
#' This matrix gives the shape of an operator applied
#' to the first set of ten states.
#' The returned matrix is 0 or 1 values in this shape.
#' \deqn{
#' \left[\begin{array}{cccccccccc}
#' . & . & . & . & . & . & . & . & . & . \\
#' 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
#' 1 & . & 1 & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
#' 1 & . & . & 1 & 1 & 1 & 1 & 1 & 1 & 1 \\
#' 1 & . & . & . & 1 & 1 & 1 & 1 & 1 & 1 \\
#' 1 & . & . & . & . & 1 & 1 & 1 & 1 & 1 \\
#' 1 & . & . & . & . & . & 1 & 1 & 1 & 1 \\
#' 1 & . & . & . & . & . & . & 1 & 1 & 1 \\
#' 1 & . & . & . & . & . & . & . & 1 & 1 \\
#' 1 & . & . & . & . & . & . & . & . & 1
#' \end{array}\right]
#' }
#' This matrix encodes which states can transition to
#' which other states, but not their rates.
#'
#' @export
SuSblock_Adam <- function(par) {
  with(par, {
    SuS <- bandSparse(N, k = ((1:N) - 1))
    SuS <- rbind(0, cbind(1, SuS))
    return(SuS)
  })
}


#' Creates a chemoprotection matrix that shows stages of infection.
#'
#' @param params Parameters with N, r, nu.
#' @return A square matrix of size (N+1), where N=9.
#'
#' This matrix is the shape of an operator to apply to
#' the second set of ten states.
#' The transitions here are in terns of \eqn{s=r(1-\nu)}
#' and \eqn{v=r\nu}.
#' \deqn{
#' \left[\begin{array}{cccccccccc}
#' 1 & 1-s-v & 1-s-v & 1-s-v & 1-s-v & 1-s-v & 1-s-v & 1-s-v & 1-s-v & 1-s-v \\
#' . & s & . & . & . & . & . & . & . & . \\
#' . & v & s & . & . & . & . & . & . & . \\
#' . & . & v & s & . & . & . & . & . & . \\
#' . & . & . & v & s & . & . & . & . & . \\
#' . & . & . & . & v & s & . & . & . & . \\
#' . & . & . & . & . & v & s & . & . & . \\
#' . & . & . & . & . & . & v & s & . & . \\
#' . & . & . & . & . & . & . & v & s & . \\
#' . & . & . & . & . & . & . & . & v & s
#' \end{array}\right]
#' }
#' This isn't, by itself, a transition matrix because the last
#' column doesn't sum to 1.
#'
#' @export
calP_SoI <- function(par) {
  with(par, {
    s <- r * (1 - nu)
    v <- r * nu
    # Remain in state with probability s, go to next state probability v
    PSoI <- bandSparse(N, k = 0) * s + t(bandSparse(N, k = 1) * v)
    PSoI <- rbind(1 - r, PSoI)
    PSoI <- cbind(0, PSoI)
    # Uninfected remain uninfected.
    PSoI[1, 1] <- 1
    return(PSoI)
  })
}


#' Chemoprotective transition matrix for the compartmental state.
#'
#' @param alpha Attack rate
#' @param pit Tracking variables
#' @param V I don't know.
#' @param par Parameters for the Adam model.
#' @return a square matrix of length 22.
calP_Adam <- function(alpha, pit, V, par) {
  with(par, {
    if (!exists("SuS")) {
      SuS <<- SuSblock_Adam(par)
      warning("Assigning to global SuS")
    }
    if (!exists("PSoI")) {
      PSoI <<- calP_SoI(par)
      warning("Assigning to global PSoI")
    }
    alpha <- min(alpha, 1)

    # xi is fraction of new infections entering ith state where i is 1:9.
    xi <- xiF(pit, V)
    # Multiply rows by the values in xi. Equivalent to diag(c(0, xi)) %*% SuS.
    E1 <- SuS * c(0, xi)
    # The next 3 lines ensure each column total sums to 1. It normalizes the rates
    # given by xi = P[X_i|U], when not every destination is available.
    e1 <- colSums(E1)
    e1[which(e1 == 0)] <- 1
    hatxi <- E1 %*% diag(1 / e1)

    # rho is base fever rate, rhoXstage is fever by stage,
    # and etaXstage is treatment rate. etaXstage and rhoXstage are length(N+1).
    # Deta = rho*etaF(V)*etaXstage*rhoXstage
    Deta <- pmin(rho * etaXstage * rhoXstage, 1)
    # Dd is background rate of taking drugs, as a probability.
    vecD <- Deta + (1 - Deta) * Dd
    # Not treating means you remain in the state. This is (N+1)x(N+1)
    notTreat <- diag(1 - vecD) # 10 x 10
    treat <- c(vecD, vecD) # length 20

    B1 <- PSoI %*% notTreat
    B2 <- hatxi %*% notTreat
    chi <- c(1, cumsum(xi))
    calP <- rbind(
      cbind((1 - alpha * chi) * B1, (1 - alpha) * B2),
      cbind(alpha * chi * B1, alpha * B2),
      treat,
      0 * treat
    )

    calP <- cbind(calP, 0)
    calP <- cbind(calP, 0)
    calP[22, 21] <- 1
    calP[1, 22] <- 1

    return(as.matrix(calP))
  })
}



#' Generate the next state given the current state.
#'
#' @param alpha The attack rate.
#' @param vecX A state vector, as a list, with 22 entries in X,
#'     the PIT variables V1 and V2, and attack rates alpha1 and alpha2
#'     for one and two steps behind the current values.
#' @param par Parameters, not to be confused with the `par()` command.
#' @return A new state vector.
#' @export
PxX_Adam <- function(alpha, vecX, par) {
  with(c(par, vecX), {
    ageInDays <- ageInDays + tau
    w <- wF(ageInDays / 365)
    Xt <- calP_Adam(w * alpha, pit, V, par) %*% X

    # These lines define the PIT variables.
    V <- alpha + V * (1 - mu1)
    pit[1] <- pit[1] * (1 - mu1) + alpha2 # The same as V.
    pit[2] <- pit[2] * (1 - mu2) + (calDimm %*% X)
    pit[3] <- pit[3] * (1 - mu3) + alpha2 * ageXimmune(ageInDays)

    alpha2 <- alpha1
    alpha1 <- w * alpha
    aid <- ageInDays
    return(list(X = Xt, pit = pit, ageInDays = aid, V = V, alpha2 = alpha2, alpha1 = alpha1))
  })
}


#' Creates an empty epidemiological state vector.
#' @return A list containing X with 22 entries, V1, V2, alpha1, alpha2.
#' @export
emptyX_Adam <- function() {
  X <- matrix(0, 22, 1)
  X[1] <- 1
  return(list(X = X, pit = c(0, 0, 0), V = 0, ageInDays = 0, alpha2 = 0, alpha1 = 0))
}


#' emptyXX creates an empty array of epidemiological state vectors.
#' @param age_group_cnt The number of age groups, of various sizes.
#' @return System state as a list.
#' @export
emptyXX_Adam <- function(L = 66) {
  X <- matrix(0, 22, L)
  X[1, ] <- 1
  pit <- matrix(0, 3, L)
  ageInDays <- rep(0, L)
  return(list(X = X, pit = pit, ageInDays = ageInDays, V = 0, alpha2 = 0, alpha1 = 0))
}


#' Simulates a cohort from birth over
#' time steps and returns it as an XX object
#' @param alpha A single constant attack rate.
#' @param params Parameters
#' @time_step_cnt The number of time steps to take.
#' @return A state at the end.
#' @export
cohortXX_Adam <- function(alpha, par, mx = 2920) {
  vecX <- emptyX_Adam()
  XX <- emptyXX_Adam(mx)
  vecX$V <- alpha / par$mu1
  for (i in 1:mx) {
    XX$X[, i] <- vecX$X
    XX$pit[, i] <- vecX$pit
    XX$ageInDays[i] <- vecX$ageInDays
    vecX <- PxX_Adam(alpha, vecX, par)
  }
  XX$alpha2 <- vecX$alpha2
  XX$alpha1 <- vecX$alpha1
  XX$V <- vecX$V
  return(XX)
}

#' PfPR by age and exposure
pfprXage <- function(alpha, rho, par, rTime = 3 * 365) {
  par$rho <- rho
  cXX <- cohortXX_Adam(alpha, par, rTime)
  x <- par$calDlm %*% cXX$X

  plot(cXX$ageInDays / 365, x, type = "l", col = "blue", xlab = "Age (Years)")
  segments(10, 0, 10, 1, col = "red")
  return(cXX)
}


PtimesX_Adam <- function(alpha, XX, par) {
  with(c(par, XX), {
    ageInDays <- (ageInDays + tau) %*% calO
    Xt <- X
    for (a in 1:66) {
      Xt[, a] <- calP_Adam(alpha * w[a], pit[, a], V, par) %*% X[, a]
    }
    Xt <- Xt %*% calO
    Xt[1, 1] <- 1

    V <- alpha + V * (1 - mu1)
    pit[1, ] <- pit[1, ] * (1 - mu1) + alpha
    pit[2, ] <- pit[2, ] * (1 - mu2) + (calDimm %*% X)
    pit[3, ] <- pit[3, ] * (1 - mu3) + alpha2 * aXi

    pit <- pit %*% calO
    ageInDays <- ageInDays %*% calO
    alpha2 <- alpha1
    alpha1 <- alpha

    return(list(X = Xt, pit = pit, ageInDays = ageInDays, V = V, alpha2 = alpha2, alpha1 = alpha1))
  })
}


#' Pr for ages 2 to 9.
XX2pr29_Adam <- function(XX, par) {
  with(par, {
    as.numeric(calDlm %*% XX$X %*% demog29)
  })
}


kappaF_Adam <- function(XX, par) {
  with(par, {
    as.numeric((calK %*% XX$X) %*% (w * demog))
  })
}



# ` Compute PfPR quickly using the cohort utility
fastPfPR_Adam <- function(alpha, par) {
  XX <- cohortXX_Adam(alpha, par, 365)
  x <- par$calDlm %*% XX$X
  mean(x[74:365])
}


#' cohort2ages bins data from a cohort matrix and bins it into ages
cohort2ages_Adam <- function(cXX, par) {
  with(par, {
    aa <- par$tau * (1:dim(cXX$X)[2]) / 365
    XX <- emptyXX_Adam()
    XX$V <- cXX$V
    XX$ageInDays <- rep(0, 66)
    ix <- which(aa < ages[1])
    if (length(ix) == 1) {
      XX$X[, 1] <- cXX$X[, ix]
      XX$pit[, 1] <- cXX$pit[, ix]
      XX$ageInDays[1] <- cXX$ageInDays[ix]
    }
    if (length(ix) > 1) {
      XX$X[, 1] <- rowSums(cXX$X[, ix]) / length(ix)
      XX$pit[, 1] <- rowSums(cXX$pit[, ix]) / length(ix)
      XX$ageInDays[1] <- mean(cXX$ageInDays[ix])
    }
    for (i in 2:length(ages)) {
      ix <- which(aa >= ages[i - 1] & aa < ages[i])
      if (length(ix) == 1) {
        XX$X[, i] <- cXX$X[, ix]
        XX$pit[, i] <- cXX$pit[, ix]
        XX$ageInDays[i] <- cXX$ageInDays[ix]
      }
      if (length(ix) > 1) {
        XX$X[, i] <- rowSums(cXX$X[, ix]) / length(ix)
        XX$pit[, i] <- rowSums(cXX$pit[, ix]) / length(ix)
        XX$ageInDays[i] <- mean(cXX$ageInDays[ix])
      }
    }
    return(XX)
  })
}



ar2stableWave_Adam <- function(alpha, par, tol = 0.01) {
  cXX <- cohortXX_Adam(alpha, par, 2920)
  XX <- cohort2ages_Adam(cXX, par)
  for (i in 1:100) XX <- PtimesX_Adam(alpha, XX, par)
  x <- XX2pr29_Adam(XX, par)
  df <- 1
  while (df > tol) {
    XXl <- XX
    XX <- PtimesX_Adam(alpha, XX, par)
    x <- XX2pr29_Adam(XX, par)
    df <- sum(abs(XXl$X - XX$X))
  }
  XX
}


ar2pr_Adam <- function(alpha, par) {
  with(par, {
    XX <- ar2stableWave_Adam(alpha[1], par)
    XX2pr29_Adam(XX, par)
  })
}



#' Calculate PfPR from attack rate for a sequence of attack rates.
#'
#' @param ar A vector of attack rates
#' @param params An Adam parameters list.
#' @return A vector of pr values for 2-9 year olds.
#' @export
ar2pr_series_Adam <- function(ar, params) {
  if (!is.finite(ar[1])) {
    stop(paste("The first ar is not a finite value.", ar[1], "\n"))
  }
  pr <- ar2stableWave_Adam(ar[1], params)
  youth_initial_pr <- XX2pr29_Adam(pr, params)
  youth_pr <- rep(youth_initial_pr, length(ar))
  for (time_idx in 2:length(ar)) {
    pr <- PtimesX_Adam(ar[time_idx], pr, params)
    youth_pr[time_idx] <- XX2pr29_Adam(pr, params)
  }
  youth_pr
}


ar2stableDiagnostic <- function(alpha, par, XX0 = NULL) {
  if (is.null(XX0)) {
    XX <- emptyXX_Adam()
  } else {
    XX <- XX0
  }
  for (i in 1:365) {
    XX <- PtimesX_Adam(alpha, XX, par)
  }
  x <- XX2pr29_Adam(XX, par)
  list(x = x, XX = XX)
}


calY <- function(V, kappa, YY, par) {
  with(YY, {
    Et <- q * (1 - q) * Z + p * V * kappa
    Zt <- q * Z + (1 - p) * V * kappa
    list(E = T, Z = Zt)
  })
}


zetaF <- function(YY, par) {
  YY$E / par$totalH
}


#' Find the attack rate for a given PfPR, using algebra.
#'
#' @param x2 PfPR in 2-9 year olds
#' @param XX system state for an age-grouped system.
#' @param par Parameters for an Adam model.
#' @return A single value for the attack rate, alpha.
peekAhead_Algebra <- function(x2, XX, par) {
  with(c(XX, par), {
    X0 <- X
    X1 <- X
    dummy <- 0.1
    age_cnt <- 66 # The number of age groups.
    for (a in 1:age_cnt) {
      calP <- calP_Adam(alpha2, c(V1[a], V2[a]), par)
      calP0 <- calP_Adam(0, c(V1[a], V2[a]), par)
      calPalpha <- (calP - calP0) / dummy
      X0[, a] <- calP0 %*% X[, a]
      X1[, a] <- w[a] * calPalpha %*% X[, a]
    }
    V1 <- V1 * (1 - mu1) + calDimm %*% X0
    V2 <- V2 * (1 - mu2) + w * alpha2
    X0 <- X0 %*% calO
    X1 <- X1 %*% calO
    X00 <- X0 * 0
    X10 <- X0 * 0
    V1 <- V1 %*% calO
    V2 <- V2 %*% calO

    for (a in 1:age_cnt) {
      calP0 <- calP_Adam(0, c(V1[a], V2[a]), par)
      X00[, a] <- calP0 %*% X0[, a]
      X10[, a] <- calP0 %*% X1[, a]
    }
    X00 <- X00 %*% calO
    X10 <- X10 %*% calO
    alpha <- (x2 - calDlm %*% X00 %*% demog29) / (((calDlm %*% X10)) %*% demog29)

    return(as.numeric(alpha))
  })
}


#' Find the attack rate that produces the given PR in 2-9 year olds using optimization.
#'
#' @param x2 PfPR in 2-9 year olds.
#' @param XX State of the system, by age. Not the cohort state.
#' @param par Parameters for the Adam model.
#' @return A single value for the minimum attack rate, alpha.
peekAhead_Adam <- function(x2, XX, par) {
  with(c(XX, par), {
    errF <- function(alpha) {
      XX1 <- PtimesX_Adam(alpha, XX, par)
      XX2 <- PtimesX_Adam(alpha, XX1, par)
      (x2 - XX2pr29_Adam(XX2, par))^2
    }
    alpha0 <- 2 * alpha1 - alpha2
    inits <- pmax(0, pmin(c(alpha0 / 3, alpha0 * 3), 1))
    optimize(errF, inits)$min
  })
}



#' Find ar for this PR, assuming Steady State (SS).
pr2arSS_Adam <- function(x, par) {
  errFgood <- function(alpha) {
    (ar2pr_Adam(alpha, par) - x)^2
  }
  errFfast <- function(alpha) {
    (fastPfPR_Adam(alpha, par) - x)^2
  }
  alpha0 <- optimize(errFfast, c(0, 1))$min
  # optimize(errFgood,c(alpha0/2,alpha0*2))$min
  return(alpha0)
}


ar2SW_Adam <- function(alpha, par, mx = 2920) {
  XX <- emptyXX_Adam()
  for (i in 1:mx) XX <- PtimesX_Adam(alpha, XX, par)
  XX
}


#' Given a sequence of PfPR, infer what attack rate would have generated that PfPR.
#'
#' @param xx A vector of PfPR, greater than 3.
#' @param par a list of parameters.
#' @return A list with attack rate and XX.
#' @export
pr2ar_Adam <- function(xx, par, alpha0 = NULL, XX0 = NULL, tol = 0.01) {
  if (is.null(alpha0)) {
    alpha <- pr2arSS_Adam(xx[1], par)
  } else {
    alpha <- alpha0
  }
  if (is.null(XX0)) {
    XX <- ar2stableWave_Adam(alpha, par, tol)
  } else {
    XX <- XX0
  }
  alphaT <- alpha
  for (t in 3:length(xx)) {
    alpha <- peekAhead_Adam(xx[t], XX, par)
    XX <- PtimesX_Adam(alpha, XX, par)
    alphaT <- c(alphaT, alpha)
  }
  return(list(alpha = alphaT, XX = XX))
}


#' Calculate burden from the state of the system.
#'
#' @param XX a system state for the Adam model.
#' @param par parameters for the Adam model.
#' @return A named vector with health states, such as fever or severe malaria.
#' @export
XX2burden_Adam <- function(XX, par) {
  with(par, {
    fever <- calBfever %*% XX$X
    severe <- calBsevere %*% XX$X
    trtdmal <- calBtrtdmal %*% XX$X
    untrmal <- calBuntrmal %*% XX$X
    trtsev <- calBtrtdsev %*% XX$X
    untsev <- calBuntrsev %*% XX$X

    return(c(
      fever = sum(fever * demog),
      U5fever = sum(fever * demog * U5),
      O5fever = sum(fever * demog * O5),
      treatedFever = sum(trtdmal * demog),
      U5trtdFever = sum(trtdmal * demog * U5),
      O5trtdFever = sum(trtdmal * demog * O5),
      untreatedFever = sum(untrmal * demog),
      U5untrdFever = sum(untrmal * demog * U5),
      O5untrdFever = sum(untrmal * demog * O5),
      severeMalaria = sum(severe * demog),
      U5sevMal = sum(severe * demog * U5),
      O5sevMal = sum(severe * demog * O5),
      treatedSevere = sum(trtsev * demog),
      U5trSevMal = sum(trtsev * demog * U5),
      O5trSevMal = sum(trtsev * demog * O5),
      untreatedSevere = sum(untsev * demog),
      U5unSevMal = sum(untsev * demog * U5),
      O5unSevMal = sum(untsev * demog * O5)
    ))
  })
}


#' Calculate burden, given attack rate.
#'
#' @param alpha Attack rate on a ten-day time step. This vector has an entry
#'     for each time step.
#' @param par parameters for the adam model.
#' @param XX0 An initial pfpr. If you don't supply this, the function will estimate
#'     a value from steady state beforehand. That estimate can take a while.
#' @return A data.frame with burden estimates. The same length as input alpha.
#' @export
ar2burden_Adam <- function(alpha, par, XX0 = NULL) {
  with(par, {
    if (is.null(XX0)) {
      XX <- ar2stableWave_Adam(alpha[1], par, .01)
    } else {
      XX <- XX0
    }

    burden <- XX2burden_Adam(XX, par)
    for (t in 2:length(alpha)) {
      XX <- PtimesX_Adam(alpha[t], XX, par)
      burden <- rbind(burden, XX2burden_Adam(XX, par))
    }
    return(data.frame(burden))
  })
}


ar2History <- function(alpha, par, XX = NULL) {
  if (is.null(XX)) {
    XX <- ar2stableWave_Adam(alpha, par)
  }
  eir <- ar2eir(alpha)
  kappa <- kappaF_Adam(XX, par)
  kappaT <- kappa
  for (i in 1:length(alpha)) {
    XX <- PtimesX_Adam(alpha[i], XX, par)
    kappa <- kappaF_Adam(XX, par)
    kappaT <- c(kappaT, kappa)
  }
  V <- eir / kappa
  history <- list(x = x, alpha = alpha, eir = eir, kappa = kappaT, V = V, XX0 = XX)
  return(history)
}


pr2History <- function(x, par) {
  init <- pr2ar_Adam(x, par)
  XX <- init$XX
  alpha <- init$alpha
  eir <- ar2eir(alpha)
  kappa <- kappaF_Adam(XX, par)
  kappaT <- kappa
  for (i in 1:length(alpha)) {
    XX <- PtimesX_Adam(alpha[i], XX, par)
    kappa <- kappaF_Adam(XX, par)
    kappaT <- c(kappaT, kappa)
  }
  V <- eir / kappa
  history <- list(x = x, alpha = alpha, eir = eir, kappa = kappaT, V = V, XX0 = XX)
  return(history)
}


pr2SS_Adam <- function(x, par, tol = 0.001) {
  alpha <- pr2arSS_Adam(x, par)
  XX <- ar2stableWave_Adam(alpha, par)
  xi <- XX2pr29_Adam(XX, par)
  df <- 1
  i <- 1
  while (df > tol) {
    kappa <- kappaF_Adam(XX, par)
    eir <- ar2eir(alpha)
    VC <- eir / kappa
    XX <- PtimesX_Adam(alpha, XX, par)
    xi <- XX2pr29_Adam(XX, par)
    df <- (x - xi)^2
  }
  return(list(XX = XX, kappa = kappa, eir = eir, alpha = alpha, x = x, V = VC))
}


vc2SS_Adam <- function(VC, par, XX0 = NULL) {
  if (is.null(XX0)) {
    XX <- emptyXX_Adam()
    for (i in 1:10) XX <- PtimesX_Adam(.1, XX, par)
  } else {
    XX <- XX0
  }

  for (i in 1:2920) {
    kappa <- kappaF_Adam(XX, par)
    eir <- VC * kappa
    alpha <- eir2ar(eir)
    XX <- PtimesX_Adam(alpha, XX, par)
    x <- XX2pr29_Adam(XX, par)
  }

  df <- 1
  XXl <- XX
  xi <- 1
  while (df > .000001) {
    kappa <- kappaF_Adam(XX, par)
    eir <- VC * kappa
    alpha <- eir2ar(eir)
    XX <- PtimesX_Adam(alpha, XX, par)
    x <- XX2pr29_Adam(XX, par)
    print(c(x = x, xi = xi))
    df <- (x - xi)^2
    xi <- x
  }
  return(list(XX = XX, kappa = kappa, eir = eir, alpha = alpha, x = x, V = VC))
}


#' XXX missing vc2ss implementation.
vcSim_Adam <- function(VC, par, XX0 = NULL) {
  if (is.null(XX0)) {
    XX <- vc2ss(VC[1])
  }
  kT <- eT <- aT <- xT <- c()
  for (i in 1:length(VC)) {
    kappa <- kappaF_Adam(XX, par)
    eir <- VC[i] * kappa
    alpha <- eir2ar(eir)
    XX <- PtimesX_Adam(alpha, XX, par)
    x <- XX2pr29_Adam(XX, par)
    kT <- c(kT, kappa)
    eT <- c(eT, eir)
    xT <- c(xT, x)
    aT <- c(aT, alpha)
  }
  return(list(x = xT, kappa = kT, eir = eT, alpha = aT, XX = XX))
}



#' Make parameters for the Adam model.
#' @param rVals
#' @param nuVals
#' @param rho Fever treatment rate, also below as rhoXstage.
#' @param delta Background rate of protective drugs.
#' @param xiF Fraction of infections entering each stage, as a function of immunity.
#'     Some versions of this function expect one argument while some expect two.
#' @param etaXstage Fever rate in each stage.
#' @param rhoXstage Fever treatment rate in each stage.
#' @param dPIT Cumulative number of days infected.
#' @param iPIT Cumulative incidence of attacks.
#' @param calBfever Burden for fever.
#' @param calBsever Burden for sever malaria.
#' @param calD detection or diagnostic operator for light microscopy.
#' @param calT Tracking of detection of immune relevance of infections.
#' @param pTau Parameter for the vector model.
#' @param qTau Parameter for the vector model.
#' @param N The number of states for active infection.
#' @param tau Days in a time step, which is 10.
#' @return A list with parameters, which are both values and functions.
#' @export
makeParameters_Adam <- function(
      tau = 10,
      N = 9,
      nuVals = c(1 / 2, 1 / 6, 1 / 7, 1 / 8, rep(1 / 9, 3), rep(1 / 12), 0),
      rVals = c(0, 1 / 1000, rep(1 / 250, 6), 1 / 400),
      calDtrue_SoI = c(0, rep(1, 9)),
      calDlm_SoI = c(0, 0.99, .95, 1 - c(3:7)^3 / 8.6^3, .05, .001),
      calDimm_SoI = c(0, 0.95, .9, .6, .3, rep(0, 5)),
      calK_SoI = c(0, 0, c(8:1) / 10),
      calBsevere_SoI = c(0, .3, rep(0, 8)),
      calBfever_SoI = c(0, 1, .7, exp(-c(1:6)), 0),
      rho = 0.2,
      etaXstage = c(0, 0.99, .95, 1 - c(3:7)^3 / 8.6^3, .05, .001)^2,
      rhoXstage = c(0, 1, 1, 1, .8, .8, .6, .6, .4, .4),
      ageForImmune = 8,
      delta = 1 / 730,
      mu1 = 0.001,
      mu2 = 0.01,
      mu3 = 0.001,
      wF = wXageF.0,
      agesF = ages_1066,
      calOF = calO_1066) {
  # Probability of taking drugs, given _any_ state. It's a background rate.
  Dd <- 1 - exp(-delta / tau)
  r <- exp(-rVals * tau)
  PSoI <- calP_SoI(list(N = N, tau = tau, nu = nuVals, r = r))
  SuS <- SuSblock_Adam(list(N = N))
  ages <- agesF()
  ageXimmune <- function(ageInDays) {
    sign(ageInDays > 8 * 365)
  }
  aXi <- ageXimmune(ages * 365)
  demog <- rep(1, 66) * c(diff(ages), 15)
  demog <- demog / sum(demog)
  ix29 <- which(ages >= 2 & ages < 10)
  demog29 <- ages >= 2 & ages < 10
  demog29 <- demog29 / sum(demog29)
  ix5 <- which(ages < 5)
  U5 <- O5 <- demog
  O5[ix5] <- 0
  U5[-ix5] <- 0
  calO <- calOF()
  w <- wF(ages)
  calBfever <- c(calBfever_SoI, calBfever_SoI, 0, 0)
  trt <- rho * rhoXstage * etaXstage
  calBtrtdmal <- c(trt * calBfever_SoI, trt * calBfever_SoI, 0, 0)
  calBuntrmal <- c((1 - trt) * calBfever_SoI, (1 - trt) * calBfever_SoI, 0, 0)
  calBsevere <- c(calBsevere_SoI, calBsevere_SoI, 0, 0)
  calBtrtdsev <- c(trt * calBsevere_SoI, trt * calBsevere_SoI, 0, 0)
  calBuntrsev <- c((1 - trt) * calBsevere_SoI, (1 - trt) * calBsevere_SoI, 0, 0)
  calK <- c(calK_SoI, calK_SoI, 0, 0)
  calK[12] <- calK[3]

  # These are different from functions of the same name above.
  xiF_h <- function(pit, V, p.xi1 = 5, p.xi2 = 1.5, N = 9) {
    xiF_0 <- function(x, p.xi = 2, N = 9, C = 10, D = 0) {
      diff(pbeta(c(0:N) / N, (1 - exp(-p.xi * x)) * C, exp(-p.xi * x) * C + D))
    }
    eps <- 0.001
    if (V > 0) {
      v <- 1 - exp(-(pit[2] + p.xi2 * pit[3]) / V)
    } else {
      v <- 1
    }
    xi1 <- xiF_0(v, p.xi1, N)
    xi <- xi1 + eps
    return(xi / sum(xi))
  }

  list(
    tau = tau,
    N = N,
    calDtrue = c(calDtrue_SoI, rep(1, 10), 0, 0),
    calDlm = c(calDlm_SoI, calDlm_SoI, 0, 0),
    calDimm = c(calDimm_SoI, calDimm_SoI, 0, 0),
    calBfever = calBfever,
    calBtrtdmal = calBtrtdmal,
    calBuntrmal = calBuntrmal,
    calBsevere = calBsevere,
    calBtrtdsev = calBtrtdsev,
    calBuntrsev = calBuntrsev,
    calK = calK,
    nu = nuVals,
    r = r,
    Dd = Dd,
    rho = rho,
    etaXstage = etaXstage,
    rhoXstage = rhoXstage,
    PSoI = PSoI,
    SuS = SuS,
    wF = wF,
    w = w,
    ages = ages,
    ageXimmune = ageXimmune,
    aXi = aXi,
    demog = demog,
    demog29 = demog29,
    U5 = U5,
    O5 = O5,
    calO = calO,
    mu1 = mu1,
    mu2 = mu2,
    mu3 = mu3,
    xiF = xiF_h
  )
}


#' Find integrated force of infection for years given small time steps.
#'
#' @param A the attack rate, equal to 1 - exp(-integrated FOI).
#' @param agg_days The number of days in a year.
#' @param step_days The number of days in a time step.
#' @return A vector of the integrated force of infection for a year.
#' @export
aggregate_foi <- function(A, agg_days, step_days) {
  total_days <- length(A) * step_days
  year_cnt <- total_days %/% agg_days

  if (year_cnt < 1) {
    return(numeric(0))
  }
  aggregated <- numeric(year_cnt)
  for (year_idx in 1:year_cnt) {
    agg_begin <- round((year_idx - 1) * agg_days) + 1
    agg_end <- round(year_idx * agg_days)
    step_begin <- (agg_begin - 1) %/% step_days + 1
    step_end <- (agg_end - 1) %/% step_days + 1

    foi <- 0.0
    for (step_idx in step_begin:step_end) {
      day_begin <- max((step_idx - 1) * step_days + 1, agg_begin)
      day_end <- min(step_idx * step_days, agg_end)
      duration <- day_end - day_begin + 1
      stopifnot(duration <= step_days)
      stopifnot(duration > 0)
      foi <- foi - log(1 - A[step_idx]) * duration / step_days
    }
    aggregated[year_idx] <- foi
  }
  aggregated
}
