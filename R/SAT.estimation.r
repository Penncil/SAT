#' @title
#' SAT estimation based on pooled subsample
#'
#' @description
#' This function implements the final SAT estimation based on data pooled from two stages of sampling.
#' A weighted logistic regression is conducted.
#'
#' @usage
#' SAT.estimation(S, X, beta.pilot, stage1.index, stage2.index, stage1.weights,
#'                stage1.y, stage2.y, method = "SAT-S")
#'
#'
#' @param S a binary vector of length n. Surrogate observations for all samples.
#' @param X a matrix of dimension n times p (the first column needs to be 1). The covariate matrix contains observations for all n samples.
#' @param beta.pilot the pilot estimator.
#' @param stage1.index a vector of length r1. The index of pilot sampled patients.
#' @param stage2.index a vector of length r. The index of second-stage sampled patients.
#' @param stage1.weights a vector of weights for patients who are selected in pilot sampling.
#' @param stage1.y a binary vector of length r1. The manual chart review results for patients in \code{stage1.index}.
#' @param stage2.y a binary vector of length r. The manual chart review results for patients in \code{stage2.index}.
#' @param method two methods are available: SAT-S or SAT-cY.
#'
#'
#'
#' @return The function returns the final SAT estimates of the association coefficients.
#'
#'
#'
#' @references Liu, X., Chubak, J., Hubbard, R. A. & Chen, Y. (2021).
#' SAT: a Surrogate Assisted Two-wave case boosting sampling method, with application to EHR-based association studies.
#'
#' @examples
#' library(SAT)
#' set.seed(0)
#'
#' colnames(lung_cancer)
#' X <- cbind(1, lung_cancer[,3:5])
#' Y <- lung_cancer[,1]
#' S <- lung_cancer[,2]
#'
#' # pilot sampling
#' stage1.index <- SAT.stage1.sampling(r1 = 400, n = 1e5, S, Rpar = 0.5)
#' # true phenotype collection
#' stage1.y <- Y[stage1.index]
#'
#' # second stage sampling
#' stage2 <- SAT.stage2.sampling(r1 = 400, n = 1e5, S, Rpar = 0.5, r = 800,
#'                               stage1.index, stage1.y, X, method = "SAT-cY")
#' # true phenotype collection
#' stage2.y <-  Y[stage2$stage2.index]
#'
#' # final estimation
#' SAT.est <- SAT.estimation(S, X, beta.pilot = stage2$beta.pilot, stage1.index = stage1.index,
#'                stage2.index = stage2$stage2.index,
#'                stage1.weights = stage2$stage1.weights,
#'                stage1.y = stage1.y, stage2.y = stage2.y,
#'                method = "SAT-cY")
#' SAT.est$SAT.estimate
#'
#' @export


SAT.estimation <- function(S, X, beta.pilot, stage1.index, stage2.index, stage1.weights, stage1.y, stage2.y, method = "SAT-S"){

  y.prop <- stage1.y
  y.mMSE <- stage2.y
  beta.prop <- beta.pilot
  idx.prop <- stage1.index
  idx.mMSE <- stage2.index
  pinv.prop <- stage1.weights
  x.prop <- X[idx.prop, ]
  s.prop <- S[idx.prop]


  P.prop <- 1 - 1/(1 + exp(X %*% beta.prop))


  # use pilot to compute n*r1*M_x
  p.prop <- P.prop[idx.prop]
  w.prop <- p.prop * (1 - p.prop)
  W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))


  x.mMSE <- X[c(idx.mMSE, idx.prop), ]
  y.mMSE <- c(y.mMSE, y.prop)



  if (method == "SAT-S"){

    # compute SSP with S
    PI.mMSE.S <- sqrt((S - P.prop)^2 * rowSums((X %*% W.prop)^2))
    PI.mMSE.S <- PI.mMSE.S/sum(PI.mMSE.S)

    # get combined weights
    pinv.mMSE <- c(1/PI.mMSE.S[idx.mMSE], pinv.prop)
    fit.mMSE <- getMLE(x = x.mMSE, y = y.mMSE, w = pinv.mMSE)
    if (fit.mMSE$message == "Successful convergence") {
      beta.mMSE <- fit.mMSE$par
    }else{
      stop(paste(rr, "SAT-S Estimation:", fit.mMSE$message))
    }

  } else if (method == "SAT-cY") {

    # compute P(Yi=1|Si=1) and P(Yi=1|Si=0)
    a1.hat <- sum((s.prop[y.prop==1]==1))/sum(y.prop==1)
    condY <- pmin(a1.hat*P.prop/mean(S),1)
    condY[S==0] <- pmin((1-a1.hat)*P.prop[S==0]/(1-mean(S)),1)

    # compute SSP with E(Y|S)
    PI.mMSE.cY <- sqrt((condY - P.prop)^2 * rowSums((X %*% W.prop)^2))
    PI.mMSE.cY <- PI.mMSE.cY/sum(PI.mMSE.cY)

    # get combined weights
    pinv.mMSE <- c(1/PI.mMSE.cY[idx.mMSE], pinv.prop)
    fit.mMSE <- getMLE(x = x.mMSE, y = y.mMSE, w = pinv.mMSE)
    if (fit.mMSE$message == "Successful convergence") {
      beta.mMSE <- fit.mMSE$par
    }else{
      stop(paste(rr, "SAT-cY Estimation:", fit.mMSE$message))
    }
  }
  return(list(SAT.estimate = beta.mMSE))
}
