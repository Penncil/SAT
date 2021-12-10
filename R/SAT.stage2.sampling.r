#' @title
#' Second stage sampling with SAT.
#'
#' @description
#' This function implements the second stage sampling in SAT by SAT-S or SAT-cY.
#'
#' @usage
#' SAT.stage2.sampling(r1, n, S, Rpar = 0.5, r, stage1.index, stage1.y, X, method = "SAT-S")
#'
#' @param r1 pilot subsample size.
#' @param n total sample size.
#' @param S a binary vector of length n. Surrogate observations for all samples.
#' @param Rpar case proportion parameter. Should be the same as in \code{SAT.stage1.sampling}.
#' @param r second stage subsample size.
#' @param stage1.index a vector of length r1. The output of \code{SAT.stage1.sampling}, i.e., the index of pilot sampled patients.
#' @param stage1.y a binary vector of length r1. The manual chart review results for patients in \code{stage1.index}.
#' @param X a matrix of dimension n times p (the first column needs to be 1). The covariate matrix contains observations for all n samples.
#'
#' @return The function returns a list:
#'   \item{beta.pilot}{the pilot estimator.}
#'   \item{stage1.index}{a vector of index for patients who are selected in pilot sampling.}
#'   \item{stage2.index}{a vector of index for patients who are selected in the second stage sampling.}
#'   \item{stage1.weights}{a vector of weights used in fitting weighted logistic regression for patients who are selected in pilot sampling.}
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
#'                               stage1.index, stage1.y, X, method = "SAT-S")
#' # true phenotype collection
#' stage2.y <-  Y[stage2$stage2.index]
#'
#' stage2$beta.pilot
#'
#' @export

SAT.stage2.sampling <- function(r1, n, S, Rpar, r, stage1.index, stage1.y, X, method = "SAT-S"){

  # note: idx.prop could be a bit different than the output of SAT.stage1.sampling()
  idx.prop <- stage1.index
  y.prop <- stage1.y
  x.prop <- X[idx.prop, ]
  s.prop <- S[idx.prop]


  pi_s1 <- Rpar*r1/mean(S)/n
  pi_s0 <- (1-Rpar)*r1/(1-mean(S))/n
  PI.prop <- pi_s1*(S==1) + pi_s0*(S==0)
  PI.prop <- PI.prop/sum(PI.prop)
  pinv.prop <- n
  pinv.prop <- 1/PI.prop[idx.prop]


  # weighted estimation to get pilot estimator
  fit.prop <- getMLE(x = x.prop, y = y.prop, w = pinv.prop)
  beta.prop <- fit.prop$par
  if (is.na(beta.prop[1])) {
    stop(paste("pilot Estimation:", fit.prop$msg))
  }
  P.prop <- 1 - 1/(1 + exp(X %*% beta.prop))


  # use pilot to compute n*r1*M_x
  p.prop <- P.prop[idx.prop]
  w.prop <- p.prop * (1 - p.prop)
  W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))


  if (method == "SAT-S"){

    #--------SAT-S--------#

    # compute SSP with S
    PI.mMSE.S <- sqrt((S - P.prop)^2 * rowSums((X %*% W.prop)^2))
    PI.mMSE.S <- PI.mMSE.S/sum(PI.mMSE.S)
    idx.mMSE <- sample(1:n, r, T, PI.mMSE.S)


  } else if (method == "SAT-cY") {

    # compute P(Yi=1|Si=1) and P(Yi=1|Si=0)
    a1.hat <- sum((s.prop[y.prop==1]==1))/sum(y.prop==1)
    condY <- pmin(a1.hat*P.prop/mean(S),1)
    condY[S==0] <- pmin((1-a1.hat)*P.prop[S==0]/(1-mean(S)),1)

    #--------SAT-cY--------#

    # compute SSP with E(Y|S)
    PI.mMSE.cY <- sqrt((condY - P.prop)^2 * rowSums((X %*% W.prop)^2))
    PI.mMSE.cY <- PI.mMSE.cY/sum(PI.mMSE.cY)
    idx.mMSE <- sample(1:n, r, T, PI.mMSE.cY)

  }

  return(list(beta.pilot = beta.prop, stage1.index = idx.prop, stage2.index = idx.mMSE, stage1.weights = pinv.prop))
}
