#' @title
#' Pilot sampling by SGS
#'
#' @description
#' This function implements the stage 1 subsampling of SAT by SGS method.
#'
#' @usage
#' SAT.stage1.sampling(r1, n, S, Rpar = 0.5)
#'
#' @param r1 pilot subsample size.
#' @param n total sample size.
#' @param S a binary vector of length n. Surrogate observations for all samples.
#' @param Rpar case proportion parameter. The recommended range is (0.3, 0.6), and default is 0.5.
#'
#' @details The region of Rpar that corresponds to lower MSEs is (0.3, 0.5) for case prevalence (i.e., P(Y=1|X)) around 4%.
#'          To avoid failures in the estimation when the case prevalence is low and r1 is small, a slightly
#'          larger Rpar in (0.5, 0.6) can be used without compromising the performance of SAT. Using Rpar=0.5
#'          is a safe choice for most situations.
#'
#' @return The function returns a vector of patient index for whom the manual chart reviews are going to be collected.
#'
#' @references Liu, X., Chubak, J., Hubbard, R. A. & Chen, Y. (2021).
#' SAT: a Surrogate Assisted Two-wave case boosting sampling method, with application to EHR-based association studies.
#'
#' @examples
#' library(SAT)
#' set.seed(0)
#' n <- 1e5
#' beta0  <- c(1/5, 0, 0, 1/2, rep(1/2, 4))
#' d <- length(beta0)
#'
#' X <- rnorm(n*(d-1), -1.5, 1)
#' X <- matrix(X, nrow = n, ncol = d - 1)
#' X <- cbind(1, X)
#'
#' P  <- 1 - 1 / (1 + exp(X %*% beta0))
#' Y  <- rbinom(n, 1, P)
#'
#' a1 <- 0.85 # sensitivity
#' a2 <- 0.95 # specificity
#' pr_s <- vector(mode = "numeric", length = n)
#' pr_s <- a1*(Y==1) + (1-a2)*(Y==0)
#' S <- rbinom(n, 1, pr_s)
#'
#' stage1.index <- SAT.stage1.sampling(r1 = 400, n = 1e5, S, Rpar = 0.5)
#' length(stage1.index)
#'
#' @export


SAT.stage1.sampling <- function(r1, n, S, Rpar){
  pi_s1 <- Rpar*r1/mean(S)/n
  pi_s0 <- (1-Rpar)*r1/(1-mean(S))/n
  PI.prop <- pi_s1*(S==1) + pi_s0*(S==0)
  PI.prop <- PI.prop/sum(PI.prop)
  idx.prop <- sample(1:n, r1, T, PI.prop)
  return(idx.prop)
}
