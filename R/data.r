#' A synthetic lung cancer data.
#'
#' Survival status of patients with lung cancer. The study includes 10,000
#' patients with their survival status, surrogate phenotypes, age at
#' diagnosis, tumor stage and gender. We use it as an illustration
#' example to show the usage of SAT functions. In most situations, the status
#' of patients only available for a small amount of patients, and here we
#' make this observation to be available for all patients to simulate the proposed
#' sampling procedure.
#'
#'
#' @format A data frame with 10,000 rows and 5 columns:
#' \describe{
#'   \item{status}{survival status with 1 = dead and 0 = alive.}
#'   \item{surrogate}{a surrogate of survival status with 1 = dead and 0 = alive.
#'                    The sensitivity is set to be 0.85 and the specificity is
#'                    set to be 0.95.}
#'   \item{age}{standardized patient age at diagnosis.}
#'   \item{tumor stage}{tumor stage with stage I = 1 and stage II = 0.}
#'   \item{gender}{patient gender with female = 1 and male = 0.}
#' }
"lung_cancer"
