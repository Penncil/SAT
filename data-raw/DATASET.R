## code to prepare `DATASET` dataset goes here
set.seed(0)
n <- 1e5
beta0  <- c(-3.2, 0.5, 0.3, 0.5)
d <- length(beta0)

# age
X1 <- sample(c(40:80), n, replace = TRUE)
X1 <- scale(X1)

# sex: female = 1, male = 0
X2 <- rbinom(n, 1, 0.4)

# tumor stage: stage I = 1, stage II = 0
X3 <- rbinom(n, 1, 0.3)

# covariate matrix
X <- cbind(1, X1, X2, X3)

# outcome: status dead = 1, alive = 0
P  <- 1 - 1 / (1 + exp(X %*% beta0))
Y  <- rbinom(n, 1, P)

a1 <- 0.85 # sensitivity
a2 <- 0.95 # specificity
pr_s <- vector(mode = "numeric", length = n)
pr_s <- a1*(Y==1) + (1-a2)*(Y==0)
S <- rbinom(n, 1, pr_s)

lung_cancer <- cbind(Y, S, X[,-1])
colnames(lung_cancer) <- c("status", "surrogate", "age", "sex", "tumor_stage")

usethis::use_data(lung_cancer, overwrite = TRUE)
