# SAT: a Surrogate Assisted Two-wave case boosting sampling method


## Outline

1. Description 
2. SAT workflow
3. Package installation and example

## Description

This README is prepared for the package SAT of implementing the surrogate assisted two-wave (SAT) case boosting sampling method. SAT is proposed to simultaneously solve two commonly co-existing issues in Electronic health records (EHR) based association studies, i.e., the bias in coefficient estimation brought by potentially error-prone EHR-derived phenotypes (i.e., surrogates) and the efficiency loss caused by low prevalence phenotypes (e.g., rare disease). Target on improving estimation accuracy, SAT selects a subsample for outcome validation through manual chart review subject to budget constraints, and a model is then fitted based on the subsample only.



## SAT Workflow
![](Figure0.png)

At stage 1, for enriching cases we apply surrogate-guided sampling (SGS) to obtain a pilot subsample of size n1, for which we collect the true phenotypes and use them to obtain a pilot estimator for the association parameters. At stage 2, based on the pilot estimator, the covariates, and the surrogate phenotype, we compute the SAT sampling probabilities which target controlling MSE, and then apply them to obtain the second-stage subsample of size n2, for which the true phenotypes are collected. The final estimator is obtained by fitting a weighted logistic regression on the pooled subsample of size n1 + n2 with their true phenotypes.

## Package installation and example
What we need:
1. A dataset including risk factors and surrogates for the full cohort;
2. Download and install R package [SAT](https://github.com/Penncil/SAT).

The example we used here is a synthetic dataset of lung cancer. The true phenotype that is of our interest is the survival status of patients with lung cancer. The study includes 10,000 patients with their survival status, surrogate phenotypes, age at diagnosis, tumor stage and gender. Here we note that in practice the true phenotype $Y$ (i.e., the survival status of patients in this lung cancer example) of the full cohort is not available, and here we make these observations be available for all patients to simulate the proposed sampling procedure.

Step 0: Install the R package SAT

```{r}
# install via github:
install.packages("devtools")
library(devtools)
devtools::install_github("penncil/SAT")
library(SAT)
```
With the installed SAT package, next we use SAT functions to select "informative" samples for which we collect true phenotypes by manual chart review. 

Step 1: 1st stage subsampling
```{r}
set.seed(0)

colnames(lung_cancer)
X <- cbind(1, lung_cancer[,3:5])
Y <- lung_cancer[,1]
S <- lung_cancer[,2]

# 1st stage sampling to get pilot sample
stage1.index <- SAT.stage1.sampling(r1 = 400, n = 1e5, S, Rpar = 0.5)

# true phenotype collection
# Note: in practice stage1.y are collected by manual chart review
stage1.y <- Y[stage1.index]
```


Stage 2: 2nd stage subsampling and final model fitting 
```{r}
# second stage sampling
stage2 <- SAT.stage2.sampling(r1 = 400, n = 1e5, S, Rpar = 0.5, r = 800,
                              stage1.index, stage1.y, X, method = "SAT-cY")
# true phenotype collection
# Note: in practice stage2.y are collected by manual chart review
stage2.y <-  Y[stage2$stage2.index]      


# final model fitting with combined samples
SAT.est <- SAT.estimation(S, X, beta.pilot = stage2$beta.pilot, 
               stage1.index = stage1.index,
               stage2.index = stage2$stage2.index,
               stage1.weights = stage2$stage1.weights,
               stage1.y = stage1.y, stage2.y = stage2.y,
               method = "SAT-cY")
SAT.est$SAT.estimate
```

