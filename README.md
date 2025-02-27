The glfBLUP (genetic latent factor BLUP) package facilitates compression of secondary features in a genomic prediction context to a smaller number of latent genetic factors.

### Installation
```
devtools::install_github("KillianMelsen/glfBLUP")
```

### Example
Assume that the dataframe `d` contains rows representing plots and columns representing the genotype ID `G`, secondary features, and focal trait `Y` in that order.
We first separate training and test data, assuming that the focal trait data in `d` has been set to NA:
```R
d.train <- droplevels(na.omit(d))
sec <- names(d[2:(ncol(d) - 1)])
foc <- names(d)[ncol(d)]
```
We then redundancy filter the data to get rid of pairs of extremely highly correlated features:
```R
temp <- glfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = F)
d.train.RF <- cbind(temp$data.RF, d.train[foc])
d.RF <- d[names(d.train.RF)]
sec.RF <- names(d.RF[2:(ncol(d.RF) - 1)])
```
Next is regularization of the genetic and residual covariance matrices of the secondary features:
```R
folds <- glfBLUP::createFolds(genos = unique(as.character(d.train.RF$G)))
tempG <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "genetic", dopar = TRUE)
tempE <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "residual", dopar = TRUE)
Rg.RF.reg <- tempG$optCor
```
We can then fit the factor model using the regularized genetic correlation matrix:
```R
# data is only used to determine the sample size for the MP-bound. what specifies that it's a genetic correlation matrix, so the
# number of training genotypes should be used, and not the number of training individuals (= genotypes * replicates).
FM.fit <- glfBLUP::factorModel(data = d.train.RF[c("G", sec.RF)], cormat = Rg.RF.reg, what = "genetic")
```
Then we estimate factor scores:
```R
# Loadings and uniquenesses were estimated on the correlation scale, but should be on the covariance scale for
# genetic-thomson scores:
D <- sqrt(diag(tempG$Sg)) # Getting standard deviations
L.cov <- diag(D) %*% FM.fit$loadings
PSI.cov <- outer(D, D) * FM.fit$uniquenesses
      
# CV1 Factor scores:
CV1.d.RF <- d.RF
CV1.d.RF[which(is.na(CV1.d.RF$Y)), 2:ncol(CV1.d.RF)] <- NA
CV1.F.scores <- glfBLUP::factorScores(data = CV1.d.RF[c("G", sec.RF)],
                                      loadings = L.cov,
                                      uniquenesses = PSI.cov,
                                      m = FM.fit$m,
                                      type = "genetic-thomson-repdiv",
                                      Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
      
CV1.d.final <- cbind(CV1.F.scores, CV1.d.RF$Y)
names(CV1.d.final)[ncol(CV1.d.final)] <- "Y"
names(CV1.d.final)[1] <- "G"
      
# CV2 Factor scores:
# First recenter/rescale the training and test secondary data together:
d.RF[sec.RF] <- sapply(d.RF[sec.RF], scale)
CV2.F.scores <- glfBLUP::factorScores(data = d.RF[c("G", sec.RF)],
                                      loadings = L.cov,
                                      uniquenesses = PSI.cov,
                                      m = FM.fit$m,
                                      type = "genetic-thomson-repdiv",
                                      Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
      
CV2.d.final <- cbind(CV2.F.scores, d.RF$Y)
names(CV2.d.final)[ncol(CV2.d.final)] <- "Y"
names(CV2.d.final)[1] <- "G"
```
And finally we can select the relevant factors and calculate the CV1 and CV2 multivariate genomic predictions using the kinship matrix `K`:
```R
selection <- glfBLUP::factorSelect(CV1.d.final, procedure = "leaps", verbose = F)
CV1.temp <- glfBLUP::glfBLUP(data = CV1.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
CV2.temp <- glfBLUP::glfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
```
