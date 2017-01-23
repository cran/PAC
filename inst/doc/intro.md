---
title: "Using the PAC package"
date: "2017-01-22"
output: rmarkdown::html_vignette
bibliography: references.bib
vignette: >
  %\VignetteIndexEntry{Using the PAC package}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


This vignette illustrates the basic usage of the PAC package for R. 

Synthetic data
--------------
For simplicity, we will generate the synthetic data from a mixture of gaussian. The details of the construction are unimportant, but are provided below for completeness.


```r
# Problem parameters
n = 5e3                       # number of observations
p = 1                         # number of dimensions
K = 3                         # number of clusters
w = rep(1,K)/K                # component weights
mu <- c(0,2,4)                # component means
sd <- rep(1,K)/K              # component standard deviations

# Generate K mixture of gaussian 
g <- sample(1:K,prob=w,size=n,replace=TRUE)   # ground truth for clustering
X <- as.matrix(rnorm(n=n,mean=mu[g],sd=sd[g]))
```



# Basic usage

The main entry point for the PAC package is the function `PAC`. To begin, we call `PAC` with all of the default options.

```r
library(PAC)
y <- PAC(X, K)
```

```
## Input Data: 5000 by 1
## Partition method: Discrepancy based partition
## Maximum level: 40
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

We can evaluate the accuracy of the clustering result by computing the fmeasure between it and the ground truth:
of the default options.

```r
print(fmeasure(g,y))
```

```
## [1] 0.9940302
```

# Changing the default parameters

We now show how to change the default value of PAC's has two optional arguments: the maximum level of the partition("maxlevel") and the partition method("method"). 

For example, one can use a larger number of partitions (the default is 40):

```r
y.deep <- PAC(X, K, maxlevel = 100)
```

```
## Input Data: 5000 by 1
## Partition method: Discrepancy based partition
## Maximum level: 100
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

```r
print(fmeasure(g,y.deep))
```

```
## [1] 0.9801951
```

Or, one can change the partition method to be bsp with limited lookahead See [@lu2013multivariate], [@wong2010optional]. The default is dsp [@yang2014density]. 

```r
y.ll <- PAC(X, K, method = 'll')
```

```
## Input Data: 5000 by 1
## Partition method: BSP with limited-lookahead
## Maximum level: 40
## partition completed
## [1] "Initial Clustering..."
## [1] "Merging..."
```

```r
print(fmeasure(g,y.ll))
```

```
## [1] 0.9956008
```
# References
