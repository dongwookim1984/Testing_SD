# Testing_SD
# R - adaptive test of stochastic dominance 
Authors : Dongwoo Kim 

This project provides an R function for a nonparametric test of stochastic monotonicity proposed by Chetverikov, Wilhelm and Kim (2018). The new test adapts to the unknown smoothness of the conditional distribution of interest and possesses desirable asymptotic properties. The test  asymptotically controls size at a polynomial rate, is non-conservative, and detects local alternatives that converge to the null with the fastest possible rate. The test is based on a data-driven bandwidth value and the critical value for the test takes this randomness into account. Simulations indicate that the test performs well in finite samples and may be significantly more powerful than existing alternative procedures.

Files contained in this package:

- The file `Test_SD_functions.R` contains necessary functions for the test. The function "CKW_test" provides test results including the test statistic, the bootstrap critical value, and the value of bandwidth at which the test statistic is maximised.
- The file `Test_SD_example.R` contains examples which show how to use the function on simulated data. 


## Installation
Test_SD_functions.R file needs to be in the working directory to use. The function requires the R package "foreach". The user can install this package by typing

```
install.packages("foreach")
```



## Syntax
The function CKW.test provides the test of the null of stochastic monotonicity between two continuous random variables Y and X, both supported on [0, 1], 

```
H_0 : F_{Y|X}(y|x') >= F_{Y|X}(y|x'') for all y, x', x'' in (0, 1) with x' <= x''
```

against the alternative

```
H_a : F_{Y|X}(y|x') < F_{Y|X}(y|x'') for some y, x', x'' in (0, 1) with x' <= x''
```

Syntax:

```
CKW.test(Y, X, brep, num.grid, u, delta, alpha)
```

where
- `Y` and `X` are random variables that need to be tested.
- `brep` is the number of bootstrap replications.
- `num.grid` is the number of grid points on the empirical support of Y. The paper suggests the sample size the user uses. For computational purposes, the user can specify the fewer number of grid points to yield results more quickly. 
- `u` and `delta` are tuning parameters which control the number of bandwidth values for the test statistic.  The default is set such that u = 2/3 and delta = 1/2. 
- `alpha` is the confidence level. The default is 0.95.


## Output

The function CKW.test generates a column vector containing four values. 

- The first element indicates whether the test rejects the null. 1 if the null is rejected and 0 otherwise.
- The second element shows the test statistic.
- The third element shows the critical value computed by multiplier bootstrap.
- The fourth element shows the bandwidth value at which the test statistic is maximised.


## Examples

Test with default options (brep = 200, num.grid = 100, u = 2/3, delta = 1/2, alpha = 0.95):
```
res <- CKW.test(Y, X)
```

Test with options (brep = 500, num.grid = sample size, u = 1/2, delta = 1/2, alpha = 0.95):
```
res <- CKW.test(Y, X, brep = 500, num.grid = length(Y), u = 2/3, delta = 1/2, alpha = 0.95)
```

# Reference
Chetverikov, D. Wilhelm, D. and Kim, D. (2018), "An adaptive test of stochastic monotonicity"
