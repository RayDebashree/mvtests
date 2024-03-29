### Description
`mvtests` will be a suite of functions implementing cross-phenotype association tests. Currently, it implements the POM-LRT approach that tests association of one or more traits with a single genetic marker using a proportional odds regression model of genotype on traits. It uses individual-level phenotype genotype data on unrelated individuals. The R function `pom()` implements this association test. It also includes functions to implement Nyholt-Šidák correction for multiple tests. Please refer/cite:

Ray, D. and Chatterjee, N. "Effect of Non-Normality and Low Count Variants on Cross-Phenotype Association Tests in GWAS". *European Journal Human Genetics*, 28(3):300-312, 2020. https://www.nature.com/articles/s41431-019-0514-2

**Key Words:** Cross-phenotype association; GWAS; Joint modeling; Multiple traits; Multivariate analysis; Non-normal traits; Proportional Odds Model 


### Requirements
R (>= 3.0.1), MASS, lmtest

### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/mvtests/blob/master/mvtests_v0.3.R?raw=TRUE")
```
It is recommended to download/copy the stand-alone R program in this repository, save it in your local directory of choice and `source()` it from your local directory. When a new version of the software is available, older versions may be removed from this repository, and the above `devtools::source_url()` technique may not work.


### Changes
Version 0.3 - June 9, 2020
> An updated version of the software to include Nyholt-Šidák correction for multiple tests.

Version 0.2 - April 1, 2019
> First public release of the software.


### Notes
1. The proportional odds framework of POM allows one or more phenotypes, which may be binary and/or continuous. Likelihood ratio test (LRT) based on POM has been previously implemented by [MultiPhen](https://rdrr.io/cran/MultiPhen/), which has been archived by CRAN in May 2018 (as of 1-Apr-2019). There are some key advantages of our software:
    * We did not find an option to modulate parameters in MultiPhen that control optimization of the likelihood function. We faced frequent non-convergence of the optimization routine for rare genetic variants. Our `pom()` function is flexible in that aspect. The likelihood maximization in `pom()` is performed using `R` base function `optim()` and the user can conveniently
control parameters related to the optimization.
    * Unlike MultiPhen, if there is non-convergence of the likelihood optimization, `NA` is returned as output by `pom()` apart from error messages.
    * When the likelihood optimization converged, POM-LRT outputs using `pom()` with `test.method="LRT"` and using `mPhen()` from MultiPhen are identical.
    * Unlike MultiPhen, `pom()` additionally implements the Wald test (although we found robust performance of POM-LRT for rare or low-frequency variants compared to POM-Wald)
 
2. The method POM and its software is designed for unrelated individuals. If you have two cohorts with overlapping samples and you want to analyze the combined sample, it is desirable to exclude the overlapping individuals, and any related individuals. 

3. The genotype X takes values 0, 1 or 2. For a given sample, if X has only two possible values, a logistic regression (`glm`) is used instead of proportional odds regression (`polr`).

4. Any individual with at least one missing observation in the data is removed before association test.
