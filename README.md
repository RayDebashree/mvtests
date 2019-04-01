### Description
`mvtests` will be a suite of functions implementing cross-phenotype association tests. Currently, it only implements the POM-LRT approach that tests association of one or more traits with a single genetic marker using a proportional odds regression model of genotype on traits. It uses individual-level phenotype genotype data on unrelated individuals. The R function `pom()` implements this association test. Please refer/cite:

Ray, D. and Chatterjee, N. "Effect of Non-Normality and Low Count Variants on Cross-Phenotype Association Tests in GWAS". *In revision*, 2019.

**Key Words:** Cross-phenotype association; GWAS; Joint modeling; Multiple traits; Multivariate analysis; Non-normal traits; Proportional Odds Model; 


### Requirements
R (>= 3.0.1), MASS, lmtest

### Changes
Version 0.2 - April 1, 2019
> First public release of the software.


### Notes
1. The proportional odds framework of POM allows one or more phenotypes, which may be binary and/or continuous. POM based on likelihood ratio test (LRT) has been previously implemented by [MultiPhen](https://rdrr.io/cran/MultiPhen/).
    * Caution: As of 1-Apr-2019, MultiPhen has been archived by CRAN in May 2018.
    * We did not find an option to modulate parameters in MultiPhen that control optimization of the likelihood function. We faced frequent non-convergence of the optimization routine for rare genetic variants. Our `pom()` function is flexible in that aspect. The likelihood maximization in `pom()` is performed using `R` base function `optim()` and the user can conveniently
control parameters related to the optimization.
    * Unlike MultiPhen, if there is non-convergence of the likelihood optimization, `NA` is returned as output by `pom()`.
    * When the likelihood optimization converged, POM-LRT outputs using `pom()` with `test.method="LRT"` and using `mPhen()` from MultiPhen are identical.
    * Unlike MultiPhen, `pom()` additionally implements the Wald test (although we found robust performance of POM-LRT under different situations of non-normality compared to POM-Wald)
 
2. The method POM and its software is designed for unrelated individuals. If you have two cohorts with overlapping samples and you want to analyze the combined sample, it is desirable to exclude the overlapping individuals, and any related individuals. 

3. The genotype X takes values 0, 1 or 2. For a given sample, if X has only two possible values, a logistic regression (`glm`) is used instead of proportional odds regression (`polr`).

4. Any individual with at least one missing observation in the data is removed.  
