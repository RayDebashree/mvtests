###--------------------- A suite of codes for testing association of multiple phenotypes with a single SNP -------------------------
#
# Cite: "Effect of Non-Normality and Low Count Variants on Cross-Phenotype Association Tests in GWAS".
#
##--------------------------------------------- Version 0.3 (dated June 9, 2020) --------------------------------------------------
# Corresponding Author: Debashree Ray <dray@jhu.edu>

############################################
library(MASS)
library(lmtest)

message("===================================================")
message("	       mvtests v0.3 is loaded")
message("===================================================")
message("If you use this software, please cite:")
message("Ray et al.(2020) Effect of Non-Normality and Low Count Variants")
message("    on Cross-Phenotype Association Tests in GWAS.")
message("    European Journal of Human Genetics, 28(3):300-312.") 
message("---------------------------------------------------")


############################################
#---------------- function for format and dimension checks of inputs in pom
.format.check<-function(Y, X, COV=NULL, msg.mute=FALSE)
{
   # check formats
   if(class(Y)!="data.frame" | class(X)!="data.frame")
       	stop("Inputs Y(phenotype), X(SNP genotype) must be in data frame format.")

   # check dimensions
   n<-nrow(Y)
   if(nrow(X)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.")
   if(!is.null(COV)) { if(nrow(COV)!=n) stop("Sample sizes (no. of rows) in input data frames do not match.") }
   if(ncol(X)!=1) stop("Data frame X (SNP genotype) should have a single column.")
}

#---------------- function for identifying number of parameters for which starting values are required in optim in polr()
# this function may come in handy when there is any error in optim and starting values need to be changed/user-specified
.get.start.length<-function(formula, design)
{
    m <- model.frame(formula, model.frame(design), na.action = na.pass)[,1] # first column is response (factor) in this model matrix
    start.length <- length(attr(terms(formula), "term.labels")) + (length(levels(m))-1)
    return(start.length)
}

.get.start.length.polr<-function(formula, X)
{
   start.length <- length(attr(terms(formula), "term.labels")) + (length(unique(X))-1)
   return(start.length)
}

#---------------- function to grab warnings from running functions (courtesy: https://stackoverflow.com/questions/3903157/)
grabWarnings <- function(expr) {
    myWarnings <- NULL
    wHandler <- function(w) {
        myWarnings <<- c(myWarnings, list(w))
        invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(output = val, warnings = myWarnings)
}

#--------------------------------------- Main function for POM-LRT (uses polr) --------------------------------
pom<-function(Y, X, COV=NULL, test.method="LRT", msg.mute=FALSE, no.format.check=FALSE, ...)
{
   #--------------------------- CHECKS ----------------------------
   if(!no.format.check) .format.check(Y, X, COV, msg.mute)
   if(!(test.method=="LRT"|test.method=="Wald")) stop("test.method is either 'LRT' (default) or 'Wald'.")

   # check names
   n<-nrow(Y)   # sample size
   q<-0         # no. of covariates
   dataf<-cbind(Y,X)
        Yname<-colnames(Y) ;    Xname<-colnames(X) ;    COVname<-NULL
   if(!is.null(COV)) { dataf<-cbind(dataf, COV) ; COVname<-colnames(COV) ; q<-ncol(COV) }
   if( length(unique(colnames(dataf)))!=ncol(dataf) )
        stop("One or more data frame inputs have common column names. Please provide distinct column names across all input data frames.")

   # other checks
   dataf<-na.omit(dataf)
   nobs<-nrow(dataf)
   if(nobs!=n & !isTRUE(msg.mute)) message('Removing samples with missing observations...')
   if(nobs==0) stop("No observation left after removing missing observations.")
   k<-ncol(Y)
   # distribution of genotypes in the complete non-missing data
   tabx <- table(dataf[,k+1])
   if(length(unique(dataf[,k+1]))==1)
     stop('All individuals have the same genotype!')
   if(length(unique(dataf[,k+1]))>2) {
     if(!is.factor(dataf[,k+1])) dataf[,k+1]<-as.factor(dataf[,k+1])
   }else {
     if(is.factor(dataf[,k+1])) dataf[,k+1]<-as.numeric(as.character(dataf[,k+1]))
     # glm() can take values 0 or 1 for the response; need to ensure that
     if(!(identical(unique(dataf[,k+1]),c(0,1)) | identical(unique(dataf[,k+1]),c(1,0)))){
	xu<-sort(unique(dataf[,k+1]))
	dataf[which(dataf[,k+1]==xu[1]),k+1]<-0
	dataf[which(dataf[,k+1]==xu[2]),k+1]<-1
     }
   }

   #------------------------------ FORMULAE -------------------------
     preds.cov<-paste(c(Yname,COVname),collapse="+")
     formula.cov<-as.formula(paste( Xname, "~", preds.cov, sep="" ))
	# for null model
	if(!is.null(COV)) {
	   preds.cov0<-paste(c(COVname),collapse="+")
	   formula.cov0<-as.formula(paste( Xname, "~", preds.cov0, sep="" ))
	}else {
	   formula.cov0<-as.formula(paste( Xname, "~ 1", sep="" ))
	}

   #------------------------------MODEL FIT-------------------------
   if(length(unique(dataf[,k+1]))==2)
   {
     warning('Only 2 possible values for genotype X; fitting logistic model instead of proportional odds model.')
        fit<-try(glm(formula.cov, family="binomial", data=dataf, ...), silent=TRUE)
	  if(!inherits(fit, "try-error")) fit.conv<-fit$"converged" else fit.conv<-FALSE
        fit0<-try(glm(formula.cov0, family="binomial", data=dataf, ...), silent=TRUE)
	  if(!inherits(fit0, "try-error")) fit0.conv<-fit0$"converged" else fit0.conv<-FALSE
   }else
   {
        fit<-try(polr(formula.cov, method="logistic", data=dataf, Hess=TRUE, ...), silent=TRUE)
	  if(!inherits(fit, "try-error")) fit.conv<-(fit$"convergence"==0) else fit.conv<-FALSE
        fit0<-try(polr(formula.cov0, method="logistic", data=dataf, Hess=TRUE, ...), silent=TRUE)
	  if(!inherits(fit0, "try-error")) fit0.conv<-(fit0$"convergence"==0) else fit0.conv<-FALSE
   }

   # Output
   if(inherits(fit, "try-error") | !fit.conv){
        # Get length of 'start' (i.e., no. of parameters to estimate)
        start.length<-.get.start.length.polr(formula.cov, dataf[,k+1])
     if(inherits(fit, "try-error")){
        error.msg<-paste(fit[1],"-- For",Xname,"may need to change full model starting values parameter 'start' of length",start.length)     
     }else{
	error.msg<-"Algorithm did not converge for full model. NA's returned."
     }
     if(!isTRUE(msg.mute)) message(error.msg)
     betas <- se.betas <- rep(NA, start.length)
     if(test.method=="Wald"){
       stat.w <- df.w <- pval.w <- statwaldF <- pvalwaldF <- NA
     }else{
       stat <- df <- pval <- NA
     }
   }else if(inherits(fit0, "try-error") | !fit0.conv){
        # Get length of 'start' (i.e., no. of parameters to estimate)
        start.length<-.get.start.length.polr(formula.cov0, dataf[,k+1])
     if(inherits(fit0, "try-error")){
        error.msg<-paste(fit0[1],"-- For",Xname,"may need to change null model starting values parameter 'start' of length",start.length)
     }else{
	error.msg<-"Algorithm did not converge for null model. NA's returned where null model is used."
     }
     if(!isTRUE(msg.mute)) message(error.msg)
     betas <- summary(fit)$"coefficients"[,1]
     se.betas <- summary(fit)$"coefficients"[,2]
     if(test.method=="Wald"){
       #------------------------------TESTING (Only Wald; Error fitting null model)-------------------------
       # Wald Test using chisq dist
       if("polr" %in% class(fit)){
          beta<-matrix(coef(fit),ncol=1)
 	  Sig<-vcov(fit)[1:k,1:k]
       }else{
          beta<-matrix(coef(fit)[-1],ncol=1)
	  Sig<-vcov(fit)[-1,-1]
       }
       Sig.inv <- try(solve(Sig), silent=TRUE)
       if(inherits(Sig.inv, "try-error")){
	  stat.w <- df.w <- pval.w <- NA
       	  if(exists("error.msg")) error.msg <- paste(error.msg,";",Sig.inv[1]) else error.msg <- Sig.inv[1]
       }else{
	  stat.w <- drop(t(beta)%*%solve(Sig)%*%beta)
	  df.w <- k
     	  pval.w <- pchisq(stat.w, df=df.w, lower.tail=F)
       }
       # Wald Test with F dist (using lmtest package) (Not possible since it requires null model fit)
       statwaldF <- pvalwaldF <- NA
       error.msgF<-"Wald Test using F distribution uses lmtest package - cannot be implemented due to error fitting null model."
       if(exists("error.msg")) error.msg <- paste(error.msg,";",error.msgF) else error.msg <- error.msgF
     }else{
       stat <- df <- pval <- NA
     }     
   }else{
     betas <- summary(fit)$"coefficients"[,1]
     se.betas <- summary(fit)$"coefficients"[,2]
     #------------------------------TESTING (Both Wald & LRT)-------------------------
     if(test.method=="Wald"){
       # Wald Test using chisq dist
       if("polr" %in% class(fit)){
          beta<-matrix(coef(fit),ncol=1)
	  Sig<-vcov(fit)[1:k,1:k]
       }else{
          beta<-matrix(coef(fit)[-1],ncol=1)
	  Sig<-vcov(fit)[-1,-1]
       }
       Sig.inv <- try(solve(Sig), silent=TRUE)
       if(inherits(Sig.inv, "try-error")){
          stat.w <- df.w <- pval.w <- NA
          if(exists("error.msg")) error.msg <- paste(error.msg,";",Sig.inv[1]) else error.msg <- Sig.inv[1]
       }else{
	  stat.w <- drop(t(beta)%*%solve(Sig)%*%beta)
	  df.w <- k
          pval.w <- pchisq(stat.w, df=df.w, lower.tail=F)
       }
       # Wald Test using F dist (using lmtest package)
       outwaldF <- try(waldtest(fit, fit0, test="F"), silent=TRUE)
       if(inherits(outwaldF, "try-error")){
          statwaldF <- pvalwaldF <- NA
          if(exists("error.msg")) error.msg <- paste(error.msg,";",outwaldF[1]) else error.msg <- outwaldF[1]
       }else{
	  statwaldF <- outwaldF$F[2]
	  pvalwaldF <- outwaldF$"Pr(>F)"[2]
       }
     }else{
       # LRT
       stat<-c(2*(logLik(fit)-logLik(fit0)))
       if(stat<0){
         stat <- df <- pval <- NA
	 uns.msg <- "Unstable results for LRT - NAs returned."
         if(exists("error.msg")) error.msg <- paste(error.msg,";",uns.msg) else error.msg <- uns.msg
       }else{
         if("polr" %in% class(fit)){
	   df<-fit$"edf"-fit0$"edf"
         }else df<-fit$"rank"-fit0$"rank"
         pval <- pchisq(q=stat, df=df, ncp=0, lower.tail=FALSE)
       }
     }
     if(!exists("error.msg")) error.msg <- "OK"
   }
   if(test.method=="Wald"){
     return(list(coef=betas, SE.coef=se.betas, stat.wald.chisq=stat.w, df.wald.chisq=df.w, pvalue.wald.chisq=pval.w, stat.wald.F=statwaldF, pvalue.wald.F=pvalwaldF, n.obs=nobs, geno.dist=tabx, error.msg=error.msg))
   }else{
     return(list(coef=betas, SE.coef=se.betas, stat.lrt=stat, df.lrt=df, pvalue.lrt=pval, n.obs=nobs, geno.dist=tabx, error.msg=error.msg))
   }
}


############################################
#---------------- function for determining the effective number of independent traits
.nyholt <- function(R){
  k <- nrow(R)
  evals <- eigen(R, only.values=TRUE)$values
  k.eff <- 1 + (k-1)*(1 - var(evals)/k)
  return(k.eff)
}

#--------------------------------------- Main function for minP (Nyholt-Šidák corrected) --------------------------------
minP.sidak <- function(P, R=NULL, method="Nyholt"){
  if(method=="Nyholt" | method=="independent"){
    k <- length(P)
    if(k==1) stop("Only 1 single-trait p-value provided in input parameter 'P'.")
    if(method=="Nyholt"){
      if(is.null(R)) stop("Input trait-correlation matrix 'R' when choosing Nyholt method.")
      if(nrow(R)!=ncol(R) | sum(diag(R))!=nrow(R)) stop("Input 'R' must be a correlation matrix.")
      if(nrow(R)!=k) stop("Mismatch in dimensions of 'P' and 'R'. 'P' is the vector of >1 single-trait genetic association p-values for a given genetic variant while 'R' is the trait correlation matrix.")
      detR <- det(R)
      if(detR<0 | detR>1) stop("Input 'R' must be a correlation matrix.")
      k.nyholt <- .nyholt(R)
      minp <- min(1-(1-P)^k.nyholt)
      return( list(minP=minp, K=k, K.indep=k.nyholt, method=method) )
    }else{
      minp <- min(1-(1-P)^k)
      return( list(minP=minp, K=k, method=method) )
    }
  }else{
    stop("Choice for method includes 'Nyholt' (for correlated traits) and 'independent' (for uncorrelated traits).")
  }
}

# function for obtaining Nyholt-Sidak corrected significance threshold
level.sidak <- function(level=5e-8, R=NULL, method="Nyholt"){
  if(method=="Nyholt" | method=="independent"){
    if(is.null(R)) stop("Input trait-correlation matrix 'R' when choosing Nyholt method. If choosing 'independent', input an identity matrix with K rows, K columns - K being the total number of traits.")
    k <- ncol(R)
    if(method=="Nyholt"){
      if(nrow(R)!=k | sum(diag(R))!=k) stop("Input 'R' must be a correlation matrix.")
      detR <- det(R)
      if(detR<0 | detR>1) stop("Input 'R' must be a correlation matrix.")
      k.nyholt <- .nyholt(R)
      levelc <- 1-(1-level)^(1/k.nyholt)
      return( list(level.corrected=levelc, K=k, K.indep=k.nyholt, method=method) )
    }else{
      levelc <- 1-(1-level)^(1/k)
      return( list(level.corrected=levelc, K=k, method=method) )
    }
  }else{
    stop("Choice for method includes 'Nyholt' (for correlated traits) and 'independent' (for uncorrelated traits).")
  }
}

