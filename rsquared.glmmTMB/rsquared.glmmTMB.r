rsquared.glmmTMB <- function(model, method = "trigamma"){
  
  if(model$modelInfo$allForm$ziformula != formula("~0")){
    stop("Zero-inflation not supported yet")
  }
  
  if(model$modelInfo$allForm$dispformula != formula("~1")){
    stop("Explicit dispersion formula not supported yet")
  }
  
  family <- family(model)$family
  
  link <- family(model)$link
  
  X <- model.matrix(model)
  
  sigma <- VarCorr(model)$cond #estimated standard deviation for RE and residuals
  
  sigmaF <- var(as.numeric(fixef(model)$cond %*% t(X))) #variance in fixed effects
  
  data <- model$frame #data used in model
  
  if(family == "gaussian"){
    
    sigmaL <- sum(sapply(1:length(sigma), function(i) {
      
      sigma. <- sigma[[i]]
      
      Z <- as.matrix(X[, rownames(sigma.), drop = FALSE])
      
      sum(rowSums(Z %*% sigma.) * Z) / nrow(X)
      
    } ) )
    
    sigmaE <- attr(sigma, "sc")^2
    
    mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)
    
    con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)
    
    list(family = "gaussian", link = "identity", method = "none", Marginal = mar, Conditional = con)
  }
  
  if(family == "poisson"){
    
    if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")
    
    sigmaL <- sum(GetOLRE.glmmTMB(sigma, model, X, data, RE = "RE"))
    
    sigmaE <- sum(GetOLRE.glmmTMB(sigma, model, X, data, RE = "OLRE"))
    
    rand <- onlyBars(formula(model))
    
    #f <- paste(all.vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
    
    nullmodel <- suppressWarnings(update(model, formula(paste(". ~ 1 +", onlyBars(formula(model))))))
    
    # lambda <- attr(VarCorr(model), "sc")^2
    
    lambda <- exp(fixef(nullmodel)$cond[1] + (sigmaL + sigmaE)/2)
    
    omega <- 1
    
    if(link == "mu^0.5") sigmaD <- 0.25 * omega else {
      
      if(link == "log") {
        
        nu <- omega / lambda
        
        if(method == "delta") sigmaD <- nu
        
        if(method == "lognormal") sigmaD <- log(1 + nu)
        
        if(method == "trigamma") sigmaD <- trigamma(1/nu)
        
      } else stop("Unsupported link function!")
      
    }
  }
  
  if(family == "binomial"){
    
      if(method == "trigamma") method <- "delta"
      
      if(!method %in% c("theoretical", "delta")) stop("Unsupported method!")
      
      sigmaL <- sum(GetOLRE.glmmTMB(sigma, model, X, data, RE = "all"))
      
      sigmaE <- 0
      
      if(method == "theoretical") sigmaD <- pi^2/3
      
      if(method == "delta") {
        
        #rand <- onlyBars(formula(model))
        
        #f <- paste(all.vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
        
        nullmodel <- suppressWarnings(update(model, formula(paste(". ~ 1 +", onlyBars(formula(model))))))
        
        vt <- sum(unlist(VarCorr(nullmodel)))
        
        pmean <- plogis(as.numeric(fixef(nullmodel)$cond) - 0.5 * vt *
                          tanh(as.numeric(fixef(nullmodel)$cond) * (1 + 2 * exp(-0.5 * vt))/6))
        
        sigmaD <- 1/(pmean * (1 - pmean))
        
      }
      
  }
  
  #note for betabinomial only delta method with logit link currently working
  if(family == "betabinomial"){
    
    sigmaL <- sum(GetOLRE.glmmTMB(sigma, model, X, data, RE = "all"))
    
    #rand <- onlyBars(formula(model))
      
    #f <- paste(all.vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
      
    nullmodel <- suppressWarnings(update(model, formula(paste(". ~ 1 +", onlyBars(formula(model))))))
      
      
    lambda <- as.numeric(exp(fixef(nullmodel)$cond + 0.5 * sum(unlist(VarCorr(nullmodel)))))
    
    if(method == "theoretical") sigmaD <- pi^2/3
    
    if(method == "delta") {
      
     
      
    }
      
    theta <- summary(model)$sigma #not so nice, to be replaced by proper parameter extraction
      
      if(link == "logit") {
        
        if(!method %in% c("delta")) stop("Unsupported method!")
        
        nu <- (1/lambda) + (1/theta)
        
        sigmaE <- nu
        
        vt <- sum(unlist(VarCorr(nullmodel)))
        
        pmean <- plogis(as.numeric(fixef(nullmodel)$cond) - 0.5 * vt *
                          tanh(as.numeric(fixef(nullmodel)$cond) * (1 + 2 * exp(-0.5 * vt))/6))
        
        sigmaD <- 1/(pmean * (1 - pmean))
        
      }
      
  }
  
  #similar to negative binomial
  if(family == "compois"){
    if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")
    
    sigmaL <- sum(GetOLRE.glmmTMB(sigma, model, X, data, RE = "RE"))
    
    sigmaE <- sum(GetOLRE.glmmTMB(sigma, model, X, data, RE = "OLRE"))
    
    rand <- onlyBars(formula(model))
    
    f <- paste(all.vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
    
    nullmodel <- suppressWarnings(glmmTMB(formula(f), family = poisson(link = link), data = data))
    
    # lambda <- attr(VarCorr(model), "sc")^2
    
    lambda <- exp(fixef(nullmodel)$cond[1] + (sigmaL + sigmaE)/2)
    
    theta <- summary(model)$sigma
    
    omega <- 1
    
    if(link == "mu^0.5") sigmaD <- 0.25 * omega else {
      
      if(link == "log") {
        
        nu <- (omega / lambda) + (1 / theta)
        
        if(method == "delta") sigmaD <- nu
        
        if(method == "lognormal") sigmaD <- log(1 + nu)
        
        if(method == "trigamma") sigmaD <- trigamma(1/nu)
        
      } else stop("Unsupported link function!")
      
    }
  }
  
  else{
    print("Family not supported yet")
  }
  
  
  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaD + sigmaE)
  
  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaD + sigmaE)
  
  list(family = family, link = link, method = method, Marginal = mar, Conditional = con)
  
}


#' Obtain (observation-level) random effects from a generalized linear mixed model
#' 
#' RE = "all" all random effects are reported
#' RE = "RE" just group effects are reported
#' RE = "OLRE" just observation-level effects are reported
#' 
#' @keywords internal
#' 
GetOLRE.glmmTMB <- function(sigma, model, X, data, RE = c("all", "RE", "OLRE")) {
  
  if(class(model) %in% c("glmmTMB")) {

    rand <- sapply(lme4::findbars(formula(model)), function(x) as.character(x)[3])
    
    rand <- rand[!duplicated(rand)] 
    
  } 
  
  # else if(class(model) %in% c("lme", "glmmPQL")) { }
  
  idx <- sapply(sapply(strsplit(rand, "\\:"), function(x) gsub("\\(|\\)", "", x)), function(x) {
    
    length(unique(data[, x])) == nrow(data)
    
  } )
  
  sigma.names <- unlist(strsplit(names(sigma), "\\."))
  
  idx. <- sapply(sigma.names, function(x) !any(x %in% rand[idx]))
  
  if(RE == "RE") 
    
    sapply(sigma[idx.], function(i) {
      
      Z <- as.matrix(X[, rownames(i), drop = FALSE])
      
      sum(rowSums(Z %*% i) * Z) / nrow(X)
      
    } ) else if(RE == "OLRE") {
      
      if(all(idx == FALSE)) 0 else {
        
        sapply(sigma[idx], function(i) {
          
          Z <- as.matrix(X[, rownames(i), drop = FALSE])
          
          sum(rowSums(Z %*% i) * Z) / nrow(X)
          
        } ) } } else if(RE == "all")
          
          sapply(sigma, function(i) {
            
            Z <- as.matrix(X[, rownames(i), drop = FALSE])
            
            sum(rowSums(Z %*% i) * Z) / nrow(X)
            
          } )
}

#' Get random effects from merMod
#' 
#' @keywords internal
#' 
onlyBars <- function(formula., slopes = TRUE) {
  
  f <- lme4::findbars(formula.)
  
  if(slopes == TRUE) paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ") else {
    
    f <- f[sapply(f, function(x) grepl("1\\||1 \\|", deparse(x)))]
    
    paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ")
    
  }
  
}

#' Get vector of transformed variables
#' 
#' @keywords internal
#' 
all.vars_trans <- function(formula.) {
  
  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(class(formula.) == "formula") {
    
    if(formula.[[3]] == 1) deparse(formula.[[2]]) else {
      
      if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)
      
      c(rownames(attr(terms(formula.), "factors"))[1], labels(terms(formula.)))
      
    }
    
  } else unlist(strsplit(formula., " ~~ "))
  
}
