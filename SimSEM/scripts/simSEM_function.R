####################

# helper functions to run the simSEM script

# author: Lionel Hertzog

# date: 17.10.2018

######


# create matrix with X variables with some relations indicated as ones
relation <- function(X=5,C=0.3){
  out <- matrix(0,nrow=X,ncol=X) #fill the matrix with 0s
  lwr <- X*(X-1)/2 #number of cells in the lower triangle of the matris
  nb.edge <- C*X**2 #number of 1s to add to the matrix
  if(nb.edge>lwr){
    stop("C is set too high to create a DAG, reduce it!")
  }
  while(!all(rowSums(out)[-1]!=0)){#while loop to ensure that all variables have at least one relation
    out <- matrix(0,nrow=X,ncol=X)
    out[lower.tri(out)][sample(x=1:lwr,size=nb.edge,replace=FALSE)] <- 1
    diag(out) <- 0 #ensure that no 1s are in the diagonal
  }
  rownames(out)<-paste0("X",1:X)
  return(out)
}



#turn a matrix with 0s and 1s into a series of R formulas
#rows are ys and columns are xs
formula_list <- function(mat){
  vars <- rownames(mat)
  formL <- list()
  for(i in 2:dim(mat)[1]){
    lhs <- vars[i] #the left-hand side of the formula
    rhs <- vars[mat[i,]==1] #the right-hand side of the formula
    if(length(rhs) == 0) { next }
    else{
      formL[[(i-1)]] <- formula(paste0(lhs,"~",paste(rhs,collapse="+")))
    }
  }
  return(formL[!sapply(formL,is.null)])
}

to_dagitty <- function(mat){
  dims <- dim(mat)
  vars <- rownames(mat)
  ids <- arrayInd(which(mat==1),.dim = dims)
  infs <- NULL
  for(i in 1:nrow(ids)){
    infs <- paste0(infs,paste0(vars[ids[i,2]]," -> ",vars[ids[i,1]]," "))
  }
  dd <- dagitty(paste0("dag{",infs,"}"))
  return(dd)
}

#turn a list of formulas into a list of models
model_list <- function(formL,data,fn=NULL){
  if(is.null(fn)){
    modL <- lapply(formL,function(x) lm(x,data))
  }
  else{
    stop("Not implemented yet") #one could allow the user to give a vector of model types
  }
  return(modL)
}

#create a function to get important infos from a sem object
grab_val <- function(modL,data,pkg="pcsem"){
  if(pkg == "lavaan"){
    semObj <- tryCatch(sem.lavaan(modL,data),error = function(e) NA)
    if(! is.na(semObj)){
      sigM <-ifelse(lavInspect(semObj,"converged"),ifelse(lavInspect(semObj,"fit")["pvalue"]>=0.05,1,0),0) #get the SEM p-value, giving a value of 0 if the model did not converged
      semCoeff <- parameterestimates(semObj)[grep("^~$",parameterestimates(semObj)$op),] #get the model coefficients
      nbPaths <- nrow(semCoeff)
      nbSigPaths <- length(which(semCoeff$pvalue<=0.05)) #how many significant paths
      avgRsq <- mean(inspect(semObj,"r2"),na.rm=TRUE)
      aic <- ifelse(lavInspect(semObj,"converged"),AICc(semObj),NA)
      bic <- ifelse(lavInspect(semObj,"converged"),BIC(semObj),NA)
      hbic <- ifelse(lavInspect(semObj,"converged"),fitMeasures(semObj)["logl"] + log(nrow(data) / (2*pi)) * fitMeasures(semObj)["npar"],NA)
      #avgSigPaths <- mean(abs(semCoeff$est[semCoeff$pvalue<=0.05])) #the average values of significant paths
    }
    else{
      sigM <- NA
      nbPaths <- NA
      nbSigPaths <- NA
      avgRsq <- NA
      aic <- bic <- hbic <- NA
    }
  }
  else{
    semObj <- tryCatch(sem.fit(modL,data,.progressBar = FALSE,conditional = TRUE),error=function(e) NA) #evaluate the independence claims
    if(!is.na(semObj)){
      sigM <-ifelse(semObj$Fisher.C$p.value>=0.05,1,0) #get the SEM p-value
      semCoeff <- sem.coefs(modL,data) #get the model coefficients
      nbPaths <- nrow(semCoeff)
      nbSigPaths <- length(which(semCoeff$p.value<=0.05)) #how many significant paths
      avgRsq <- mean(laply(modL,function(x) summary(x)$adj.r.squared))
      #avgSigPaths <- mean(abs(semCoeff$estimate[semCoeff$p.value<=0.05])) #the average values of significant paths
      aic <- sem.aic(modL,data,.progressBar = FALSE)$AICc
      bic <- semObj$Fisher.C[1] + semObj$AIC["K"] * log(semObj$AIC["n"])
      hbic <- semObj$Fisher.C[1] + log(nrow(data) / (2*pi)) * semObj$AIC["K"]
    }
    else{
      sigM <- NA
      nbPaths <- NA
      nbSigPaths <- NA
      avgRsq <- NA
      aic <- bic <- hbic <- NA
    }
  }
  out <- data.frame(N=nrow(data),X=ncol(data),sigM=sigM,nbPaths=nbPaths,nbSigPaths=nbSigPaths,avgRsq=avgRsq,AIC=aic,BIC=as.numeric(bic),hbic=as.numeric(hbic))
  return(out)
}

#a function to generate data based on diverse fidelity from the simulated list of models
generate_data <- function(relL,N=20,X=5,type="random",p_var=0.25,sd_eff=2.5,sd_res=1, out_dat = TRUE){
  formL <- formula_list(relL)
  if(type=="random"){
    dat <- as.data.frame(matrix(rnorm(X*N,0,sd_res),ncol=X,dimnames = list(1:N,paste0("X",1:X))))
    relL_m <- relL
  }
  else if(type=="exact"){
    dat <- data.frame(X1=runif(N,-2,2))
    for(f in formL){
      allv <- all.vars(f)
      y <- allv[1]
      xs <- dat[,allv[2:length(allv)]]
      xs <- as.matrix(cbind(rep(1,N),xs))
      coeff <- rnorm(n = ncol(xs),0,sd_eff)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,sd_res))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
      relL_m <- relL
    }
  }
  #reverse a certain proportion of the causal relation, such as X1 -> X2 becomes X2 -> X1
  else if(type=="shuffled"){
    relL_m <- relL
    #how many relation to shuffle?
    #default to 0.25
    n_s <- ifelse(ceiling(sum(relL)*p_var)<1,1,ceiling(sum(relL)*p_var))
    #grab cell index to shuffle
    to_s <- sample(which(relL == 1),n_s,replace=FALSE)
    #grab row and column id for these cells
    id_s <- arrayInd(to_s,dim(relL))
    #transpose the 1s up the diagonal
    relL_m[id_s] <- 0
    relL_m[id_s[,c(2,1)]] <- 1
    #add column names
    dimnames(relL_m)[[2]] <- dimnames(relL_m)[[1]]
    if(!is_dag(graph_from_adjacency_matrix(relL_m))){
      while(!is_dag(graph_from_adjacency_matrix(relL_m))){
        #grab cell index to shuffle
        to_s <- sample(which(relL_m == 1),n_s,replace=FALSE)
        #grab row and column id for these cells
        id_s <- arrayInd(to_s,dim(relL))
        #transpose the 1s up the diagonal
        relL_m[id_s] <- 0
        relL_m[id_s[,c(2,1)]] <- 1
        #add column names
        dimnames(relL_m)[[2]] <- dimnames(relL_m)[[1]]
      }
    }
    #sort back the DAG
    od <- topo_sort(graph_from_adjacency_matrix(relL_m),mode="in")
    relL_m <- relL_m[od,od]
    formL_m <- formula_list(relL_m)
    #a random vector for the exogenous variable(s) to start of the data creation
    exo <- which(rowSums(relL_m) == 0)
    dat <- as.data.frame(matrix(runif(N * length(exo),-2,2),ncol=length(exo)))
    names(dat) <- names(exo)
    for(f in 1:length(formL_m)){
      allv <- all.vars(formL_m[[f]])
      y <- allv[1]
      xs <- dat[,allv[2:length(allv)]]
      xs <- as.matrix(cbind(rep(1,N),xs))
      coeff <- rnorm(n = ncol(xs),0,sd_eff)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,sd_res))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
    }
  }
  else if(type == "overspecified"){ #drop some relation in the data creation, ie model will fit unneedded relations
    relL_m <- relL
    #how many relation to drop?
    #default to 0.25
    n_d <- ifelse(ceiling(sum(relL)*p_var)<1,1,ceiling(sum(relL)*p_var))
    #drop relation ensuring that all variables X2-Xn still have some relations
    #grab cell index to drop
    to_d <- sample(which(relL == 1),n_d,replace=FALSE)
    #grab row and column id for these cells
    id_d <- arrayInd(to_d,dim(relL))
    relL_m[id_d] <- 0
    while(any(rowSums(relL_m[-1,]) == 0)){
      relL_m <- relL
      #grab cell index to drop
      to_d <- sample(which(relL == 1),n_d,replace=FALSE)
      #grab row and column id for these cells
      id_d <- arrayInd(to_d,dim(relL))
      relL_m[id_d] <- 0
    }
    #generate the data
    formL_m <- formula_list(relL_m)
    dat <- data.frame(X1=runif(N,-2,2))
    for(f in formL_m){
      allv <- all.vars(f)
      y <- allv[1]
      xs <- dat[,allv[2:length(allv)]]
      xs <- as.matrix(cbind(rep(1,N),xs))
      coeff <- rnorm(n = ncol(xs),0,sd_eff)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,sd_res))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
    }
  }
  else if(type == "underspecified"){ #add some more relation than previously fitted
    relL_m <- relL
    #how many relation to add?
    #default to 0.25
    n_a <- ifelse(ceiling(sum(relL)*p_var)<1,1,ceiling(sum(relL)*p_var))
    #add relation ensuring that X1 still is an exogenous variable
    #grab cell index to drop
    to_a <- sample(which(relL[lower.tri(relL)] == 0),n_a,replace=FALSE)
    #grab row and column id for these cells
    relL_m[lower.tri(relL)][to_a] <- 1
    #generate the data
    formL_m <- formula_list(relL_m)
    dat <- data.frame(X1=runif(N,-2,2))
    for(f in formL_m){
      allv <- all.vars(f)
      y <- allv[1]
      xs <- dat[,allv[2:length(allv)]]
      xs <- as.matrix(cbind(rep(1,N),xs))
      coeff <- rnorm(n = ncol(xs),0,sd_eff)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,sd_res))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
    }
  }
  else{
    stop("Not implemented yet") #one could create datasets with non-random structure, to be implemented
  }
  if(out_dat){
    return(dat)
  }
  else{
    return(formula_list(relL_m))
  }
}


#put it all in one function to replicate over
sim_sem <- function(N=20,X=5,C=0.3,type="random",p_var=0.25,pkg="lavaan",sd_eff = 2.5,sd_res = 1){
  relL <- relation(X=X,C=C)
  formL <- formula_list(relL)
  dagg <- to_dagitty(relL) #get DAG
  dat <- generate_data(N=N,X=X,relL=relL,type=type,p_var=p_var,sd_eff=sd_eff,sd_res=sd_res)
  modL <- model_list(formL,dat)
  out <- grab_val(modL,dat,pkg=pkg)
  #get proportion of conditional independence violated
  ll <- localTests(dagg,dat)
  out$Local <- sum(ll$p.value < 0.05) / length(ll$p.value)
  out$Sd_eff <- sd_eff
  out$Sd_res <- sd_res
  out$pkg <- pkg
  out$type <- type
  out$C <- C
  return(out)
}

# a new simsem function to truly compare models, step 1: generate exact data, step 2: create shuffled / under or overspecified model
# potentially with varying p_var, step 3: grab_val on actual model_list and on model_list generated from shuffled / under/over specified model

sim_sem_c <- function(N = 20, X = 5, C = 0.3, p_var = 0.25, pkg = "pcsem", sd_eff = 2.5, sd_res = 1){
  relL <- relation(X=X,C=C)
  formL <- formula_list(relL)
  
  # exact data generation
  dat <- generate_data(N=N,X=X,relL=relL,type="exact",p_var=p_var,sd_eff=sd_eff,sd_res=sd_res)
  # generate wrong formLs
  formL_s <- generate_data(N=N,X=X,relL=relL,type="shuffled",p_var=p_var,sd_eff=sd_eff,sd_res=sd_res,out_dat=FALSE)
  formL_u <- generate_data(N=N,X=X,relL=relL,type="underspecified",p_var=p_var,sd_eff=sd_eff,sd_res=sd_res,out_dat=FALSE)
  formL_o <- generate_data(N=N,X=X,relL=relL,type="overspecified",p_var=p_var,sd_eff=sd_eff,sd_res=sd_res,out_dat=FALSE)
  # true modL
  modL_t <- model_list(formL,dat)
  #wrong modL
  modL_s <- model_list(formL_s, dat)
  modL_u <- model_list(formL_u, dat)
  modL_o <- model_list(formL_o, dat)
  
  # grab vals
  out_t <- grab_val(modL_t,dat,pkg=pkg)
  out_s <- grab_val(modL_s,dat,pkg=pkg)
  out_u <- grab_val(modL_u,dat,pkg=pkg)
  out_o <- grab_val(modL_o,dat,pkg=pkg)
  
  # put IC together
  out_f <- data.frame(N=N,X=X,sd_eff=sd_eff,sd_res=sd_res,C=C,p_var=p_var,
                      AIC_t=out_t$AIC,AIC_s=out_s$AIC,AIC_u=out_u$AIC,AIC_o=out_o$AIC,
                      BIC_t=out_t$BIC,BIC_s=out_s$BIC,BIC_u=out_u$BIC,
                      BIC_o=out_o$BIC,hbic_t=out_t$hbic,hbic_s=out_s$hbic,
                      hbic_u=out_u$hbic,hbic_o=out_o$hbic,stringsAsFactors = FALSE)
  return(out_f)
}

#a function to loop through the parameters
sim_for <- function(sims,n_rep=100){
  res <- NULL
  for(i in 1:nrow(sims)){
    tmp <- replicate(n_rep,sim_sem(N = sims[i,"N"],X=sims[i,"X"],
                                   C=sims[i,"C"],type=sims[i,"type"],
                                   pkg=sims[i,"pkg"], sd_eff = sims[i,"sd_eff"], 
                                   sd_res = sims[i,"sd_res"]))
    out <- as.data.frame(matrix(unlist(tmp),nrow = n_rep,byrow=TRUE,dimnames=list(1:n_rep,attr(tmp,"dimnames")[[2]])))
    out$Exp <- i
    res <- rbind(res,out)
    print(paste0(sims[i,"type"]," : ",i, " out of ",nrow(sims)))
  }
  return(res)
}

# a second sim_for to re-do the model comparison
# the sims object as this time no type column
sim_for_c <- function(sims,n_rep=5){
  res <- NULL
  for(i in 1:nrow(sims)){
    tmp <- replicate(n_rep,sim_sem_c(N = sims[i,"N"],X=sims[i,"X"],
                                     C=sims[i,"C"],p_var=sims[i,"p_var"],
                                     sd_eff = sims[i,"sd_eff"], 
                                     sd_res = sims[i,"sd_res"]))
    out <- as.data.frame(matrix(unlist(tmp),nrow = n_rep,byrow=TRUE,dimnames=list(1:n_rep,attr(tmp,"dimnames")[[2]])),stringsAsFactors = FALSE)
    out$Exp <- i
    res <- rbind(res,out)
    print(paste0("Row : ",i, " out of ",nrow(sims)))
  }
  return(res)
}
