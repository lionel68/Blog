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
grab_val <- function(modL,data,C,lavaan=FALSE){
  if(lavaan){
    lv <- "lavaan"
    semObj <- sem.lavaan(modL,data)
    sigM <-ifelse(!semObj@optim$converged,0,ifelse(semObj@test[[1]]$pvalue>=0.05,1,0)) #get the SEM p-value, giving a value of 0 if the model did not converged
    semCoeff <- parameterestimates(semObj)[grep("^~$",parameterestimates(semObj)$op),] #get the model coefficients
    nbPaths <- nrow(semCoeff)
    nbSigPaths <- length(which(semCoeff$pvalue<=0.05)) #how many significant paths
    avgRsq <- mean(laply(modL,function(x) summary(x)$adj.r.squared))
    #avgSigPaths <- mean(abs(semCoeff$est[semCoeff$pvalue<=0.05])) #the average values of significant paths
  }
  else{
    lv <- "pcSEM"
    semObj <- sem.fit(modL,data,.progressBar = FALSE,conditional = TRUE) #evaluate the independence claims
    sigM <-ifelse(semObj$Fisher.C$p.value>=0.05,1,0) #get the SEM p-value
    semCoeff <- sem.coefs(modL,data) #get the model coefficients
    nbPaths <- nrow(semCoeff)
    nbSigPaths <- length(which(semCoeff$p.value<=0.05)) #how many significant paths
    avgRsq <- mean(laply(modL,function(x) summary(x)$adj.r.squared))
    #avgSigPaths <- mean(abs(semCoeff$estimate[semCoeff$p.value<=0.05])) #the average values of significant paths
  }
  out <- data.frame(N=nrow(data),X=ncol(data),C=C,type=lv,sigM=sigM,nbPaths=nbPaths,nbSigPaths=nbSigPaths,avgRsq=avgRsq)
  return(out)
}

#a function to generate data based on diverse fidelity from the simulated list of models
generate_data <- function(N,X,relL,formL,type,p_var){
  if(type=="random"){
    dat <- as.data.frame(matrix(rnorm(X*N,0,1),ncol=X,dimnames = list(1:N,paste0("X",1:X))))
  }
  else if(type=="exact"){
    dat <- data.frame(X1=runif(N,-2,2))
    for(f in formL){
      allv <- all.vars(f)
      y <- allv[1]
      xs <- dat[,allv[2:length(allv)]]
      xs <- as.matrix(cbind(rep(1,N),xs))
      coeff <- rcauchy(n = ncol(xs),0,2.5)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,1))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
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
    relL_mo <- relL_m[od,od]
    formL_m <- formula_list(relL_mo)
    #a random vector for the exogenous variable(s) to start of the data creation
    exo <- which(rowSums(relL_mo) == 0)
    dat <- as.data.frame(matrix(runif(N * length(exo),-2,2),ncol=length(exo)))
    names(dat) <- names(exo)
    for(f in 1:length(formL_m)){
      allv <- all.vars(formL_m[[f]])
      y <- allv[1]
      xs <- dat[,allv[2:length(allv)]]
      xs <- as.matrix(cbind(rep(1,N),xs))
      coeff <- rcauchy(n = ncol(xs),0,2.5)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,1))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
    }
  }
  else if(type == "complex"){ #drop some relation in the data creation, ie model will fit unneedded relations
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
      coeff <- rcauchy(n = ncol(xs),0,2.5)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,1))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
    }
  }
  else if(type == "simple"){ #add some more relation than previously fitted
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
      coeff <- rcauchy(n = ncol(xs),0,2.5)
      dat<-cbind(dat,rnorm(N,xs%*%coeff,1))
      colnames(dat)[ncol(dat)]<-y
      dat <- as.data.frame(apply(dat,2,scale))
    }
  }
  else{
    stop("Not implemented yet") #one could create datasets with non-random structure, to be implemented
  }
  return(dat)
}


#put it all in one function to replicate over
sim_sem <- function(N=20,X=5,C=0.3,type="random",p_var=0.25,lv=FALSE){
  relL <- relation(X=X,C=C)
  formL <- formula_list(relL)
  dagg <- to_dagitty(relL) #get DAG
  dat <- generate_data(N,X,relL,formL,type,p_var)
  modL <- model_list(formL,dat)
  out <- grab_val(modL,dat,C,lavaan=lv)
  #get proportion of conditional independence violated
  ll <- localTests(dagg,dat)
  out$Local <- sum(ll$p.value < 0.05) / length(ll$p.value)
  return(out)
}

#a function to loop through the parameters
sim_for <- function(sims,type="random"){
  res <- NULL
  for(i in 1:nrow(sims)){
    x <- as.numeric(sims[i,])
    tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3],type=type,lv=x[4]))
    out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
    out$Exp <- i
    out$type <- type
    out$pkg <- ifelse(x[4],"lavaan","pcsem")
    res <- rbind(res,out)
    print(paste0(type," : ",i, " out of ",nrow(sims)))
  }
  return(res)
}
