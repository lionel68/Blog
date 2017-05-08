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
    formL[[(i-1)]] <- formula(paste0(lhs,"~",paste(rhs,collapse="+")))
  }
  return(formL)
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
grab_val <- function(modL,data){
  semObj <- sem.fit(modL,data,.progressBar = FALSE) #evaluate the independence claims
  sigM <-ifelse(semObj$Fisher.C$p.value>=0.05,1,0) #get the SEM p-value
  semCoeff <- sem.coefs(modL,data) #get the model coefficients
  nbPaths <- nrow(semCoeff)
  nbSigPaths <- length(which(semCoeff$p.value<=0.05)) #how many significant paths
  avgSigPaths <- mean(semCoeff$estimate[semCoeff$p.value<=0.05]) #the average values of significant paths
  out <- data.frame(N=nrow(data),X=ncol(data),C=0.3,sigM=sigM,nbPaths=nbPaths,nbSigPaths=nbSigPaths,avgSigPaths=avgSigPaths)
  return(out)
}
 
#put it all in one function to replicate over
sim_sem <- function(N=20,X=5,C=0.3,type="random"){
  formL <- formula_list(relation(X=X,C=C))
  if(type=="random"){
    dat <- as.data.frame(matrix(rnorm(X*N,0,1),ncol=X,dimnames = list(1:N,paste0("X",1:X))))
  }
  else{
    stop("Not implemented yet") #one could create datasets with non-random structure, to be implemented
  }
  modL <- model_list(formL,dat)
  out <- grab_val(modL,dat)
  return(out)
}

#test over different sample size (N) and number of variables (X)
sims <- expand.grid(N=seq(20,100,10),X=5:10,C=0.3)
res <- NULL
for(i in 1:54){
  x <- as.numeric(sims[i,])
  tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3]))
  out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
  out$Exp <- i
  res <- rbind(res,out)
  print(i)
}

#average over the different replicates
res%>%
  group_by(Exp,N,X)%>%
  summarise(PropSigM=sum(sigM)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),AvgSigPaths=mean(avgSigPaths,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths)->res_dd

#some plots
ggplot(res_dd,aes(x=N,y=PropSigM))+geom_path()+facet_wrap(~X)
ggplot(res_dd,aes(x=N,y=PropSigPaths))+geom_path()+facet_wrap(~X)
ggplot(res_dd,aes(x=N,y=AvgSigPaths))+geom_path()+facet_wrap(~X)
