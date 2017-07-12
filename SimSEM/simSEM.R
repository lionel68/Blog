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

#look at a couple of graphs that are created by the functions
relL <- list(relation(X=5,C=0.3),relation(X=5,C=0.2),relation(X=10,C=0.3),relation(X=10,C=0.4))
datL <- list(as.data.frame(matrix(rnorm(5*50,0,1),ncol=5,dimnames = list(1:50,paste0("X",1:5)))),
             as.data.frame(matrix(rnorm(5*50,0,1),ncol=5,dimnames = list(1:50,paste0("X",1:5)))),
             as.data.frame(matrix(rnorm(10*50,0,1),ncol=10,dimnames = list(1:50,paste0("X",1:10)))),
             as.data.frame(matrix(rnorm(10*50,0,1),ncol=10,dimnames = list(1:50,paste0("X",1:10)))))
formLL <- lapply(relL,function(x) formula_list(x))
modLL <- mapply(formL=formLL,data=datL,function(formL,data) model_list(formL,data))
#the figure
par(mfrow=c(2,2))
mapply(modL=modLL,data=datL,function(modL,data) sem.plot(modL,data,scaling=NA))

#create a function to get important infos from a sem object
grab_val <- function(modL,data,lavaan=FALSE){
  if(lavaan){
    lv <- "lavaan"
    semObj <- sem.lavaan(modL,data)
    sigM <-ifelse(semObj@test[[1]]$pvalue>=0.05,1,0) #get the SEM p-value
    semCoeff <- parameterestimates(semObj)[grep("^~$",parameterestimates(semObj)$op),] #get the model coefficients
    nbPaths <- nrow(semCoeff)
    nbSigPaths <- length(which(semCoeff$pvalue<=0.05)) #how many significant paths
    avgSigPaths <- mean(abs(semCoeff$est[semCoeff$pvalue<=0.05])) #the average values of significant paths
  }
  else{
    lv <- "pcSEM"
    semObj <- sem.fit(modL,data,.progressBar = FALSE) #evaluate the independence claims
    sigM <-ifelse(semObj$Fisher.C$p.value>=0.05,1,0) #get the SEM p-value
    semCoeff <- sem.coefs(modL,data) #get the model coefficients
    nbPaths <- nrow(semCoeff)
    nbSigPaths <- length(which(semCoeff$p.value<=0.05)) #how many significant paths
    avgSigPaths <- mean(abs(semCoeff$estimate[semCoeff$p.value<=0.05])) #the average values of significant paths
  }
  out <- data.frame(N=nrow(data),X=ncol(data),C=0.3,type=lv,sigM=sigM,nbPaths=nbPaths,nbSigPaths=nbSigPaths,avgSigPaths=avgSigPaths)
  return(out)
}
 
#put it all in one function to replicate over
sim_sem <- function(N=20,X=5,C=0.3,type="random",lv=FALSE){
  relL <- relation(X=X,C=C)
  formL <- formula_list(relL)
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
  #all causal directions (expect those on the first variable) are reversed
  #TODO: only reverse a smaller portion of the causal directions
  else if(type=="shuffled"){
    relL_m <- relL
    to_s <- ifelse(floor(sum(relL)*0.2)<1,1,floor(sum(relL)*0.2))
    #grab values to shuffle
    
    relL_m[-1,-1] <- t(relL_m[-1,-1])
    formL_m <- formula_list(relL_m)
    dat <- data.frame(X1=runif(N,-2,2))
    for(f in length(formL_m):1){
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
  else{
    stop("Not implemented yet") #one could create datasets with non-random structure, to be implemented
  }
  modL <- model_list(formL,dat)
  out <- grab_val(modL,dat,lavaan=lv)
  return(out)
}

#test over different sample size (N) and number of variables (X)
sims <- expand.grid(N=seq(20,100,10),X=5:10,C=0.3,lv=FALSE)
res <- NULL
for(i in 1:54){
  x <- as.numeric(sims[i,])
  tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3],lv=x[4]))
  out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
  out$Exp <- i
  res <- rbind(res,out)
  print(i)
}
#now same stuff for lavaan
sims$lv<-TRUE
reslv <- NULL
for(i in 1:54){
  x <- as.numeric(sims[i,])
  tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3],lv=x[4]))
  out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
  out$Exp <- i
  reslv <- rbind(reslv,out)
  print(i)
}
res$type<-"pcsem"
reslv$type<-"lv"
resa<-rbind(res,reslv)
resa%>%
  group_by(Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),AvgSigPaths=mean(avgSigPaths,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths)->res_dd

#some plots
ggplot(res_dd,aes(x=N,y=PropSigM,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Proportion of model accepted (p>0.05)")
ggplot(res_dd,aes(x=N,y=PropSigPaths,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Proportion of significant paths (p<0.05)")
ggplot(res_dd,aes(x=N,y=AvgSigPaths,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Average value of significant paths")


#now simulate the exact distribution for the fitted model
sims <- expand.grid(N=seq(20,100,10),X=5:10,C=0.3,lv=FALSE)
res <- NULL
for(i in 1:54){
  x <- as.numeric(sims[i,])
  tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3],lv=x[4],type="exact"))
  out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
  out$Exp <- i
  res <- rbind(res,out)
  print(i)
}
#now same stuff for lavaan
sims$lv<-TRUE
reslv <- NULL
for(i in 1:54){
  x <- as.numeric(sims[i,])
  tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3],lv=x[4],type="exact"))
  out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
  out$Exp <- i
  reslv <- rbind(reslv,out)
  print(i)
}
res$type<-"pcsem"
reslv$type<-"lv"
resa<-rbind(res,reslv)
resa%>%
  group_by(Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),AvgSigPaths=mean(avgSigPaths,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths)->res_dd

#some plots
ggplot(res_dd,aes(x=N,y=PropSigM,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Proportion of model accepted (p>0.05)")
ggplot(res_dd,aes(x=N,y=PropSigPaths,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Proportion of significant paths (p<0.05)")
ggplot(res_dd,aes(x=N,y=AvgSigPaths,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Average value of significant paths")

#now simulate with shuffling the causal directions
sims <- expand.grid(N=seq(20,100,10),X=5:10,C=0.3,lv=FALSE)
res <- NULL
for(i in 1:54){
  x <- as.numeric(sims[i,])
  tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3],lv=x[4],type="shuffled"))
  out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
  out$Exp <- i
  res <- rbind(res,out)
  print(i)
}
#now same stuff for lavaan
sims$lv<-TRUE
reslv <- NULL
for(i in 1:54){
  x <- as.numeric(sims[i,])
  tmp <- replicate(100,sim_sem(N = x[1],X=x[2],C=x[3],lv=x[4],type="exact"))
  out <- as.data.frame(matrix(unlist(tmp),nrow = 100,byrow=TRUE,dimnames=list(1:100,attr(tmp,"dimnames")[[2]])))
  out$Exp <- i
  reslv <- rbind(reslv,out)
  print(i)
}
res$type<-"pcsem"
reslv$type<-"lv"
resa<-rbind(res,reslv)
resa%>%
  group_by(Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),AvgSigPaths=mean(avgSigPaths,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths)->res_dd

#some plots
ggplot(res_dd,aes(x=N,y=PropSigM,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Proportion of model accepted (p>0.05)")
ggplot(res_dd,aes(x=N,y=PropSigPaths,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Proportion of significant paths (p<0.05)")
ggplot(res_dd,aes(x=N,y=AvgSigPaths,color=type))+geom_path()+facet_wrap(~X)+
  labs(x="Sample size",y="Average value of significant paths")

