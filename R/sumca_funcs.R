##------------------------------------------------------------
###EBLUP and MSPE estimation using Sumca method 
###Authors: Yuanyuan Li, Jiming Jiang
##------------------------------------------------------------
#Compute EBLUP and MSPE for general linear mixed model using Sumca method. Compatible to
#models fitted using `lme4` or `nlme` packages.
#Let y=X\beta+Z\alpha+\varepsilon, where \alpha \sim N(0, G), \varepsilon \sim N(0,R).
#Let \theta=K\beta + L\alpha be the mixed effect of interest. If K and L are
#not provided, the default values are X and Z, \theta becomes response means.

library(mvnfast)
library(broom.mixed)
# Get an EBLUP component A; EBLUP=K\beta+LA
A_func=function(G, R, Z, y, X, beta){
  V=Z %*% G %*% t(Z)+ R#variaince of y
  Vinv=solve(V)
  A=G%*% t(Z) %*% Vinv%*% (y-as.numeric(X %*% beta))
  return(list(A=A, Vinv=Vinv))
}

#Get the central matrix in Var(\theta|y)
B_func=function(G, Z, Vinv){
  B = G- G%*% t(Z) %*% Vinv%*% Z %*% G
  return(B)
}

#Generate new response from fitted model.
resampling=function(beta, G, R, X, Z){
  alpha = rmvn(1, mu=rep(0, nrow(G)), G)[1,]
  error = rmvn(1, mu=rep(0, nrow(R)), R)[1,]
  y = as.numeric(X %*% beta + Z %*% alpha) + error
  return(y)
}

#Main function
#Model is a lme4 or nlme type model object
MSPEsumca=function(model, K=NULL,L=NULL,resamp_k=NULL,data=NULL){
  designM=extract_design(model,data=data)
  attach(designM)
  est=extract_est(model)
  attach(est)
  y=DF$y
  n = length(y)
  if(is.null(K)){K=X}
  if(is.null(L)){L=Z}
  if(is.null(resamp_k)){resamp_k=ngrp}
  A_org = A_func(G, R, Z, y, X, beta)
  B_org = B_func(G, Z, A_org$Vinv)
  V_org = diag(L %*% B_org %*% t(L))#=Var(\theta|y)
  EBLUP= as.numeric(K %*% beta + L %*% A_org$A)
  dphi=matrix(0, nrow = nrow(K), ncol = resamp_k)
  beta_k=matrix(0,nrow=length(beta),ncol=resamp_k)
  sd_k=matrix(0,nrow=length(est$sigma),ncol=resamp_k)
  df_k=DF
  for (k in 1: resamp_k){
    yk=resampling(beta, G, R, X, Z)
    df_k$y=yk
    mk=update(model,data=df_k)
    est_k=extract_est(mk)
    beta_k[,k]=est_k$beta
    sd_k[,k]=est_k$sigma
    G_k= est_k$G
    Rk=est_k$R#error variance
    A_k = A_func(G_k, Rk, Z, yk, X, beta_k[,k])
    B_k = B_func(G_k, Z, A_k$Vinv)
    ak =((K %*% beta_k[,k] + L %*% A_k$A)-
           (K %*% beta+ L %*% A_func(G, R, Z, yk, X, beta)$A))^2+ V_org
    V_k= diag(L %*% B_k %*% t(L))
    dphi[,k]= as.numeric(ak-V_k)
  }
  #est_boot=flatten(est_k, resamp_k)
  sd_beta= apply(beta_k,1,sd)
  names(sd_beta)=names(beta)
  sd_sigma=apply(sd_k, 1,sd)
  names(sd_sigma)=names(sigma)
  MSPE = V_org + apply(dphi, 1, mean)
  return(list(EBLUP=EBLUP, MSPE=MSPE,sd_beta=sd_beta, sd_sigma=sd_sigma))
}


extract_design=function(model,data){
  type=class(model)[1]
  result=list()
  if (type %in%c("lmerMod","glmerMod")){
  if (is.null(data)){
    stop("Please provide the data used to fit this 
         lme4-type model by adding 'data='!")
  }
  result$DF=data
  result$X=as.matrix(getME(model,"X"))
  #Zm=getME(model,"mmList")
  result$Z=as.matrix(getME(model,"Z"))
  result$ngrp=summary(model)$ngrps[1]
  }else{
    result$DF=model$data
    result$X=model.matrix(model,data = model$data)
    Zm= model.matrix(formula(model$modelStruct$reStr)[[1]],data=model$data)
    f <- model$groups[,1]
    Ji<-t(as(f,Class="sparseMatrix"))
    result$Z=t(KhatriRao(t(Ji),t(Zm)))
    result$ngrp=model$dims$ngrps[1]
  }
  return(result)
}

getErrorV=function(phi,form,data,var_e){
cs1AR1 <- corAR1(phi, form=form)
cs1AR1. <- Initialize(cs1AR1, data = data)
vlist=corMatrix(cs1AR1.)
Ematrix=as.matrix(.bdiag(vlist))*var_e
return(Ematrix)
}


extract_est=function(model){
  result=list()
  result$beta=fixef(model)
  type=class(model)[1]
  sum_stats=tidy(model)
  result$sigma=sum_stats$estimate[sum_stats$effect=="ran_pars"]
  var_e=tail(result$sigma,1)^2
  if (type %in%c("lmerMod","glmerMod")){
    sigma_est=VarCorr(model)
    mm=getME(model,"l_i")
    ml=lapply(1:length(mm),function(x)kronecker(diag(mm[x]), sigma_est[[x]]))
    result$G= .bdiag(ml)
    n=getME(model,"devcomp")$dims["n"]
    result$R=var_e*diag(n)
  }else if(type=="lme"){
    sigma_est=getVarCov(model)
    SD=VarCorr(model)
    result$G=kronecker(diag(model$dims$ngrps[1]), sigma_est)
    if(is.null(model$modelStruct$corStruct)){
      result$R=var_e*diag(model$dims$N)
    }else{
      phi=coefficients(model$modelStruct$corStruct,unconstrained = F)
      result$R=getErrorV(phi,form=formula(model$modelStruct$corStruct),
                         model$data,tail(var_e,1))
    }
  }else{
    stop("Please provide a lme4 or nlme type model object!")
  }
  return(result)
}


#extract_est=function(model){
#  result=list()
#  result$beta=fixef(model)
#  type=class(model)[1]
#  if (type %in%c("lmerMod","glmerMod")){
#  sigma_est=VarCorr(model)
#  df.sigma= as.data.frame(sigma_est)
#  result$sigma = df.sigma$sdcor
#  names(result$sigma)=paste(df.sigma$grp,df.sigma$var1,df.sigma$var2)
#  mm=getME(model,"l_i")
#  ml=lapply(1:length(mm),function(x)kronecker(diag(mm[x]), sigma_est[[x]]))
#  result$G= .bdiag(ml)
#  var_e=tail(df.sigma$vcov,1)
#  n=getME(model,"devcomp")$dims["n"]
#  result$R=var_e*diag(n)
#  }else if(type=="lme"){
#    sigma_est=getVarCov(model)
#    SD=VarCorr(model)
#    if (ncol(SD)>2){
#    result$sigma = as.numeric(append(SD[,"StdDev"],
#     na.omit(suppressWarnings(as.numeric(SD[,"Corr"]))),
#      after=nrow(SD)-1))
#    }else{result$sigma = as.numeric(SD[,"StdDev"])}
#    names(result$sigma)=row.names(SD)
#    result$G=kronecker(diag(model$dims$ngrps[1]), sigma_est)
#    if(is.null(model$modelStruct$corStruct)){
#      result$R=tail(result$sigma,1)^2*diag(model$dims$N)
#    }else{
#    phi=coefficients(model$modelStruct$corStruct,unconstrained = F)
#    result$R=getErrorV(phi,form=formula(model$modelStruct$corStruct),
#                      model$data,tail(result$sigma,1))
#    }
#  }else{
#  stop("Please provide a lme4 or nlme type model object!")
#}
#  return(result)
#}


