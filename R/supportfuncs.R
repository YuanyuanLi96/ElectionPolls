##------------------------------------------------------------
###EBLUP and MSPE estimation using Sumca method 
###Authors: Yuanyuan Li, Jiming Jiang
##------------------------------------------------------------

A_func=function(G, var_e, Z, y, X, beta,n){
  A=G%*% t(Z) %*% solve(Z %*% G %*% t(Z)+ 
                          var_e * diag(n))%*% (y-X %*% beta)
  return(A)
}
B_func=function(G, var_e, Z, n){
  B = G- G%*% t(Z) %*% solve(Z %*% G %*% t(Z)+ 
                               var_e * diag(n))%*% Z %*% G
  return(B)
}
resampling=function(beta, G, var_e, X, Z,n){
  v = mvrnorm(1, mu=rep(0, nrow(G)),Sigma = G)
  y = X %*% beta + Z %*% v + rnorm(n,sd=sqrt(var_e))
  return(y)
}
flatten = function(alist,K){
  newlist=alist[[1]]
  if (K >1){
    for(j in 2:K){
      newlist=Map(rbind,newlist,alist[[j]])
    }
  }
  return(newlist)
}
#theta=K%*%beta + L%*% v
MSPE_k=function(y,X, Z, K,L,resamp_k, Zm, grp,est){
  n = length(y)
  beta=est$beta
  G= est$G
  var_e = est$var_e
  A_org = A_func(G, var_e, Z, y, X, beta,n)
  B_org = B_func(G, var_e, Z, n)
  V_org = diag(L %*% B_org %*% t(L))
  EBLUP= as.numeric(K %*% beta + L %*% A_org)
  dphi=matrix(0, nrow = nrow(K), ncol = resamp_k)
  rl=list()
  for(i in 1:length(Zm)) { 
    nam_z <- paste("Zm", i, sep = "")
    nam_g= paste("grp", i, sep = "")
    assign(nam_z, Zm[[i]])
    assign(nam_g, grp[[i]])
    rl[[i]]=paste("(0+", nam_z, "|", nam_g,")")
  }
  fr= paste(rl,collapse = " + ")
  fk=as.formula(paste("yk~0+X +", fr))
  est_k= list()
  for (k in 1: resamp_k){
    yk=resampling(beta, G, var_e, X, Z,n)
    mk=lmer(fk)
    est_k[[k]]=extract_est(mk)
    beta_k=est_k[[k]]$beta
    G_k= est_k[[k]]$G
    var_e_k=est_k[[k]]$var_e
    A_k = A_func(G_k, var_e_k, Z, yk, X, beta_k,n)
    B_k = B_func(G_k, var_e_k, Z, n)
    ak =((K %*% beta_k+ L %*% A_k)-(K %*% beta+ L %*% 
                                      A_func(G, var_e, Z, yk, X, beta,n)))^2 + V_org
    V_k= diag(L %*% B_k %*% t(L))
    dphi[,k]= as.numeric(ak-V_k)
  }
  est_boot=flatten(est_k, resamp_k)
  #sd_beta= apply(est_boot$beta,2,sd)
  #sd_sigma=apply(est_boot$sigma, 2,sd)
  MSPE = V_org + apply(dphi, 1, mean)
  return(list(EBLUP=EBLUP, MSPE=MSPE))
}

#grp could include more than one term
lmm_random=function(y,model,grp,resamp_k){
  X=as.matrix(getME(model,"X"))
  Z=as.matrix(getME(model,"Z"))
  Zm=getME(model,"mmList")
  est_org=extract_est(model)
  lmm_result=MSPE_k(y,X,Z=Z,K=X,L=Z,resamp_k=resamp_k, 
                    Zm=Zm, grp=grp,est=est_org)
  return(lmm_result)
}

extract_est=function(model){
  result=list()
  result$beta=fixef(model)
  sigma_est=VarCorr(model)
  df.sigma= as.data.frame(sigma_est)
  mm=getME(model,"l_i")
  ml=lapply(1:length(mm),function(x)kronecker(diag(mm[x]), sigma_est[[x]]))
  result$G= .bdiag(ml)
  result$sigma = df.sigma$sdcor
  result$var_e=tail(df.sigma$vcov,1)
  return(result)
}
plotall=function(lmm_result, data,labels="state"){
  #labels includes the text to print in EBLUP figures
  result=list()
  ####Prediction interval
  data$fitted=lmm_result$EBLUP
  data$MSPE=lmm_result$MSPE
  result$fit=data.frame(EBLUP=data$fitted, MSPE=data$MSPE)
  result$mspe_box=boxplot(data$MSPE, main="Distribution of MSPE")
  ct=qnorm(0.95)
  ##EBLUP
  result$eblup_pp=ggpaired(data, x = "candidate", y = "fitted",
                           xlab = FALSE,
                           ylab = "Predicted outcome",
                           color = "candidate", line.color = "gray", line.size = 0.4,
                           facet.by = "year", short.panel.labs = F,label=labels,
                           font.label = list(size = 8),
                           label.select = list(criteria="(`y` >0.2  & `x` == 'R')"))+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))
  ##Plot by states
  mixed=data[c("year","state","candidate","fitted",
               "MSPE")]
  mixed=unique(mixed)
  result$mspe_states=ggplot(mixed, aes(year, fitted)) +
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color = candidate), position = position_dodge(0.3)) +
    facet_wrap(~state)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07")) 
  mixed[,"state"]=as.numeric(mixed[,"state"])
  result$mspe_cycle=ggplot(mixed, aes(state, fitted)) +
    geom_errorbar(
      aes(ymin = fitted-ct*sqrt(MSPE), ymax = fitted+ct*sqrt(MSPE), color = candidate,
          linetype=candidate),
      position = position_dodge(0.3), width = 0.2)+
    geom_point(aes(color =candidate, shape=candidate), position = position_dodge(0.3)) +
    #facet_grid(rows = vars(year))+
    facet_wrap(~year)+
    scale_color_manual(values = c("#0073C2FF", "#FC4E07"))+
    ylab("Prediction Interval")
  return(result)
}

