focusedMR <- function(bx,by,sx,sy,V,k0,alpha,gamma){
  alpha.opts = c(alpha/4,2*alpha/4,3*alpha/4)
  delta.opts = rep(alpha,3)-alpha.opts
  
  # re-order summary data
  Bx <- c(bx[V],bx[-V]); By <- c(by[V],by[-V]); Sx <- c(sx[V],sx[-V]); Sy <- c(sy[V],sy[-V]); V <- 1:length(V)
  
  if(k0>1){
    # candidate sets of additional instruments
    sel <- kmeans(By[-V]/Bx[-V], k0, iter.max = 10, nstart = 1); sel <- sel$cluster
    Rs <- vector(,length=k0) 
    for (k in 1:k0){Rs[k] <- length(combn(1:k0, k, simplify=FALSE))}
    Rs1 <- vector(,length=k0)
    for (k in 1:k0){Rs1[k] <- sum(Rs[1:k])}
    SA <- list()
    for (k in 1:Rs[1]){SA[[k]] <- combn(1:k0, 1, simplify=FALSE)[[k]]}
    for (l in 2:k0){for (k in 1:Rs[l]){SA[[k+Rs1[(l-1)]]] <- combn(1:k0, l, simplify=FALSE)[[k]]}}
    S <- list() 
    for (s in 1:k0){S[[s]] <- length(V)+which(sel==s)}
    for (s in (k0+1):length(SA)){S[[s]] <- unlist(S[SA[[s]]])}
    K=length(S)
    rm(sel, SA, Rs, Rs1)
  } else {
    S <- list(); S[[1]] <- (length(V)+1):length(Bx); K=length(S)
  }
  
  # valid estimator
  Q <- function(tet){sum(((By[V]-(Bx[V]*tet))^2)/(Sy[V]^2 + (tet^2)*Sx[V]^2))}
  init.val <- seq(-1,1,0.2)
  Q.init <- vector(,length=length(init.val))
  for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value}
  tet_v <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
  eta_v <- sum(((Bx[V]^2)-Sx[V]^2)/((Sy[V]^2)+((tet_v^2)*(Sx[V]^2))))
  ci_v <- sum((Sx[V]^2)*(Sy[V]^2)/(((Sy[V]^2)+((tet_v^2)*(Sx[V]^2)))^2))
  var_v <- ((1/eta_v)+(ci_v/eta_v^2))
  eta_v <- sum(((Bx[V]^2)-Sx[V]^2)/((Sy[V]^2)+((((tet_v^2)-var_v))*(Sx[V]^2))))
  ci_v <- sum((Sx[V]^2)*(Sy[V]^2)/(((Sy[V]^2)+((((tet_v^2)-var_v))*(Sx[V]^2)))^2))
  var_v <- ((1/eta_v)+(ci_v/eta_v^2))
  
  # estimators using additional instruments
  tet_s <- vector(,length=K)
  for (k in 1:K){
    Q <- function(tet){sum(((By[c(V,S[[k]])]-(Bx[c(V,S[[k]])]*tet))^2)/(Sy[c(V,S[[k]])]^2 + (tet^2)*Sx[c(V,S[[k]])]^2))}
    init.val <- seq(-1,1,0.2)
    Q.init <- function(l){optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value}
    Q.init <- sapply(1:length(init.val),Q.init)
    tet_s[k] <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
  }
  
  K0 <- 1:K
  
  eta_s <- vector(,length=K); ci_s <- vector(,length=K)
  for(k in 1:K){
    eta_s[k] <- sum(((Bx[S[[k]]]^2)-Sx[S[[k]]]^2)/((Sy[S[[k]]]^2)+((((tet_v^2)-var_v))*(Sx[S[[k]]]^2))))
    ci_s[k] <- sum((Sx[S[[k]]]^2)*(Sy[S[[k]]]^2)/(((Sy[S[[k]]]^2)+((((tet_v^2)-var_v))*(Sx[S[[k]]]^2)))^2))
  }
  
  # bias estimator
  b_s <- function(k){sum((Bx[S[[k]]]*By[S[[k]]]-tet_v*((Bx[S[[k]]]^2)-(Sx[S[[k]]]^2)))/((Sy[S[[k]]]^2)+((((tet_v^2)-var_v))*(Sx[S[[k]]]^2))))}
  xi_s <- function(k){2*(((tet_v^2)-var_v))*sum((Sx[S[[k]]]^4)/(((Sy[S[[k]]]^2)+((((tet_v^2)-var_v))*(Sx[S[[k]]]^2)))^2))}
  Vb_s <- function(k){(eta_s[k]+ci_s[k]+xi_s[k])+((eta_s[k]^2)*(eta_v+ci_v)/(eta_v^2))}
  Vb_s.fisc <- function(k){(eta_s[k])+((eta_s[k]^2)*(eta_v)/(eta_v^2))}
  b_s <- sapply(1:K,b_s); xi_s <- sapply(1:K,xi_s); Vb_s <- sapply(1:K,Vb_s); Vb_s.fisc <- sapply(1:K,Vb_s.fisc)
  
  # AMSE estimator
  amse <- function(k){max((((b_s[k]/(eta_v+eta_s[k]))^2)-((1/((eta_v+eta_s[k])^2))*Vb_s[k])),0)+((1/(eta_v+eta_s[k]))+((ci_v+ci_s[k])/((eta_v+eta_s[k])^2)))}
  amse <- sapply(1:K,amse)
  
  # Post-selection estimator
  sel <- which.min(c(var_v,amse))
  tet_ps <- c(tet_v,tet_s)[sel]
  if(sel==1){sel0 <- 1} else {sel0 <- (1-ifelse(sel==1,1,0))*(1+K0[(sel-1)])}
  
  # Covariance estimator
  E11 <- var_v
  E21 <- function(k){(eta_v+ci_v)/(eta_v*(eta_v+eta_s[k]))}; E21 <- sapply(1:K,E21)
  E31 <- function(k){-eta_s[k]*(eta_v+ci_v)/((eta_v^2)*(eta_v+eta_s[k]))}; E31 <- sapply(1:K,E31)
  E22 <- matrix(,nrow=K,ncol=K)
  for (k in 1:K){for (l in 1:K){
    if(sum(S[[k]]%in%S[[l]])>0){E22[k,l] <- (eta_v+ci_v+eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))])/((eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    if(sum(S[[k]]%in%S[[l]])==0){E22[k,l] <- (eta_v+ci_v)/((eta_v+eta_s[k])*(eta_v+eta_s[l]))}
  }}
  E32 <- matrix(,nrow=K,ncol=K)
  for (k in 1:K){for (l in 1:K){
    if(sum(S[[k]]%in%S[[l]])>0){E32[k,l] <- ((eta_v*(eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]))-(eta_s[k]*(eta_v+ci_v)))/(eta_v*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    if(sum(S[[k]]%in%S[[l]])==0){E32[k,l] <- -(eta_s[k]*(eta_v+ci_v))/(eta_v*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
  }}
  E33 <- matrix(,nrow=K,ncol=K)
  for (k in 1:K){for (l in 1:K){
    if(sum(S[[k]]%in%S[[l]])>0){E33[k,l] <- (((eta_v^2)*(eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+xi_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]))+((eta_s[k]*eta_s[l])*(eta_v+ci_v)))/((eta_v^2)*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    if(sum(S[[k]]%in%S[[l]])==0){E33[k,l] <- ((eta_s[k]*eta_s[l])*(eta_v+ci_v))/((eta_v^2)*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
  }}
  E1 <- c(E11,E21,E31); E2 <- rbind(t(E21),E22,E32); E3 <- rbind(t(E31),t(E32),E33); E <- cbind(E1,E2,E3)
  makeSymm <- function(m){
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  E <- makeSymm(E); colnames(E) <- NULL
  
  # confidence intervals
  var_ests <- function(k){(1/(eta_v+eta_s[k]))+((ci_v+ci_s[k])/((eta_v+eta_s[k])^2))}; var_ests <- sapply(1:K,var_ests)
  naive.var <- c(var_v,var_ests)[sel]
  
  R=5e2; R1 <- 4e2
  suppressWarnings(U <- rmvnorm(R,rep(0,(2*K)+1),E))
  bs <- function(k){b_s[k]/(eta_v+eta_s[k])}; bs <- sapply(1:K,bs)
  
  # (1-alpha/2) x 100% bias confidence region 
  b_int <- matrix(,nrow=K,ncol=2)
  for (k in 1:K){
    b_int[k,1] <- bs[k]-quantile(U[,1+K+k],1 - alpha/4)[[1]]
    b_int[k,2] <- bs[k]-quantile(U[,1+K+k], alpha/4)[[1]]
  }
  B0 <- matrix(NA,nrow=K,ncol=(R1-1))
  for (k in 1:K){B0[k,] <- runif((R1-1),b_int[k,1],b_int[k,2])}
  B0 <- cbind(B0,bs); colnames(B0)<-NULL 
  
  # sampling from distribution of tet_ps - tet0
  A0 <- function(b){
    A <- matrix(,nrow=R,ncol=K); for (k in 1:K){A[,k] <- pmax(((U[,1+K+k]+b[k])^2) - E33[k,k],0) + E22[k,k]}
    selA <- function(r){which.min(c(E11,A[r,]))}; selA <- sapply(1:R,selA)
    A1 <- vector(,length=R)
    for (r in 1:R){
      C0 <- vector(,length=K)
      for (k in 1:K){C0[k] <- (ifelse(selA[r]==(k+1),1,0)*(U[r,(k+1)]+b[k]))}
      A1[r] <- (ifelse(selA[r]==1,1,0)*U[r,1])+sum(C0)
    }
    return(A1)
  }
  
  # one-step interval
  onestep.ci <- c(tet_ps-quantile(A0(bs),1-alpha/2)[[1]],tet_ps-quantile(A0(bs),alpha/2)[[1]])
  
  # two-step interval 
  A2 <- matrix(,nrow=R,ncol=R1) # each column corresponds to a particular bias level
  for (r1 in 1:R1){
    A2[,r1] <- A0(B0[,r1])
  }
  
  qA2 <- matrix(,nrow=2,ncol=R1)
  for (r1 in 1:R1){
    qA2[1,r1] <- quantile(A2[,r1],alpha/4)[[1]]
    qA2[2,r1] <- quantile(A2[,r1],1 - alpha/4)[[1]]
  }
  twostep.ci <- c(tet_ps-max(qA2[2,]),tet_ps-min(qA2[1,]))
  
  onestep.min <- function(l0){
  # (1-delta) x 100% bias confidence region 
  b_int <- matrix(,nrow=K,ncol=2)
  for (k in 1:K){
    b_int[k,1] <- bs[k]-quantile(U[,1+K+k],1 - delta.opts[l0]/2)[[1]]
    b_int[k,2] <- bs[k]-quantile(U[,1+K+k], delta.opts[l0]/2)[[1]]
  }
  B0 <- matrix(NA,nrow=K,ncol=(R1-1))
  for (k in 1:K){B0[k,] <- runif((R1-1),b_int[k,1],b_int[k,2])}
  B0 <- cbind(B0,bs); colnames(B0)<-NULL 
  
  # sampling from distribution of tet_ps - tet0
  A0 <- function(b){
    A <- matrix(,nrow=R,ncol=K); for (k in 1:K){A[,k] <- pmax(((U[,1+K+k]+b[k])^2) - E33[k,k],0) + E22[k,k]}
    selA <- function(r){which.min(c(E11,A[r,]))}; selA <- sapply(1:R,selA)
    A1 <- vector(,length=R)
    for (r in 1:R){
      C0 <- vector(,length=K)
      for (k in 1:K){C0[k] <- (ifelse(selA[r]==(k+1),1,0)*(U[r,(k+1)]+b[k]))}
      A1[r] <- (ifelse(selA[r]==1,1,0)*U[r,1])+sum(C0)
    }
    return(A1)
  }
  
  A2 <- matrix(,nrow=R,ncol=R1) # each column corresponds to a particular bias level
  for (r1 in 1:R1){
      A2[,r1] <- A0(B0[,r1])
  }
  
 # shortest one-step interval according to a minimum coverage constraint
  qA2 <- matrix(,nrow=2,ncol=R1)
  for (r1 in 1:R1){
    qA2[1,r1] <- quantile(A2[,r1],alpha.opts[l0]/2)[[1]]
    qA2[2,r1] <- quantile(A2[,r1],1 - alpha.opts[l0]/2)[[1]]
  }

  # tip: use ecdf(A2[,r1])(qA2[2,]) not ecdf(A2[,r1])(qA2[2,k])
  coverage <- matrix(,nrow=R1,ncol=R1)
  for (r1 in 1:R1){
      coverage[,r1] <- ecdf(A2[,r1])(qA2[2,])-ecdf(A2[,r1])(qA2[1,])
  }
  res.cov <- function(k){sum(coverage[k,] >= 1-alpha.opts[l0]-gamma)==R1}; res.cov <- sapply(1:R1,res.cov)
  if(sum(!res.cov)==R1){onestep.min.ci <- c(NA,NA)} else{
  sel.cov <- which(res.cov)
  sel.cov1 <- which.min(qA2[2,sel.cov]-qA2[1,sel.cov]) # need to match new sel.cov onto old 1:400 indices
  sel.cov <- sel.cov[sel.cov1]
  onestep.min.ci <- c(tet_ps-qA2[2,sel.cov],tet_ps-qA2[1,sel.cov])
  }
  return(onestep.min.ci)
  }
  onestep.min <- sapply(1:length(alpha.opts),onestep.min)
  onestep.min <- c(onestep.min[1,which.min(onestep.min[2,]-onestep.min[1,])],onestep.min[2,which.min(onestep.min[2,]-onestep.min[1,])])
  res.list <- list("valid"=tet_v, "valid.se"=var_v, "valid.ci"=c(tet_v-qnorm(1 - alpha/2)*sqrt(var_v),tet_v+qnorm(1 - alpha/2)*sqrt(var_v)), "focused"=tet_ps, "onestep.ci"=onestep.ci, "onestep.min.ci"=onestep.min, "twostep.ci"=twostep.ci, "additional.ivs"=((ifelse(sel==1,1,0)*0)+(ifelse(sel>1,1,0)*length(S[[which.min(amse)]]))), "valid.CP"=CP_V <- sum(pmax((Bx[V]^2)-(Sx[V]^2),0)/(Sx[V]^2))/length(V), "focused.CP"= sum(pmax((Bx[S[[which.min(amse)]]]^2)-(Sx[S[[which.min(amse)]]]^2),0)/(Sx[S[[which.min(amse)]]]^2))/length(S[[which.min(amse)]]), "winner"=sel0, "naive.var"=naive.var, "default"=0)
 return(res.list)
}
