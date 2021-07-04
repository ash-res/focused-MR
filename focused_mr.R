################################################################
### FOCUSED INSTRUMENT SELECTION FOR MENDELIAN RANDOMIZATION ###
################################################################

# estimation and inference at 95% asymptotic confidence level

# inputs: 
# bx = vector of p genetic associations with the exposure 
# sx = vector of p standard errors for genetic-exposure associations
# by = vector of p genetic associations with the outcome
# sy = vector of p standard errors for genetic-outcome associations
# V = indices of valid instruments
# tau (optional) : upper bound on the absolute value of pleiotropic effects

# output list:
# valid = Valid estimator using only the instruments V
# valid.se = standard error of the Valid estimator
# valid.ci = 95% asymptotic confidence intervals for the Valid estimator
# focused = Focused estimator which minimising estimated asymptotic MSE
# onestep.ci = 95% One-step asymptotic confidence intervals
# twostep.ci = 95% Two-step asymptotic confidence intervals
# twostepSR.ci = if tau is supplied, 95% Two-step SR asymptotic confidence intervals
# additional.ivs = number of additional IVs used by the Focused estimator
# valid.CP = concentration parameter of valid instruments V
# focused.CP = concentration parameter of additional IVs used by focused estimator
# clusters = number of clusters formed from additional IVs using NBClust


focused_mr <- function(bx,by,sx,sy,V,tau=NULL){
library(NbClust)
library(mvtnorm)
  
Bx <- c(bx[V],bx[-V])
By <- c(by[V],by[-V])
Sx <- c(sx[V],sx[-V])
Sy <- c(sy[V],sy[-V])

V <- 1:length(V)

## ADDITIONAL INSTRUMENT SETS BY K-MEANS CLUSTERING ##
sel <- NbClust(By[-V]/Bx[-V], method="kmeans",max.nc=8)$Best.partition
k0 <- max(sel)

Rs <- vector(,length=k0) 
for (k in 1:k0){
  Rs[k] <- length(combn(1:k0, k, simplify=FALSE))
}

Rs1 <- vector(,length=k0)
for (k in 1:k0){
  Rs1[k] <- sum(Rs[1:k])
}

SA <- list()
for (k in 1:Rs[1]){
  SA[[k]] <- combn(1:k0, 1, simplify=FALSE)[[k]]
}

for (l in 2:k0){
  for (k in 1:Rs[l]){
    SA[[k+Rs1[(l-1)]]] <- combn(1:k0, l, simplify=FALSE)[[k]]
  }
}

S <- list() 
for (s in 1:k0){
S[[s]] <- length(V)+which(sel==s)
}

for (s in (k0+1):length(SA)){
  S[[s]] <- unlist(S[SA[[s]]])
}

K=length(S)
rm(sel, SA, Rs, Rs1)


# VALID ESTIMATOR
Q <- function(tet){sum(((By[V]-(Bx[V]*tet))^2)/(Sy[V]^2 + (tet^2)*Sx[V]^2))}
init.val <- seq(-1,1,0.2)
Q.init <- vector(,length=length(init.val))
for(l in 1:length(init.val)){
  Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value
}
tet_v <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
eta_v <- sum(((Bx[V]^2)-Sx[V]^2)/((Sy[V]^2)+((tet_v^2)*(Sx[V]^2))))
ci_v <- sum((Sx[V]^2)*(Sy[V]^2)/(((Sy[V]^2)+((tet_v^2)*(Sx[V]^2)))^2))
var_v <- ((1/eta_v)+(ci_v/eta_v^2))

# ESTIMATION WITH ADDITIONAL INSTRUMENTS
tet_s <- vector(,length=K)
eta_s <- vector(,length=K)
ci_s <- vector(,length=K)
for (k in 1:K){
  Q <- function(tet){sum(((By[c(V,S[[k]])]-(Bx[c(V,S[[k]])]*tet))^2)/(Sy[c(V,S[[k]])]^2 + (tet^2)*Sx[c(V,S[[k]])]^2))}
  init.val <- seq(-0.8,1.2,0.2)
  Q.init <- vector(,length=length(init.val))
  for(l in 1:length(init.val)){
    Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value
  }
  tet_s[k] <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
  eta_s[k] <- sum(((Bx[S[[k]]]^2)-Sx[S[[k]]]^2)/((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2))))
  ci_s[k] <- sum((Sx[S[[k]]]^2)*(Sy[S[[k]]]^2)/(((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2)))^2))
}

if(missing(tau)) {
  # ESTIMATING THE BIASES
  b_s <- vector(,length=K)
  xi_s <- vector(,length=K)
  for (k in 1:K){
    b_s[k] <- sum((Bx[S[[k]]]*By[S[[k]]]-tet_v*((Bx[S[[k]]]^2)-(Sx[S[[k]]]^2)))/((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2))))
    xi_s[k] <- 2*(tet_v^2)*sum((Sx[S[[k]]]^4)/(((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2)))^2))
  }
  Vb_s <- vector(,length=K)
  for (k in 1:K){
    Vb_s[k] <- (eta_s[k]+ci_s[k]+xi_s[k])+((eta_s[k]^2)*(eta_v+ci_v)/(eta_v^2))
  }
  
  # AMSE ESTIMATION
  amse <- vector(,length=K)
  for (k in 1:K){
    amse[k] <- max((((b_s[k]/(eta_v+eta_s[k]))^2)-((1/((eta_v+eta_s[k])^2))*Vb_s[k])),0)+((1/(eta_v+eta_s[k]))+((ci_v+ci_s[k])/((eta_v+eta_s[k])^2)))
  }
  
  # POST-SELECTION ESTIMATOR
  sel <- which.min(c(var_v,amse))
  tet_ps <- c(tet_v,tet_s)[sel]
  
  
  # COVARIANCE ESTIMATION
  E11 <- var_v
  E21 <- vector(,length=K)
  E31 <- vector(,length=K)
  for (k in 1:K){
    E21[k] <- (eta_v+ci_v)/(eta_v*(eta_v+eta_s[k]))
    E31[k] <- -eta_s[k]*(eta_v+ci_v)/((eta_v^2)*(eta_v+eta_s[k]))
  }
  
  E22 <- matrix(,nrow=K,ncol=K)
  for (k in 1:K){
    for (l in 1:K){
      if(sum(S[[k]]%in%S[[l]])>0){E22[k,l] <- (eta_v+ci_v+eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))])/((eta_v+eta_s[k])*(eta_v+eta_s[l]))}
      if(sum(S[[k]]%in%S[[l]])==0){E22[k,l] <- (eta_v+ci_v)/((eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    }
  }
  
  
  E32 <- matrix(,nrow=K,ncol=K)
  for (k in 1:K){
    for (l in 1:K){
      if(sum(S[[k]]%in%S[[l]])>0){E32[k,l] <- ((eta_v*(eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]))-(eta_s[k]*(eta_v+ci_v)))/(eta_v*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
      if(sum(S[[k]]%in%S[[l]])==0){E32[k,l] <- -(eta_s[k]*(eta_v+ci_v))/(eta_v*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    }
  }
  
  
  E33 <- matrix(,nrow=K,ncol=K)
  for (k in 1:K){
    for (l in 1:K){
      if(sum(S[[k]]%in%S[[l]])>0){E33[k,l] <- (((eta_v^2)*(eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+xi_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]))+((eta_s[k]*eta_s[l])*(eta_v+ci_v)))/((eta_v^2)*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
      if(sum(S[[k]]%in%S[[l]])==0){E33[k,l] <- ((eta_s[k]*eta_s[l])*(eta_v+ci_v))/((eta_v^2)*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    }
  }
  
  E1 <- c(E11,E21,E31)
  E2 <- rbind(t(E21),E22,E32)
  E3 <- rbind(t(E31),t(E32),E33)
  E <- cbind(E1,E2,E3)
  makeSymm <- function(m){
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  E <- makeSymm(E)
  colnames(E) <- NULL
  
  
  # ONE-STEP INTERVALS
  R=5e2
  suppressWarnings(U <- rmvnorm(R,rep(0,(2*K)+1),E))
  A <- matrix(,nrow=R,ncol=K)
  for (k in 1:K){
    A[,k] <- ((U[,1+K+k]+(b_s[k]/(eta_v+eta_s[k])))^2) - E33[k,k] + E22[k,k]
  }
  
  selA <- vector(,length=R)
  for (r in 1:R){
    selA[r] <- which.min(c(E11,A[r,]))
  }
  
  A0 <- vector(,length=R)
  for (r in 1:R){
    C0 <- vector(,length=K)
    for (k in 1:K){
      C0[k] <- (ifelse(selA[r]==(k+1),1,0)*(U[r,(k+1)]+(b_s[k]/(eta_v+eta_s[k]))))
    }
    A0[r] <- (ifelse(selA[r]==1,1,0)*U[r,1])+sum(C0)
  }
  
  
  # TWO-STEP INTERVALS
  b_int <- matrix(,nrow=K,ncol=2)
  for (k in 1:K){
    b_int[k,1] <- (b_s[k]/(eta_v+eta_s[k]))-quantile(U[,1+K+k],0.9875)[[1]]
    b_int[k,2] <- (b_s[k]/(eta_v+eta_s[k]))-quantile(U[,1+K+k],0.0125)[[1]]
  }
  
  R2 <- 5e2
  CI_2S_L <- vector(,length=R2)
  CI_2S_U <- vector(,length=R2)
  for (r2 in 1:R2){
    B0 <- vector(,length=K)
    for (k in 1:K){
      B0[k] <- runif(1,b_int[k,1],b_int[k,2])
    }
    A1 <- vector(,length=R)
    for (r in 1:R){
      C1 <- vector(,length=K)
      for (k in 1:K){
        C1[k] <- (ifelse(selA[r]==(k+1),1,0)*(U[r,(k+1)]+B0[k]))
      }
      A1[r] <- (ifelse(selA[r]==1,1,0)*U[r,1])+sum(C1)
    }
    CI_2S_L[r2] <- tet_ps-quantile(A1,0.9875)[[1]]
    CI_2S_U[r2] <- tet_ps-quantile(A1,0.0125)[[1]]
  }
  
  res.list <- list("valid"=round(tet_v,3), "valid.se"=round(sqrt(var_v),3), "valid.ci"=round(c(tet_v-qnorm(0.975)*sqrt(var_v),tet_v+qnorm(0.975)*sqrt(var_v)),3), "focused"=round(tet_ps,3), "onestep.ci"=round(c(tet_ps-quantile(A0,0.975)[[1]],tet_ps-quantile(A0,0.025)[[1]]),3), "twostep.ci"=round(c(min(CI_2S_L),max(CI_2S_U)),3), "additional.ivs"=((ifelse(sel==1,1,0)*0)+(ifelse(sel>1,1,0)*length(S[[which.min(amse)]]))), "valid.CP"=round(CP_V <- (sum(((Bx[V]^2)-(Sx[V]^2))/(Sx[V]^2)))/length(V),0), "focused.CP"=round(CP_S <- (sum(((Bx[S[[which.min(amse)]]]^2)-(Sx[S[[which.min(amse)]]]^2))/(Sx[S[[which.min(amse)]]]^2)))/length(S[[which.min(amse)]]),0), "clusters"=k0)
}

else {
# ESTIMATING THE BIASES
b_s <- vector(,length=K)
xi_s <- vector(,length=K)
b_bl <- vector(,length=K)
b_bu <- vector(,length=K)
tau_b <- tau
for (k in 1:K){
  b_s[k] <- sum((Bx[S[[k]]]*By[S[[k]]]-tet_v*((Bx[S[[k]]]^2)-(Sx[S[[k]]]^2)))/((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2))))
  xi_s[k] <- 2*(tet_v^2)*sum((Sx[S[[k]]]^4)/(((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2)))^2))
  b_bl[k] <- ((-tau_b*sum(abs(Bx[S[[k]]])/((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2))))) - (tau_b*sqrt(2/pi)*sum(Sx[S[[k]]]/((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2))))))/(eta_v+eta_s[k])
  b_bu[k] <- ((tau_b*sum(abs(Bx[S[[k]]])/((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2))))) + (tau_b*sqrt(2/pi)*sum(Sx[S[[k]]]/((Sy[S[[k]]]^2)+((tet_v^2)*(Sx[S[[k]]]^2))))))/(eta_v+eta_s[k])
}
Vb_s <- vector(,length=K)
for (k in 1:K){
  Vb_s[k] <- (eta_s[k]+ci_s[k]+xi_s[k])+((eta_s[k]^2)*(eta_v+ci_v)/(eta_v^2))
}

# AMSE ESTIMATION
amse <- vector(,length=K)
for (k in 1:K){
  amse[k] <- max((((b_s[k]/(eta_v+eta_s[k]))^2)-((1/((eta_v+eta_s[k])^2))*Vb_s[k])),0)+((1/(eta_v+eta_s[k]))+((ci_v+ci_s[k])/((eta_v+eta_s[k])^2)))
}

# POST-SELECTION ESTIMATOR
sel <- which.min(c(var_v,amse))
tet_ps <- c(tet_v,tet_s)[sel]


# COVARIANCE ESTIMATION
E11 <- var_v
E21 <- vector(,length=K)
E31 <- vector(,length=K)
for (k in 1:K){
  E21[k] <- (eta_v+ci_v)/(eta_v*(eta_v+eta_s[k]))
  E31[k] <- -eta_s[k]*(eta_v+ci_v)/((eta_v^2)*(eta_v+eta_s[k]))
}

E22 <- matrix(,nrow=K,ncol=K)
for (k in 1:K){
  for (l in 1:K){
    if(sum(S[[k]]%in%S[[l]])>0){E22[k,l] <- (eta_v+ci_v+eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))])/((eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    if(sum(S[[k]]%in%S[[l]])==0){E22[k,l] <- (eta_v+ci_v)/((eta_v+eta_s[k])*(eta_v+eta_s[l]))}
  }
}


E32 <- matrix(,nrow=K,ncol=K)
for (k in 1:K){
  for (l in 1:K){
    if(sum(S[[k]]%in%S[[l]])>0){E32[k,l] <- ((eta_v*(eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]))-(eta_s[k]*(eta_v+ci_v)))/(eta_v*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    if(sum(S[[k]]%in%S[[l]])==0){E32[k,l] <- -(eta_s[k]*(eta_v+ci_v))/(eta_v*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
  }
}


E33 <- matrix(,nrow=K,ncol=K)
for (k in 1:K){
  for (l in 1:K){
    if(sum(S[[k]]%in%S[[l]])>0){E33[k,l] <- (((eta_v^2)*(eta_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+ci_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]+xi_s[(k*ifelse(length(S[[k]])<=length(S[[l]]),1,0))+(l*ifelse(length(S[[k]])>length(S[[l]]),1,0))]))+((eta_s[k]*eta_s[l])*(eta_v+ci_v)))/((eta_v^2)*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
    if(sum(S[[k]]%in%S[[l]])==0){E33[k,l] <- ((eta_s[k]*eta_s[l])*(eta_v+ci_v))/((eta_v^2)*(eta_v+eta_s[k])*(eta_v+eta_s[l]))}
  }
}

E1 <- c(E11,E21,E31)
E2 <- rbind(t(E21),E22,E32)
E3 <- rbind(t(E31),t(E32),E33)
E <- cbind(E1,E2,E3)
makeSymm <- function(m){
 m[upper.tri(m)] <- t(m)[upper.tri(m)]
 return(m)
}
E <- makeSymm(E)
colnames(E) <- NULL


# ONE-STEP INTERVALS
R=5e2
suppressWarnings(U <- rmvnorm(R,rep(0,(2*K)+1),E))
A <- matrix(,nrow=R,ncol=K)
for (k in 1:K){
  A[,k] <- ((U[,1+K+k]+(b_s[k]/(eta_v+eta_s[k])))^2) - E33[k,k] + E22[k,k]
}

selA <- vector(,length=R)
for (r in 1:R){
  selA[r] <- which.min(c(E11,A[r,]))
}

A0 <- vector(,length=R)
for (r in 1:R){
  C0 <- vector(,length=K)
  for (k in 1:K){
    C0[k] <- (ifelse(selA[r]==(k+1),1,0)*(U[r,(k+1)]+(b_s[k]/(eta_v+eta_s[k]))))
  }
  A0[r] <- (ifelse(selA[r]==1,1,0)*U[r,1])+sum(C0)
}


# TWO-STEP INTERVALS
b_int <- matrix(,nrow=K,ncol=2)
for (k in 1:K){
  b_int[k,1] <- (b_s[k]/(eta_v+eta_s[k]))-quantile(U[,1+K+k],0.9875)[[1]]
  b_int[k,2] <- (b_s[k]/(eta_v+eta_s[k]))-quantile(U[,1+K+k],0.0125)[[1]]
}

R2 <- 5e2
CI_2S_L <- vector(,length=R2)
CI_2S_U <- vector(,length=R2)
CI_2S_Lb <- vector(,length=R2)
CI_2S_Ub <- vector(,length=R2)
for (r2 in 1:R2){
  B0 <- vector(,length=K)
  B1 <- vector(,length=K)
  for (k in 1:K){
    B0[k] <- runif(1,b_int[k,1],b_int[k,2])
    if(min(c(b_int[k,2],b_bu[k]))<=max(c(b_int[k,1],b_bl[k]))){B1[k] <- runif(1,b_bl[k],b_bu[k])
    } else { B1[k] <- runif(1,max(c(b_int[k,1],b_bl[k])),min(c(b_int[k,2],b_bu[k])))}
  }
  A1 <- vector(,length=R)
  A2 <- vector(,length=R)
  for (r in 1:R){
    C1 <- vector(,length=K)
    C2 <- vector(,length=K)
    for (k in 1:K){
      C1[k] <- (ifelse(selA[r]==(k+1),1,0)*(U[r,(k+1)]+B0[k]))
      C2[k] <- (ifelse(selA[r]==(k+1),1,0)*(U[r,(k+1)]+B1[k]))
    }
    A1[r] <- (ifelse(selA[r]==1,1,0)*U[r,1])+sum(C1)
    A2[r] <- (ifelse(selA[r]==1,1,0)*U[r,1])+sum(C2)
  }
  CI_2S_L[r2] <- tet_ps-quantile(A1,0.9875)[[1]]
  CI_2S_U[r2] <- tet_ps-quantile(A1,0.0125)[[1]]
  CI_2S_Lb[r2] <- tet_ps-quantile(A2,0.9875)[[1]]
  CI_2S_Ub[r2] <- tet_ps-quantile(A2,0.0125)[[1]]
}

res.list <- list("valid"=round(tet_v,3), "valid.se"=round(sqrt(var_v),3), "valid.ci"=round(c(tet_v-qnorm(0.975)*sqrt(var_v),tet_v+qnorm(0.975)*sqrt(var_v)),3), "focused"=round(tet_ps,3), "onestep.ci"=round(c(tet_ps-quantile(A0,0.975)[[1]],tet_ps-quantile(A0,0.025)[[1]]),3), "twostep.ci"=round(c(min(CI_2S_L),max(CI_2S_U)),3), "twostepSR.ci"=round(c(min(CI_2S_Lb),max(CI_2S_Ub)),3), "additional.ivs"=((ifelse(sel==1,1,0)*0)+(ifelse(sel>1,1,0)*length(S[[which.min(amse)]]))), "valid.CP"=round(CP_V <- (sum(((Bx[V]^2)-(Sx[V]^2))/(Sx[V]^2)))/length(V),0), "focused.CP"=round(CP_S <- (sum(((Bx[S[[which.min(amse)]]]^2)-(Sx[S[[which.min(amse)]]]^2))/(Sx[S[[which.min(amse)]]]^2)))/length(S[[which.min(amse)]]),0), "clusters"=k0)
}

return(res.list)
}