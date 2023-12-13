# drug target PCSK9 gene example

rm(list=ls())
library(TwoSampleMR); library(ieugwasr); library(mvtnorm); source("focusedMR.R")

# choose outcome
outcome.name <- c("coronary artery disease", "atrial fibrillation", "heart failure", "stroke (small vessels)", "stroke (large artery)", "stroke (cardioembolic)", "ischemic stroke", "stroke (any)","alzheimer's disease", "crohn's disease", "inflammatory bowel disease", "multiple sclerosis", "type 2 diabetes", "parkinson's disease", "lung cancer", "rheumatoid arthritis")
outcome.id <-  c("ebi-a-GCST003116", "ebi-a-GCST006414", "ebi-a-GCST009541", "ebi-a-GCST006909", "ebi-a-GCST006907", "ebi-a-GCST006910", "ebi-a-GCST006908", "ebi-a-GCST006906","ieu-b-2", "ieu-a-30", "ieu-a-31", "ieu-b-18", "finn-b-T2D", "ieu-b-7", "ieu-a-966", "ebi-a-GCST002318")

# repeat the following for any of the outcomes listed above
t = 1 # here we choose the first outcome: CAD risk

# extract outcome data
load("PCSK9.Rdata")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcome.id[t])
dat0 <- harmonise_data(exposure_dat, outcome_dat); rm(exposure_dat,outcome_dat)

# two-sample summary data 
bx <- dat0$beta.exposure; by <- dat0$beta.outcome; sx <- dat0$se.exposure; sy <- dat0$se.outcome
S0 <- which(55505221 <= dat0$pos & dat0$pos <= 55530525)

# run focused-MR function
res <- focusedMR(bx,by,sx,sy,S0,k0=1,alpha=0.05,gamma=0.2)

# the results are...
res

# if we used all variants as instruments...
Q <- function(tet){sum(((by-(bx*tet))^2)/(sy^2 + (tet^2)*sx^2))}
init.val <- seq(-1,1,0.2)
Q.init <- vector(,length=length(init.val))
for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value}
tet_est <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par # full estimator
eta_est <- sum(((bx^2)-sx^2)/((sy^2)+((tet_est^2)*(sx^2)))) 
ci_est <- sum((sx^2)*(sy^2)/(((sy^2)+((tet_est^2)*(sx^2)))^2))
var_est <- ((1/eta_est)+(ci_est/eta_est^2)) # variance of the full estimator
full_ci <- c(tet_est-qnorm(1 - 0.05/2)*sqrt(var_est),tet_est+qnorm(1 - 0.05/2)*sqrt(var_est)) # 95% confidence intervals for the full estimator 
full_pval <- 2*(1-pnorm(abs(tet_est)/sqrt(var_est))) # p-value corresponding to the full estimator

