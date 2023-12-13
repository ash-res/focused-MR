# Vitamin D example

rm(list=ls())
library(TwoSampleMR); library(ieugwasr); library(mvtnorm); source("focusedMR.R")

# choose outcome
outcome.name <- c("coronary artery disease", "atrial fibrillation", "heart failure", "ischemic stroke", "multiple sclerosis", "alzheimer's disease", "rheumatoid arthritis", "primary biliary cirrhosis", "eczema", "asthma", "lung cancer", "type 2 diabetes", "osteoarthritis", "anorexia nervosa", "major depressive disorder", "covid-19 infection")
outcome.id <-  c("ebi-a-GCST003116", "ebi-a-GCST006414", "ebi-a-GCST009541", "ebi-a-GCST006908", "ieu-b-18", "ieu-b-2", "ebi-a-GCST002318", "ebi-a-GCST003129", "ieu-a-996", "ebi-a-GCST006862", "ieu-a-966", "finn-b-T2D", "ebi-a-GCST005814", "ieu-a-1186", "ebi-a-GCST005904", "ebi-a-GCST011073")

# repeat the following for any of the outcomes listed above
t = 1 # here we choose the first outcome: CAD risk

# extract outcome data
exposure_dat <- extract_instruments("ieu-b-4812") # ID for vitamin D
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcome.id[t])
dat <- harmonise_data(exposure_dat, outcome_dat)
  
# selecting core IVs
sel_GC <- which(72607410-500000 <= dat$pos & dat$pos <= 72671237+500000 & dat$chr == 4)
sel_DHCR7 <- which(71145457-500000 <= dat$pos & dat$pos <= 71159439+500000  & dat$chr == 11)
sel_CYP2R1 <- which(14898986-500000 <= dat$pos & dat$pos <= 14913777+500000  & dat$chr == 11)
sel_CYP24A1 <- which(52769985-500000 <= dat$pos & dat$pos <= 52790525+500000  & dat$chr == 20)
  
S0 <- c(sel_GC,sel_DHCR7,sel_CYP2R1,sel_CYP24A1)
rm(sel_GC,sel_DHCR7,sel_CYP2R1,sel_CYP24A1,exposure_dat,outcome_dat)

# two-sample summary data 
bx <- dat$beta.exposure; by <- dat$beta.outcome; sx <- dat$se.exposure; sy <- dat$se.outcome

# run focused-MR function
res <- focusedMR(bx,by,sx,sy,S0,k0=3,alpha=0.05,gamma=0.2)

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
