rm(list=ls())
library(TwoSampleMR); library(ieugwasr); library(mvtnorm); library(ggplot2); library(ggpubr); source("focusedMR.R")

results <- list(); p1 <- list(); df <- list(); results <- list()

# choose outcome
outcome.name <- c("coronary artery disease", "atrial fibrillation", "heart failure", "stroke (small vessels)", "stroke (large artery)", "stroke (cardioembolic)", "ischmeic stroke", "stroke (any)", "parkinson's disease", "multiple sclerosis", "alzheimer's disease", "ulcerative colitis", "rheumatoid arthritis", "primary biliary cirrhosis", "inflammatory bowel disease", "eczema", "crohn's disease", "asthma", "allergies", "squamous cell lung cancer", "lung cancer", "lung adenocarcinoma", "breast cancer", "type 2 diabetes", "osteoarthritis", "ADHD", "anorexia nervosa", "autism", "depressed affect", "depression (broad)", "depressive symptoms", "major depressive disorder (ICD-coded)", "major depressive disorder", "covid-19 (infection)", "covid-19 (hospitalisation)", "covid-19 (severe)")
outcome.id <-  c("ebi-a-GCST003116", "ebi-a-GCST006414", "ebi-a-GCST009541", "ebi-a-GCST006909", "ebi-a-GCST006907", "ebi-a-GCST006910", "ebi-a-GCST006908", "ebi-a-GCST006906", "ieu-b-7", "ieu-b-18", "ieu-b-2", "ieu-a-32", "ebi-a-GCST002318", "ebi-a-GCST003129", "ieu-a-31", "ieu-a-996", "ieu-a-30", "ebi-a-GCST006862", "ebi-a-GCST005038", "ieu-a-967", "ieu-a-966", "ieu-a-965", "ieu-a-1126", "finn-b-T2D", "ebi-a-GCST005814", "ieu-a-1183", "ieu-a-1186", "ieu-a-1185", "ebi-a-GCST006475", "ebi-a-GCST005902", "ieu-a-1000", "ebi-a-GCST005903", "ebi-a-GCST005904", "ebi-a-GCST011073" , "ebi-a-GCST011081", "ebi-a-GCST011075")

rm.outcome <- c("stroke (small vessels)", "stroke (large artery)", "stroke (cardioembolic)", "stroke (any)", "parkinson's disease", "ulcerative colitis", "inflammatory bowel disease", "crohn's disease", "allergies", "squamous cell lung cancer", "lung adenocarcinoma", "breast cancer", "ADHD", "autism", "depressed affect", "depression (broad)", "depressive symptoms", "major depressive disorder (ICD-coded)", "covid-19 (hospitalisation)", "covid-19 (severe)")
rm.outcome <- which (outcome.name %in% rm.outcome)
outcome.id <- outcome.id[-rm.outcome]
outcome.name <- outcome.name[-rm.outcome]

for (t in 1:length(outcome.id)){

# extract outcome data
exposure_dat <- extract_instruments("ieu-b-4812")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcome.id[t])
dat <- harmonise_data(exposure_dat, outcome_dat)
  
# selecting valid IVs
sel_GC <- which(72607410-500000 <= dat$pos & dat$pos <= 72671237+500000 & dat$chr == 4)
sel_DHCR7 <- which(71145457-500000 <= dat$pos & dat$pos <= 71159439+500000  & dat$chr == 11)
sel_CYP2R1 <- which(14898986-500000 <= dat$pos & dat$pos <= 14913777+500000  & dat$chr == 11)
sel_CYP24A1 <- which(52769985-500000 <= dat$pos & dat$pos <= 52790525+500000  & dat$chr == 20)
  
V <- c(sel_GC,sel_DHCR7,sel_CYP2R1,sel_CYP24A1)
rm(sel_GC,sel_DHCR7,sel_CYP2R1,sel_CYP24A1,exposure_dat,outcome_dat)

# two-sample summary data 
bx <- dat$beta.exposure; by <- dat$beta.outcome; sx <- dat$se.exposure; sy <- dat$se.outcome

# run focused-MR function
res <- focusedMR(bx,by,sx,sy,V,k0=3,alpha=0.05,gamma=0.2)

# valid heterogeneity
Q.V <- sum(((by[V]-(bx[V]*res$valid))^2)/(sy[V]^2 + (res$valid^2)*sx[V]^2))
J_test_pval_V <- 1-pchisq(Q.V,length(bx[V])-1) # passes heterogeneity test

# if we used all variants as instruments...
Q <- function(tet){sum(((by-(bx*tet))^2)/(sy^2 + (tet^2)*sx^2))}
init.val <- seq(-1,1,0.2)
Q.init <- vector(,length=length(init.val))
for(l in 1:length(init.val)){Q.init[l]<-optim(init.val[l], Q, method="Brent",lower=-1e2,upper=1e2)$value}
tet_est <- optim(init.val[which.min(Q.init)[[1]]], Q, method="Brent",lower=-1e2,upper=1e2)$par
eta_est <- sum(((bx^2)-sx^2)/((sy^2)+((tet_est^2)*(sx^2))))
ci_est <- sum((sx^2)*(sy^2)/(((sy^2)+((tet_est^2)*(sx^2)))^2))
var_est <- ((1/eta_est)+(ci_est/eta_est^2))
full_ci <- c(tet_est-qnorm(1 - 0.05/2)*sqrt(var_est),tet_est+qnorm(1 - 0.05/2)*sqrt(var_est))
full_pval <- 2*(1-pnorm(abs(tet_est)/sqrt(var_est)))
J_test_pval <- 1-pchisq(Q(tet_est),length(bx)-1)

# results
valid <- list("est"=res$valid, "ci"=res$valid.ci, "cp"=res$valid.CP, "number"=length(V), "het"=J_test_pval_V)
focused <- list("est"=res$focused, "onestep"=res$onestep.ci, "onestepmin"=res$onestep.min.ci, "twostep"=res$twostep.ci, "cp"=res$focused.CP, "number"=res$additional.ivs, "naive"=c(res$focused-qnorm(0.975)*sqrt(res$naive.var),res$focused+qnorm(0.975)*sqrt(res$naive.var)))
all <- list("est"=tet_est, "ci"=full_ci, "het"=J_test_pval, "number"=length(bx))
results[[t]] <- list()
results[[t]]$valid <- valid; results[[t]]$focused <- focused; results[[t]]$all <- all

## box plot ##
boxLabels = c("Full", "2-step", "1-step", "Focused", "Valid")

df[[t]] <- data.frame(yAxis = length(boxLabels):1, 
                 boxOdds = c(results[[t]]$all$est,results[[t]]$focused$est,results[[t]]$focused$est,results[[t]]$focused$est,results[[t]]$valid$est),
                 boxCILow = c(results[[t]]$all$ci[1],results[[t]]$focused$twostep[1],results[[t]]$focused$onestep[1],results[[t]]$focused$onestepmin[1],results[[t]]$valid$ci[1]), 
                 boxCIHigh = c(results[[t]]$all$ci[2],results[[t]]$focused$twostep[2],results[[t]]$focused$onestep[2],results[[t]]$focused$onestepmin[2],results[[t]]$valid$ci[2])
)


p1[[t]] <- ggplot(df[[t]], aes(x = boxOdds, y = boxLabels)) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 2, color = "orange") +
    scale_x_continuous(limits = c(min(0,min(df[[t]]$boxCILow)),max(0,max(df[[t]]$boxCIHigh)))) +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("") +
    ggtitle(paste0(outcome.name[t])) + 
    theme(plot.title = element_text(hjust = 0.5))
}

plot1 <- ggarrange(p1[[1]],NULL,p1[[2]],NULL,p1[[3]],NULL,p1[[4]],p1[[5]],NULL,p1[[6]],NULL,p1[[7]],NULL,p1[[8]],p1[[9]],NULL,p1[[10]],NULL,p1[[11]],NULL,p1[[12]],p1[[13]],NULL,p1[[14]],NULL,p1[[15]],NULL,p1[[16]],ncol=7,nrow=4,widths=c(1,-0.02,1,-0.02,1,-0.02,1,1,-0.02,1,-0.02,1,-0.02,1,1,-0.02,1,-0.02,1,-0.02,1,1,-0.02,1,-0.02,1,-0.02,1))

IVs.numbers <- function(k){c(results[[k]]$valid$number,results[[k]]$focused$number,results[[k]]$all$number-results[[k]]$valid$number)}
sapply(1:length(outcome.id),IVs.numbers)

CP <- function(k){c(results[[k]]$valid$cp,results[[k]]$focused$cp)}
sapply(1:length(outcome.id),CP)

plot1 <- annotate_figure(plot1,
                         bottom = text_grob((expression("Estimated causal effect: "*hat(theta)*" ")), color = "black", size = 16),
)

ggsave(filename="vitaminD.png", plot=plot1, device="png",width=15,height=8,dpi=300)