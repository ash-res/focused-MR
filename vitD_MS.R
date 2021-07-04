rm(list=ls())
library(TwoSampleMR)
source("focused_mr.R")

exposure_dat <- extract_instruments("ebi-a-GCST90000615")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-b-18")
dat <- harmonise_data(exposure_dat, outcome_dat)

### SELECTING VALID IVS
sel_GC <- which(72607410-500000 <= dat$pos & dat$pos <= 72671237+500000 & dat$chr == 4)
sel_DHCR7 <- which(71145457-500000 <= dat$pos & dat$pos <= 71159439+500000  & dat$chr == 11)
sel_CYP2R1 <- which(14898986-500000 <= dat$pos & dat$pos <= 14913777+500000  & dat$chr == 11)
sel_CYP24A1 <- which(52769985-500000 <= dat$pos & dat$pos <= 52790525+500000  & dat$chr == 20)

V <- c(sel_GC,sel_DHCR7,sel_CYP2R1,sel_CYP24A1)
rm(sel_GC,sel_DHCR7,sel_CYP2R1,sel_CYP24A1)


# TWO-SAMPLE SUMMARY DATA ON GENETIC ASSOCIATIONS 
bx <- dat$beta.exposure
by <- dat$beta.outcome
sx <- dat$se.exposure
sy <- dat$se.outcome

# RUN FOCUSED-MR FUNCTION
focused_mr(bx,by,sx,sy,V)