rm(list=ls())
library(TwoSampleMR)
source("focused_mr.R")

exposure_dat <- extract_instruments("ieu-b-38")
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-297")
dat <- harmonise_data(exposure_dat, outcome_dat)


### SELECTING VALID IVS FROM GENE REGIONS CORRESPONDING TO ANTIHYPERTENSIVE DRUGS
sel_ACE <- which(61554422-500000 <= dat$pos & dat$pos <= 61575734+500000 & dat$chr == 17)
sel_AGTR1 <- which(148415690-500000 <= dat$pos & dat$pos <= 148460790+500000  & dat$chr == 3)
sel_ADRB1 <- which(115803625-500000 <= dat$pos & dat$pos <= 115806663+500000  & dat$chr == 10)
sel_SLC12A3 <- which(56899119-500000 <= dat$pos & dat$pos <= 56949762+500000  & dat$chr == 16)
sel_CACNA1D<- which(53528638-500000 <= dat$pos & dat$pos <= 53847760+500000  & dat$chr == 3)
sel_CACNA2D1<- which(81575760-500000 <= dat$pos & dat$pos <= 82073272+500000  & dat$chr == 7)
sel_CACNA2D2<- which(50400044-500000 <= dat$pos & dat$pos <= 50541675+500000  & dat$chr == 3)
sel_CACNA1S<- which(201008640-500000 <= dat$pos & dat$pos <= 201081554+500000  & dat$chr == 1)
sel_CACNB1<- which(37329706-500000 <= dat$pos & dat$pos <= 37353922+500000  & dat$chr == 17)
sel_CACNB2<- which(18429353-500000 <= dat$pos & dat$pos <= 18832486+500000  & dat$chr == 10)
sel_CACNB3<- which(49208263-500000 <= dat$pos & dat$pos <= 49222724+500000  & dat$chr == 12)
sel_CACNB4<- which(152689285-500000 <= dat$pos & dat$pos <= 152955681+500000  & dat$chr == 2)
sel_CACNG1 <- which(65040670-500000 <= dat$pos & dat$pos <= 65052913+500000  & dat$chr == 17)
sel_CACNA1C <-which(2162153-500000 <= dat$pos & dat$pos <= 2807116+500000  & dat$chr == 12)

V <- c(sel_ACE,sel_AGTR1,sel_ADRB1,sel_SLC12A3,sel_CACNA1D,sel_CACNA2D1,sel_CACNA2D2,sel_CACNA1S,sel_CACNB1,sel_CACNB2,sel_CACNB3,sel_CACNB4,sel_CACNG1,sel_CACNA1C)
rm(sel_ACE,sel_AGTR1,sel_ADRB1,sel_SLC12A3,sel_CACNA1D,sel_CACNA2D1,sel_CACNA2D2,sel_CACNA1S,sel_CACNB1,sel_CACNB2,sel_CACNB3,sel_CACNB4,sel_CACNG1,sel_CACNA1C)


# TWO-SAMPLE SUMMARY DATA ON GENETIC ASSOCIATIONS 
bx <- dat$beta.exposure
by <- dat$beta.outcome
sx <- dat$se.exposure
sy <- dat$se.outcome

# RUN FOCUSED-MR FUNCTION
focused_mr(bx,by,sx,sy,V)


