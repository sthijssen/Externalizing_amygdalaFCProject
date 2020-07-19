
################
### REVISION ###
################

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

library(irr)
library(MASS)
library(glmm)
library(nlme)
library(lme4)
library(foreign)
library(ggplot2)
library(dplyr)
library(lattice)
library(afex)
library(effects)
library(mgcv)
library(gamm4)
library(psy)
library(car)
library(mlmRev)
library(guide_legend)
library(psych)

options(scipen=99)

dataset.ext = read.spss("Z:/X.sav", use.missings=TRUE, to.data.frame=TRUE)


######### SENSITIVITY ANALYSIS LONGER SCANS #######

###ACC

#ACCs1

ACCs1Ext_LS <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)
ACCs1ExtDelta_LS <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(ACCs1Ext_LS)
stdCoef.merMod(ACCs1Ext_LS)

summary(ACCs1ExtDelta_LS)
stdCoef.merMod(ACCs1ExtDelta_LS)

#ACCs3

ACCs3Ext_LS <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)
ACCs3ExtDelta_LS <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(ACCs3Ext_LS)
stdCoef.merMod(ACCs3Ext_LS)

summary(ACCs3ExtDelta_LS)
stdCoef.merMod(ACCs3ExtDelta_LS)

#ACCs5

ACCs5Ext_LS <- lmer(biACCs5_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(ACCs5Ext_LS)
stdCoef.merMod(ACCs5Ext_LS)

#ACCs7

ACCs7Ext_LS <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)
ACCs7ExtDelta_LS <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(ACCs7Ext_LS)
stdCoef.merMod(ACCs7Ext_LS)

summary(ACCs7ExtDelta_LS)
stdCoef.merMod(ACCs7ExtDelta_LS)

#ACCi9

ACCi9Ext_LS <- lmer(biACCi9_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(ACCi9Ext_LS)
stdCoef.merMod(ACCi9Ext_LS)

## OFC

#OFC A

OFCAExt_LS <- lmer(biOFCA_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(OFCAExt_LS)
stdCoef.merMod(OFCAExt_LS)

#OFC I

OFCIExt_LS <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)
OFCIExtDelta_LS <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(OFCIExt_LS)
stdCoef.merMod(OFCIExt_LS)

summary(OFCIExtDelta_LS)
stdCoef.merMod(OFCIExtDelta_LS)

#OFC L

OFCLExt_LS <- lmer(biOFCL_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(OFCLExt_LS)
stdCoef.merMod(OFCLExt_LS)

#OFC M

OFCMExt_LS <- lmer(biOFCM_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.ext)

summary(OFCMExt_LS)
stdCoef.merMod(OFCMExt_LS)

#OFC P

OFCPExt_LS <- lmer(biOFCP_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (-1 | ID), data=dataset.ext)

summary(OFCPExt_LS)
stdCoef.merMod(OFCPExt_LS)


########### ICC ########## choose single_fixed_raters
## ACC

biACCs1_resUpgrade = read.spss("Z:/Ext_project200604_ACCs1_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biACCs1_resUpgrade)

biACCs3_resUpgrade = read.spss("Z:/Ext_project200604_ACCs3_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biACCs3_resUpgrade)

biACCs5_resUpgrade = read.spss("Z:/Ext_project200604_ACCs5_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biACCs5_resUpgrade)

biACCs7_resUpgrade = read.spss("Z:/Ext_project200604_ACCs7_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biACCs7_resUpgrade)

biACCi9_resUpgrade = read.spss("Z:/Ext_project200604_ACCi9_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biACCi9_resUpgrade)

## OFC
biOFCA_resUpgrade = read.spss("Z:/Ext_project200604_OFCA_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biOFCA_resUpgrade)

biOFCI_resUpgrade = read.spss("Z:/Ext_project200604_OFCI_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biOFCI_resUpgrade)

biOFCL_resUpgrade = read.spss("Z:/Ext_project200604_OFCL_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biOFCL_resUpgrade)

biOFCM_resUpgrade = read.spss("Z:/Ext_project200604_OFCM_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biOFCM_resUpgrade)

biOFCP_resUpgrade = read.spss("Z:/Ext_project200604_OFCP_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(biOFCP_resUpgrade)


####### ICC EXTERNALIZING ##########

EXT = read.spss("Z:/Ext_project200710_EXT_wide.sav", use.missings=TRUE, to.data.frame=TRUE)
ICC(EXT)

