#####################################################################
### ANALYSES PAPER EXTERNALIZING BEHAVIOR AND ACC/OFC-AMYGDALA FC ###
#################### ANTERIOR CINGULATE CORTEX ######################
#####################################################################

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

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

dataset.Ext = read.spss("Z:/X.sav", use.missings=TRUE, to.data.frame=TRUE)


#### INSPECTION RAW DATA

##ACC

ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs1_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biACCs1_resUpgrade, group=ID)) +
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")  +
  labs(title="Amygdala-caudal ACC functional connectivity", y="Z value", x="Age")  


ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs3_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biACCs3_resUpgrade, group=ID)) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")  +
  labs(title="Amygdala-dorsal ACC functional connectivity", y="Z value", x="Age")  
    

ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs5_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biACCs5_resUpgrade, group=ID)) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")+
  labs(title="Amygdala-rostral ACC functional connectivity", y="Z value", x="Age") 

ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs7_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biACCs7_resUpgrade, group=ID)) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")+
    labs(title="Amygdala-perigenual ACC functional connectivity", y="Z value", x="Age")   

ggplot(dataset.Ext, aes(x=AgeYear, y=biACCi9_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biACCi9_resUpgrade, group=ID)) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")  +
 labs(title="Amygdala-subgenual ACC functional connectivity", y="Z value", x="Age")  
    

#### ACC baseline model
##ACCs1: ACCs1UnGrNoSl lowest BIC

dataset.Ext <- dataset.Ext %>%
  mutate(AgeSq = AgeYear*AgeYear)

ACCs1UnGr <- lmer(biACCs1_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
ACCs1UnGrNoSl <- lmer(biACCs1_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
ACCs1UnGrNoI <- lmer(biACCs1_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
ACCs1UnGrNoSlSq <- lmer(biACCs1_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)

ACCs1UnconGrNoIAr1 <- lme(biACCs1_resUpgrade ~ 1   + AgeYear11, 
                          random = ~ 1 | ID, 
                          data = dataset.Ext, 
                          correlation=corAR1()
)

BIC(ACCs1UnGrNoSl, ACCs1UnGrNoI, ACCs1UnGr, ACCs1UnGrNoSlSq, ACCs1UnconGrNoIAr1)
AIC(ACCs1UnGrNoSl, ACCs1UnGrNoI, ACCs1UnGr, ACCs1UnGrNoSlSq, ACCs1UnconGrNoIAr1)

summary(ACCs1UnGrNoSl)
summary(ACCs1UnGrNoSlSq)

##ACCs3: ACCs3UnGrNoSl lowest BIC

ACCs3UnGr <- lmer(biACCs3_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
ACCs3UnGrNoSl <- lmer(biACCs3_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
ACCs3UnGrNoI <- lmer(biACCs3_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
ACCs3UnGrNoSlSq <- lmer(biACCs3_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


ACCs3UnconGrNoIAr1 <- lme(biACCs3_resUpgrade ~ 1   + AgeYear11, 
                          random = ~ 1 | ID, 
                          data = dataset.Ext, 
                          correlation=corAR1()
)

BIC(ACCs3UnGrNoSl, ACCs3UnGrNoI, ACCs3UnGr, ACCs3UnGrNoSlSq, ACCs3UnconGrNoIAr1)
AIC(ACCs3UnGrNoSl, ACCs3UnGrNoI, ACCs3UnGr, ACCs3UnGrNoSlSq, ACCs3UnconGrNoIAr1)

summary(ACCs3UnGrNoSl)
summary(ACCs1UnGrNoSlSq)

##ACCs5: ACCs3UnGrNoSl lowest BIC

ACCs5UnGr <- lmer(biACCs5_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
ACCs5UnGrNoSl <- lmer(biACCs5_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
ACCs5UnGrNoI <- lmer(biACCs5_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
ACCs5UnGrNoSlSq <- lmer(biACCs5_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


ACCs5UnconGrNoIAr1 <- lme(biACCs5_resUpgrade ~ 1   + AgeYear11, 
                          random = ~ 1 | ID, 
                          data = dataset.Ext, 
                          correlation=corAR1()
)

BIC(ACCs5UnGrNoSl, ACCs5UnGrNoI, ACCs5UnGr, ACCs5UnGrNoSlSq, ACCs5UnconGrNoIAr1)
AIC(ACCs5UnGrNoSl, ACCs5UnGrNoI, ACCs5UnGr, ACCs5UnGrNoSlSq, ACCs5UnconGrNoIAr1)

summary(ACCs5UnGrNoSl)
summary(ACCs5UnGrNoSlSq)

##ACCs7: ACCs7UnGrNoSl lowest BIC


summary(ACCs7UnGrNoSl)
summary(ACCs7UnGrNoSlSq)

##ACCi9: ACCi9UnGrNoSl lowest BIC

ACCi9UnGr <- lmer(biACCi9_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
ACCi9UnGrNoSl <- lmer(biACCi9_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
ACCi9UnGrNoI <- lmer(biACCi9_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
ACCi9UnGrNoSlSq <- lmer(biACCi9_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


ACCi9UnconGrNoIAr1 <- lme(biACCi9_resUpgrade ~ 1   + AgeYear11, 
                          random = ~ 1 | ID, 
                          data = dataset.Ext, 
                          correlation=corAR1()
)

BIC(ACCi9UnGrNoSl, ACCi9UnGrNoI, ACCi9UnGr, ACCi9UnGrNoSlSq, ACCi9UnconGrNoIAr1)
AIC(ACCi9UnGrNoSl, ACCi9UnGrNoI, ACCi9UnGr, ACCi9UnGrNoSlSq, ACCi9UnconGrNoIAr1)

summary(ACCi9UnGrNoSl)
summary(ACCi9UnGrNoI)
summary(ACCi9UnGrNoSlSq)

##### ACC conditional model
##ACCs1: Ext sign.


ACCs1Sex <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs1SexInt <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(ACCs1Sex, ACCs1SexInt)
AIC(ACCs1Sex, ACCs1SexInt)

summary(ACCs1Sex)
stdCoef.merMod(ACCs1Sex)

ACCs1Ext <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs1ExtInt <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs1ExtSexInt <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)

ACCs1ExtDelta <- lmer(biACCs1_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + (1 | ID), data=dataset.Ext)



BIC(ACCs1Ext, ACCs1ExtInt, ACCs1ExtSexInt)
AIC(ACCs1Ext, ACCs1ExtInt, ACCs1ExtSexInt)

summary(ACCs1Ext)
stdCoef.merMod(ACCs1Ext)

summary(ACCs1ExtDelta)
stdCoef.merMod(ACCs1ExtDelta)

ACCs1ExtDeltaAge <- lmer(biACCs1_resUpgrade ~ 1 + Age_T2 + DeltaAge + ZExt + ZExt * Age_T2 + ZExt * DeltaAge + DeltaAge * Age_T2 + ZExt * Age_T2 * DeltaAge + Sex + FD + (1 | ID), data=dataset.Ext)
summary(ACCs1ExtDeltaAge)


##ACCs3: Ext sign.

ACCs3Sex <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs3SexInt <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(ACCs3Sex, ACCs3SexInt)

summary(ACCs3Sex)
stdCoef.merMod(ACCs3Sex)

ACCs3Ext <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs3ExtInt <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs3ExtSexInt <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)

ACCs3Ext_new <- lmer(biACCs3_resUpgradeLongerScan ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs3Ext_new2 <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + LongerScan + (1 | ID), data=dataset.Ext)


BIC(ACCs3Ext, ACCs3ExtInt, ACCs3ExtSexInt)

summary(ACCs3Ext_new2)
stdCoef.merMod(ACCs3Ext)

ACCs3ExtDelta <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + (1 | ID), data=dataset.Ext)
summary(ACCs3ExtDelta)
stdCoef.merMod(ACCs3ExtDelta)

ACCs3ExtDeltaAge <- lmer(biACCs3_resUpgrade ~ 1 + Age_T2 + DeltaAge + ZExt + ZExt * Age_T2 + ZExt * DeltaAge + DeltaAge * Age_T2 + ZExt * Age_T2 * DeltaAge + Sex + FD + (1 | ID), data=dataset.Ext)
summary(ACCs3ExtDeltaAge)

summary(ACCs3ExtSexInt)
stdCoef.merMod(ACCs3ExtSexInt)

#since sex int is sgn. perform analyses for boys and girls seperately

sampleGirls <- subset(dataset.Ext, Gender > 1)
sampleGirls
ACCs3ExtGirls <- lmer(biACCs3_resUpgrade ~ 1 + ZAge_raw11 + ZExt + FD + (1 | ID), data=sampleGirls)

summary(ACCs3ExtGirls)
stdCoef.merMod(ACCs3ExtGirls)

sampleBoys <- subset(dataset.Ext, Gender < 2)
sampleBoys
ACCs3ExtBoys <- lmer(biACCs3_resUpgrade ~  + ZAge_raw11 + ZExt + FD + (1 | ID), data=sampleBoys)

summary(ACCs3ExtBoys)
stdCoef.merMod(ACCs3ExtBoys)



##ACCs5: ext marginally sign. 

ACCs5Sex <- lmer(biACCs5_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs5SexInt <- lmer(biACCs5_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(ACCs5Sex, ACCs5SexInt)

summary(ACCs5Sex)
stdCoef.merMod(ACCs5Sex)

ACCs5Ext <- lmer(biACCs5_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs5ExtInt <- lmer(biACCs5_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs5ExtSexInt <- lmer(biACCs5_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)


BIC(ACCs5Ext, ACCs5ExtInt, ACCs5ExtSexInt)

summary(ACCs5Ext)
stdCoef.merMod(ACCs5Ext)

ACCs5ExtDelta <- lmer(biACCs5_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + (1 | ID), data=dataset.Ext)
summary(ACCs5ExtDelta)

summary(ACCs5ExtSexInt)

ACCs5ExtDeltaAge <- lmer(biACCs5_resUpgrade ~ 1 + Age_T2 + DeltaAge + ZExt + ZExt * Age_T2 + ZExt * DeltaAge + DeltaAge * Age_T2 + ZExt * Age_T2 * DeltaAge + Sex + FD + (1 | ID), data=dataset.Ext)
summary(ACCs5ExtDeltaAge)


##ACCs7: Ext sign

ACCs7Sex <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs7SexInt <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

ACCs7AgeSqSex <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + AgeSq + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs7AgeSqSexInt <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + AgeSq* Sex + FD + (1 | ID), data=dataset.Ext)

BIC(ACCs7Sex, ACCs7SexInt)

summary(ACCs7Sex)
stdCoef.merMod(ACCs7Sex)

summary(ACCs7AgeSqSex)



ACCs7Ext <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs7ExtInt <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs7ExtRS <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 + Ext | ID), data=dataset.Ext)
ACCs7AgeSqExt <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11+ AgeSq + Ext + Sex + FD + (1 | ID), data=dataset.Ext)
ACCs7SexExtInt <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)


BIC(ACCs7Ext, ACCs7ExtInt, ACCs7ExtRS, ACCs7SexExtInt)

summary(ACCs7ExtRS)

summary(ACCs7Ext)
stdCoef.merMod(ACCs7Ext)
summary(ACCs7AgeSqExt)

ACCs7ExtDelta <- lmer(biACCs7_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + (1 | ID), data=dataset.Ext)
summary(ACCs7ExtDelta)
stdCoef.merMod(ACCs7ExtDelta)

ACCs7ExtDeltaAge <- lmer(biACCs7_resUpgrade ~ 1 + Age_T2 + DeltaAge + ZExt + ZExt * Age_T2 + ZExt * DeltaAge + DeltaAge * Age_T2 + ZExt * Age_T2 * DeltaAge + Sex + FD + (1 | ID), data=dataset.Ext)
summary(ACCs7ExtDeltaAge)

##ACCi9

ACCi9Sex <- lmer(biACCi9_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
ACCi9SexInt <- lmer(biACCi9_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(ACCi9Sex, ACCi9SexInt)


summary(ACCi9Sex)
stdCoef.merMod(ACCi9Sex)


ACCi9Ext <- lmer(biACCi9_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCi9ExtInt <- lmer(biACCi9_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
ACCi9ExtSexInt <- lmer(biACCi9_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)


BIC(ACCi9Ext, ACCi9ExtInt, ACCi9ExtSexInt)

summary(ACCi9Ext)
stdCoef.merMod(ACCi9Ext)

#### PLOTS ####

dataset.Ext <- dataset.Ext %>%
  mutate(predACCs1Ext = predict(ACCs1Ext))  %>%
  mutate(predACCs3Ext = predict(ACCs3Ext))  %>%
  mutate(predACCs5Ext = predict(ACCs5Ext))  %>%
  mutate(predACCs7Ext = predict(ACCs7Ext))  %>%
  mutate(predACCs3ExtSexInt = predict(ACCs3ExtSexInt))  


ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs7_resUpgrade, colour = NExt)) +
  geom_point(size=1) + theme_bw() +  
  geom_smooth(aes(y=predACCs7Ext), method="lm",se=T,size=2) +
  labs(title="Amygdala-perigenual ACC functional connectivity", y="Z value", x="Age")  + 
  scale_colour_discrete(name  ="Externalizing Behavior")


ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs5_resUpgrade, colour = NExt)) +
  geom_point(size=1) + theme_bw() +  
  geom_smooth(aes(y=predACCs5Ext), method="lm",se=T,size=2) +
  labs(title="Amygdala-rostral ACC functional connectivity", y="Z value", x="Age")  + 
  scale_colour_discrete(name  ="Externalizing Behavior")

ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs3_resUpgrade, colour = NExt)) +
  geom_point(size=1) + theme_bw() +  
  geom_smooth(aes(y=predACCs3Ext), method="lm",se=T,size=2) +
  labs(title="Amygdala-dorsal ACC functional connectivity", y="Z value", x="Age")  + 
  scale_colour_discrete(name  ="Externalizing Behavior")

ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs3_resUpgrade, colour = NExt)) +
  geom_point(size=1) + theme_bw() +  
  geom_smooth(aes(y=predACCs3ExtSexInt), method="lm",se=T,size=2) + facet_grid(. ~ Sex)+
  labs(title="Amygdala-dorsal ACC functional connectivity", y="Z value", x="Age")  + 
  scale_colour_discrete(name  ="Externalizing Behavior")


ggplot(dataset.Ext, aes(x=AgeYear, y=biACCs1_resUpgrade, colour = NExt)) +
  geom_point(size=1) + theme_bw() +  
  geom_smooth(aes(y=predACCs1Ext), method="lm",se=T,size=2) +
  labs(title="Amygdala-caudal ACC functional connectivity", y="Z value", x="Age")  + 
  scale_colour_discrete(name  ="Externalizing Behavior")

