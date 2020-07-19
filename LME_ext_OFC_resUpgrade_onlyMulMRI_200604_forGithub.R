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

#### inspection of raw data

ggplot(dataset.Ext, aes(x=AgeYear, y=Ext)) +
  geom_point(size=1) +
  geom_line(aes(y=Ext, group=ID)) + #facet_wrap(~ID,ncol=5) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "red") +
labs(title="Distribution of Externalizing Behavior", y="Externalizing behavior, percentage", x="Age")  + 
  scale_colour_discrete(name  ="Age")


##OFC

ggplot(dataset.Ext, aes(x=AgeYear, y=biOFCA_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biOFCA_resUpgrade, group=ID)) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green") +
  labs(title="Amygdala-anterior OFC functional connectivity", y="Z value", x="Age")  

ggplot(dataset.Ext, aes(x=AgeYear, y=biOFCI_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biOFCI_resUpgrade, group=ID)) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")  +
  labs(title="Amygdala-intermediate OFC functional connectivity", y="Z value", x="Age") 

ggplot(dataset.Ext, aes(x=AgeYear, y=biOFCL_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biOFCL_resUpgrade, group=ID)) + 
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")  +
  labs(title="Amygdala-lateral OFC functional connectivity", y="Z value", x="Age") 

ggplot(dataset.Ext, aes(x=AgeYear, y=biOFCM_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biOFCM_resUpgrade, group=ID)) +  
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green") +
  labs(title="Amygdala-medial OFC functional connectivity", y="Z value", x="Age")  

ggplot(dataset.Ext, aes(x=AgeYear, y=biOFCP_resUpgrade)) +
  geom_point(size=1) +
  geom_line(aes(y=biOFCP_resUpgrade, group=ID)) +  
  geom_smooth(method="lm", formula = y ~ x + I(x^2),se=T,size=2) +
  geom_smooth(method="lm",se=T,size=2, colour = "green")+
  labs(title="Amygdala-posterior OFC functional connectivity", y="Z value", x="Age")   


#### OFC baseline model

dataset.Ext <- dataset.Ext %>%
  mutate(AgeSq = AgeYear*AgeYear)

##OFC A: OFCAUnGrNoI lowest BIC (by less than 2%), but intercept explains more variance, so choose noSl

OFCAUnGr <- lmer(biOFCA_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
OFCAUnGrNoSl <- lmer(biOFCA_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
OFCAUnGrNoI <- lmer(biOFCA_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
OFCAUnGrNoSlSq <- lmer(biOFCA_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


OFCAUnconGrNoIAr1 <- lme(biOFCA_resUpgrade ~ 1   + AgeYear11, 
                         random = ~ 1 | ID, 
                         data = dataset.Ext, 
                         correlation=corAR1()
)

BIC(OFCAUnGrNoSl, OFCAUnGrNoI, OFCAUnGr, OFCAUnconGrNoIAr1, OFCAUnGrNoSlSq)
AIC(OFCAUnGrNoSl, OFCAUnGrNoI, OFCAUnGr, OFCAUnconGrNoIAr1, OFCAUnGrNoSlSq)

summary(OFCAUnGrNoSl)
summary(OFCAUnGrNoI)

stdCoef.merMod(OFCAUnGrNoSl)
##OFC P: OFCPUnGrNoI lowest BIC 
OFCPUnGr <- lmer(biOFCP_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
OFCPUnGrNoSl <- lmer(biOFCP_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
OFCPUnGrNoI <- lmer(biOFCP_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
OFCPUnGrNoSlSq <- lmer(biOFCP_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


OFCPUnconGrNoIAr1 <- lme(biOFCP_resUpgrade ~ 1   + AgeYear11, 
                         random = ~ 1 | ID, 
                         data = dataset.Ext, 
                         correlation=corAR1()
)

BIC(OFCPUnGrNoSl, OFCPUnGrNoI, OFCPUnGr, OFCPUnconGrNoIAr1, OFCPUnGrNoSlSq)
AIC(OFCPUnGrNoSl, OFCPUnGrNoI, OFCPUnGr, OFCPUnconGrNoIAr1, OFCPUnGrNoSlSq)

summary(OFCPUnGrNoI)
summary(OFCPUnGrNoSl)

stdCoef.merMod(OFCPUnGrNoI)

##OFC L: OFCAUnGrNoSl lowest BIC 
OFCLUnGr <- lmer(biOFCL_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
OFCLUnGrNoSl <- lmer(biOFCL_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
OFCLUnGrNoI <- lmer(biOFCL_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
OFCLUnGrNoSlSq <- lmer(biOFCL_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


OFCLUnconGrNoIAr1 <- lme(biOFCL_resUpgrade ~ 1   + AgeYear11, 
                         random = ~ 1 | ID, 
                         data = dataset.Ext, 
                         correlation=corAR1()
)

BIC(OFCLUnGrNoSl, OFCLUnGrNoI, OFCLUnGr, OFCLUnconGrNoIAr1, OFCLUnGrNoSlSq)
AIC(OFCLUnGrNoSl, OFCLUnGrNoI, OFCLUnGr, OFCLUnconGrNoIAr1, OFCLUnGrNoSlSq)

summary(OFCLUnGrNoSl)
summary(OFCLUnGrNoSlSq)

stdCoef.merMod(OFCLUnGrNoSl)


##OFC M: OFCAUnGrNoI lowest BIC (by 1), but RI explains more variance thus choose NoSl
OFCMUnGr <- lmer(biOFCM_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
OFCMUnGrNoSl <- lmer(biOFCM_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
OFCMUnGrNoI <- lmer(biOFCM_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
OFCMUnGrNoSlSq <- lmer(biOFCM_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


OFCMUnconGrNoIAr1 <- lme(biOFCM_resUpgrade ~ 1   + AgeYear11, 
                         random = ~ 1 | ID, 
                         data = dataset.Ext, 
                         correlation=corAR1()
)

BIC(OFCMUnGrNoSl, OFCMUnGrNoI, OFCMUnGr, OFCMUnconGrNoIAr1, OFCMUnGrNoSlSq)
AIC(OFCMUnGrNoSl, OFCMUnGrNoI, OFCMUnGr, OFCMUnconGrNoIAr1, OFCMUnGrNoSlSq)

summary(OFCMUnGrNoSl)
summary(OFCMUnGrNoI)

stdCoef.merMod(OFCMUnGrNoSl)


##OFC I: OFCAUnGrNoSl lowest BIC 
OFCIUnGr <- lmer(biOFCI_resUpgrade ~ 1 + AgeYear11  + (1 + AgeYear11 | ID), data=dataset.Ext)
OFCIUnGrNoSl <- lmer(biOFCI_resUpgrade ~ 1 + AgeYear11  + (1 | ID), data=dataset.Ext)
OFCIUnGrNoI <- lmer(biOFCI_resUpgrade ~ 1 + AgeYear11 + (-1+ AgeYear11 | ID), data=dataset.Ext)
OFCIUnGrNoSlSq <- lmer(biOFCI_resUpgrade ~ 1 + AgeYear11  + AgeSq + (1 | ID), data=dataset.Ext)


OFCIUnconGrNoIAr1 <- lme(biOFCI_resUpgrade ~ 1   + AgeYear11, 
                         random = ~ 1 | ID, 
                         data = dataset.Ext, 
                         correlation=corAR1()
)

BIC(OFCIUnGrNoSl, OFCIUnGrNoI, OFCIUnGr, OFCIUnconGrNoIAr1, OFCIUnGrNoSlSq)
AIC(OFCIUnGrNoSl, OFCIUnGrNoI, OFCIUnGr, OFCIUnconGrNoIAr1, OFCIUnGrNoSlSq)

summary(OFCIUnGrNoSl)
summary(OFCIUnGrNoSlSq)

stdCoef.merMod(OFCIUnGrNoSl)


##### OFC conditional model
##OFC A


OFCASex <- lmer(biOFCA_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
OFCASexInt <- lmer(biOFCA_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCASex, OFCASexInt)

summary(OFCASex)
stdCoef.merMod(OFCASex)


OFCAExt <- lmer(biOFCA_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCAExtInt <- lmer(biOFCA_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCAExtSexInt <- lmer(biOFCA_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCAExt, OFCAExtInt, OFCAExtSexInt)

summary(OFCAExt)
stdCoef.merMod(OFCAExt)


##OFC I: Ext sign; Int*Sex marginally sign. 


OFCISex <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
OFCISexInt <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCISex, OFCISexInt)

summary(OFCISex)
stdCoef.merMod(OFCISex)

OFCIExt <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCIExtInt <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11* ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCIExtSexInt <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11+ ZExt * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCIExt, OFCIExtInt, OFCIExtSexInt)

summary(OFCIExt)
stdCoef.merMod(OFCIExt)


OFCIExtDelta <- lmer(biOFCI_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + (1 | ID), data=dataset.Ext)
summary(OFCIExtDelta)
stdCoef.merMod(OFCIExtDelta)



 ###OFC L


OFCLSex <- lmer(biOFCL_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
OFCLSexInt <- lmer(biOFCL_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCLSex, OFCLSexInt)

summary(OFCLSex)
stdCoef.merMod(OFCLSex)

OFCLExt <- lmer(biOFCL_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCLExtInt <- lmer(biOFCL_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCLExtSexInt <- lmer(biOFCL_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCLExt, OFCLExtInt, OFCLExtSexInt)

summary(OFCLExt)
stdCoef.merMod(OFCLExt)



##OFC M: 


OFCMSex <- lmer(biOFCM_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (1 | ID), data=dataset.Ext)
OFCMSexInt <- lmer(biOFCM_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCMSex, OFCMSexInt)

summary(OFCMSex)
stdCoef.merMod(OFCMSex)

OFCMExt <- lmer(biOFCM_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCMExtInt <- lmer(biOFCM_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCMExtSexInt <- lmer(biOFCM_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCMExt, OFCMExtInt, OFCMExtSexInt)

summary(OFCMExt)
stdCoef.merMod(OFCMExt)




##OFC P: Ext sign


OFCPSex <- lmer(biOFCP_resUpgrade ~ 1 + ZAge_raw11 + Sex + FD + (-1 | ID), data=dataset.Ext)
OFCPSexInt <- lmer(biOFCP_resUpgrade ~ 1 + ZAge_raw11 * Sex + FD + (-1 | ID), data=dataset.Ext)

BIC(OFCPSex, OFCPSexInt)

summary(OFCPSex)
stdCoef.merMod(OFCPSex)

OFCPExt <- lmer(biOFCP_resUpgrade ~ 1 + ZAge_raw11 + ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCPExtInt <- lmer(biOFCP_resUpgrade ~ 1 + ZAge_raw11 * ZExt + Sex + FD + (1 | ID), data=dataset.Ext)
OFCPExtSexInt <- lmer(biOFCP_resUpgrade ~ 1 + ZAge_raw11 + ZExt * Sex + FD + (1 | ID), data=dataset.Ext)

BIC(OFCPExt, OFCPExtInt, OFCPExtSexInt)

summary(OFCPExt)
stdCoef.merMod(OFCPExt)

OFCPExtDelta <- lmer(biOFCP_resUpgrade ~ 1 + ZAge_raw11 + ExtT2 + ExtDelta + Sex + FD + (1 | ID), data=dataset.Ext)
summary(OFCPExtDelta)
stdCoef.merMod(OFCPExtDelta)


#### PLOTS ####

dataset.Ext <- dataset.Ext %>%

  mutate(predOFCPExt = predict(OFCPExt))  %>%
  mutate(predOFCIExt = predict(OFCIExt))



ggplot(dataset.Ext, aes(x=AgeYear, y=biOFCI_resUpgrade, colour = NExt)) +
  geom_point(size=1) + theme_bw() +  
  geom_smooth(aes(y=predOFCIExt), method="lm",se=T,size=2) +
  labs(title="Amygdala-intermediate OFC functional connectivity", y="Z value", x="Age")  + 
  scale_colour_discrete(name  ="Externalizing Behavior")

ggplot(dataset.Ext, aes(x=AgeYear, y=biOFCP_resUpgrade, colour = NExt)) +
  geom_point(size=1) + theme_bw() +  
  geom_smooth(aes(y=predOFCPExt), method="lm",se=T,size=2)+
  labs(title="Amygdala-posterior OFC functional connectivity", y="Z value", x="Age")  + 
  scale_colour_discrete(name  ="Externalizing Behavior")



