library(sandwich)
library(lmtest)
library(lme4)
library(lmerTest)
library(nlme)
library(CR2)
library(brms)
dat <- read.delim("~/Desktop/Documents/R/cleandata.txt")
dat <-  subset(dat, SoundType != "FAN")
dat$SoundType <- as.factor(dat$SoundType)
dat$Age <- as.factor(dat$Age)
dat$Age <- relevel(dat$Age, ref = "4")
dat$BlockTotal <- as.factor(dat$BlockTotal)
dat$SoundType_re <- ifelse(dat$SoundType == "DEV", 1,0)
dat$SoundType_re <- as.factor(dat$SoundType_re)
dat$BlockOrder <- factor(dat$BlockOrder)
# m <- RT ~ 1 + SoundType_re + Age + BlockTotal + SoundType_re:BlockTotal + SoundType_re:Age +
#   Age:BlockTotal + SoundType_re:Age:BlockTotal + (1 + SoundType_re | VP)

m <- RT ~ 1 + SoundType_re*Age*BlockOrder + (1 + SoundType_re | VP)

fit<-lmer(m, data=dat, REML = T)  
summary<- summary(fit)  

summary$varcor



fit_2<- nlme::lme(fixed = RT ~ SoundType_re*Age*BlockOrder ,
                   random =  ~ SoundType_re | VP, data=dat, 
                  weights =  nlme::varIdent(form = ~1|factor(Age)),
                  control = nlme::nlmeControl(maxIter = 2000, msMaxIter = 2000, pnlsMaxIter = 2000,
                                        niterEM = 1000, opt = "nlminb"))
summary(fit_2)

??nlme::lme



fit_3<-robust_mixed(lmer(RT ~ 1 + SoundType_re*Age*BlockOrder + 
                           (1 + SoundType_re | VP), data = dat),type = "CR2")
fit_3$results


fit_4 <-brm(bf(RT ~ 1 + SoundType_re*Age*BlockOrder + 
              (1 + SoundType_re | VP),
              sigma ~ 1 + SoundType_re*Age*BlockOrder + 
                (1 + SoundType_re | VP)),
    data = dat,
    iter = 7000,
    cores =  2,
    save_pars = save_pars(all = TRUE))

print(fit_4,digits=3)

fit_5 <-brm(bf(RT ~ 1 + SoundType_re*Age*BlockOrder + 
                 (1 + SoundType_re | VP),
               sigma ~ 1 + SoundType_re*Age*BlockOrder),
            data = dat,
            iter = 7000,
            cores =  2,
            save_pars = save_pars(all = TRUE))

print(fit_5,digits=3)


fit_4 <-brm(bf(RT ~ 1 + SoundType_re*Age*BlockOrder + 
                 (1 + SoundType_re | VP),
               sigma ~ 1 + Age,
            data = dat,
            iter = 10000,
            cores =  2,
            save_pars = save_pars(all = TRUE))

WAIC(fit_4,fit_5)


#### ####
dat4 <- subset(dat, Age == 4)
dat5 <- subset(dat, Age == 5)
dat6 <- subset(dat, Age == 6)
dat7 <- subset(dat, Age == 7)
dat8 <- subset(dat, Age == 8)
dat9 <- subset(dat, Age == 9)
dat10 <- subset(dat, Age == 10)
dat_adult <- subset(dat, Age == "Adult")

m2 <- RT ~ 1 + SoundType_re + BlockOrder + SoundType_re:BlockOrder + (1 + SoundType_re | VP)

fit4 <- lmer(m2, data = dat4 , REML = T)
fit5 <- lmer(m2, data = dat5 , REML = T)
fit6 <- lmer(m2, data = dat6 , REML = T)
fit7 <- lmer(m2, data = dat7 , REML = T)
fit8 <- lmer(m2, data = dat8 , REML = T)
fit9 <- lmer(m2, data = dat9 , REML = T)
fit10 <- lmer(m2, data = dat10 , REML = T)
fit_adult <- lmer(m2, data = dat_adult , REML = T)

summary(fit4)  
summary(fit5) 
summary(fit6) 
summary(fit7) 
summary(fit8) 
summary(fit9) 
summary(fit10) 
summary(fit_adult) 
#### ####







