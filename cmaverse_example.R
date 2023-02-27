## Load libraries
library(data.table)
library(ggplot2)
library(CMAverse)
library(paths)
library(gbm)

###############################################################################################
## 1) SET UP DATA GENERATING PROCESS 
###############################################################################################

## Edited version of cmaverse example from https://bs1125.github.io/CMAverse/articles/post_exposure_confounding.html
## Tweak these parameters of the true data generating process and see what happens!
set.seed(1)
expit <- function(x) exp(x)/(1+exp(x))
n <- 10000
C1 <- rnorm(n, mean = 1, sd = 0.1)
C2 <- rnorm(n, mean = 0.6, sd = 0.1)
A <- rbinom(n, 1, expit(0.2 + 0.5*C1 + 0.1*C2))
L <- rnorm(n, mean = 1 + 5*A - C1 - 0.5*C2, sd = 1)
M <- rnorm(n, mean = (1 + 5*A + 10*L + 1.5*C1 + 0.8*C2), sd=1)
Y <- rnorm(n, mean = (10 + 0.4*A + 1.2*M + 0.5*A*M + 10*L + 0.3*C1 + 0.6*C2), sd=20)
data <- data.frame(A, M, Y, C1, C2, L)

###############################################################################################
## 2) ONE MEDIATOR
###############################################################################################

## Estimate CDE using Baron-Kenny mediation
reg <- lm(Y ~ C1 + C2 + A*M)
out1 <- data.table(effect='cde (Baron-Kenny)',mean=coef(reg)[['A']],lower=confint(reg)['A',1],upper=confint(reg)['A',2])

## Estimate all mediation estimands using g-formula
res_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                      mediator = "M", basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                      mreg = list("linear"), yreg = "linear", postcreg = list("linear"),
                      astar = 0, a = 1, mval = list(0), 
                      estimation = "imputation", inference = "bootstrap", nboot = 100)
out2 <- summary(res_gformula)
out2 <- data.table(effect=names(out2$effect.pe),mean=out2$effect.pe,lower=out2$effect.ci.low,upper=out2$effect.ci.high)

## Plot
out <- rbind(out1,out2,fill=T)
out <- out[effect %in% c('te','cde (Baron-Kenny)','cde','rintref','rintmed','rpnie')]
out[, effect := factor(effect,levels=c('te','cde (Baron-Kenny)','cde','rintref','rintmed','rpnie'))]
ggplot(data=out,
       aes(x=effect,y=mean,ymin=lower,ymax=upper)) + 
  geom_linerange() +
  geom_point(size=3) + 
  theme_minimal()

###############################################################################################
## 3) MULTIPLE MEDIATORS
###############################################################################################

## Estimate all mediation estimands using g-formula with two mediators (L, M).
res_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                      mediator = c("L","M"), basec = c("C1", "C2"), EMint = TRUE,
                      mreg = list("linear","linear"), yreg = "linear",
                      astar = 0, a = 1, mval = list(0,0), 
                      estimation = "imputation", inference = "bootstrap", nboot = 100)
out3 <- summary(res_gformula)
out3 <- data.table(effect=names(out3$effect.pe),mean=out3$effect.pe,lower=out3$effect.ci.low,upper=out3$effect.ci.high)
out3 <- out3[effect %in% c('te','cde','intref','intmed','pnie')]
out3[, effect := factor(effect,levels=c('te','cde','intref','intmed','pnie'))]
ggplot(data=out3,
       aes(x=effect,y=mean,ymin=lower,ymax=upper)) + 
  geom_linerange() +
  geom_point(size=3) + 
  theme_minimal()

## Estimate path-specific effects of L and M using paths package. Use gradient boosting machines (could also use Bayesian additive regression trees; BART).
gbm_m1 <- gbm(Y ~ A + C1 + C2, data = data, distribution = "gaussian", interaction.depth = 3)
gbm_m2 <- gbm(Y ~ A + L + C1 + C2, data = data, distribution = "gaussian", interaction.depth = 3)
gbm_m3 <- gbm(Y ~ A + L + M + C1 + C2, data = data, distribution = "gaussian", interaction.depth = 3)
gbm_ymodels <- list(gbm_m1,gbm_m2,gbm_m3)
gbm_ps <- gbm(A ~ C1 + C2, data = data, distribution = "gaussian", interaction.depth = 3)
paths_gbm <- paths(a = "A", y = "Y", m = list('L','M'), gbm_ymodels, ps_model = gbm_ps, data = data, nboot = 100)
out4 <- as.data.table(paths_gbm$hybrid)
out4 <- out4[decomposition=='Type I']
out4[, estimand := factor(estimand,levels=c('total','direct','via M1','via M2'))]
ggplot(data=out4,
       aes(x=estimand,y=estimate,ymin=lower,ymax=upper)) + 
  geom_linerange() +
  geom_point(size=3) + 
  theme_minimal()
