
## ---------------------------------------------------------
## Gammatype variants survival analysis
## ---------------------------------------------------------


## Libs

library(data.table)
library(tidyverse)
library(ranger)
library(survival)
library(ggpubr)


## Data

# read data
gt2 <- fread('./data/GT_match_FU', data.table=F, header=T, na.strings=c("NA", "na")) 
# combine follow-up columns
gt2 <- mutate(gt2, newFU=sapply(1:nrow(gt2), function(x) ifelse(gt2$LastFU[x]=='', gt2$LastFU_2[x], gt2$LastFU[x])))
# date format
gt2$newFU <- gt2$newFU %>% as.Date(., format="%d.%m.%Y")
gt2$Txdate <- gt2$Txdate %>% as.Date(., format="%d.%m.%Y")
# remove missing
gt2 <- gt2[complete.cases(gt2), ]
# create survival time variable
gt2 <- mutate(gt2, time=as.numeric(newFU-Txdate))
# create survival status variable
gt2$status <- c(0,1)[gt2$status %>% factor]
# arrange cols
gt2 <- gt2[, c(1, 50:53, 43:47, 2:42)]
# create binary variables
gt2 <- mutate(gt2, cGvHD2=c('no','yes','yes')[cGvHD %>% factor] %>% factor)
gt2 <- mutate(gt2, ptage2=ifelse(ptage>median(ptage), '>53', '<53') %>% factor)
gt2 <- mutate(gt2, dnage2=ifelse(dnage>median(dnage), '>34', '<34') %>% factor)
gt2 <- mutate(gt2, GT_MATCH2=ifelse(GT_MATCH>24, 'full', 'mismatch') %>% factor)


## Analyses

# explore importance with ranger survival
gt2.ranger <- ranger(Surv(time, status)~., data=gt2[, c(5, 9, 11:ncol(gt2))], importance="permutation")
gt2.ranger.vi <- round(gt2.ranger$variable.importance, 4) %>% sort(., decreasing=T)

# select and arrange variables for CV
gt2.tmp <- gt2[, c(5, 9, 11:ncol(gt2))]
gt2.tmp <- gt2.tmp[, -c(44:47)]

# CV ranger survival
cv.rng <- lapply(1:nrow(gt2.tmp), function(i) {
  tmp.rng <- ranger(Surv(time, status)~., data=gt2.tmp[-i, ], mtry=41, write.forest=T, num.tree=1000, importance="permutation")
  tmp <- predict(tmp.rng, gt2.tmp[i, -c(1,2)])
  cbind(tmp$unique.death.times, tmp$survival, tmp.rng$variable.importance)
})
cv.rng.vi <- sapply(cv.rng, function(x) x[1:41, 3]) %>% t
colnames(cv.rng.vi) <- colnames(gt2.tmp)[3:43]
cv.rng.vi.top <- cv.rng.vi[, cv.rng.vi %>% colMeans %>% order(., decreasing=T)][, 1:7]

# fit Kaplan-Meier survival models for 4 variables
gt2.km1  <- survfit(Surv(time, status) ~ cGvHD2, data=gt2)
gt2.km2  <- survfit(Surv(time, status) ~ dnage2, data=gt2)
gt2.km3  <- survfit(Surv(time, status) ~ ptage2, data=gt2)
gt2.km4  <- survfit(Surv(time, status) ~ GT_MATCH2, data=gt2)

# fit Cox model
gt2.cox <- coxph(Surv(time, status) ~dnage+cGvHD2+ptage+GT_MATCH2, data=gt2)


## Plots

# variable importances; original and log
p11 <- ggplot(cv.rng.vi.top %>% data.frame %>% gather, aes(x=key, y=value %>% log)) + 
  geom_boxplot() +
  scale_x_discrete(limits=cv.rng.vi.top %>% colnames) +
  xlab('Variable') +
  ylab('log(Importance)') +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
p12 <- ggplot(cv.rng.vi.top %>% data.frame %>% gather, aes(x=key, y=value)) + 
  geom_boxplot() +
  scale_x_discrete(limits=cv.rng.vi.top %>% colnames) +
  ylab('Importance') +
  xlab('Variable') +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

pdf(height=3.5, width=7, file='./results/Survival_Importances.pdf')
ggpubr::ggarrange(p12, p11, nrow=1, ncol=2, align='h')
dev.off()

# plot KM's
p1 <- autoplot(gt2.km1) + guides(color=guide_legend(title='cGvHD'), fill=guide_legend(title='cGvHD'))
p2 <- autoplot(gt2.km2) + guides(color=guide_legend(title='Donor age'), fill=guide_legend(title='Donor age'))              
p3 <- autoplot(gt2.km3) + guides(color=guide_legend(title='Patient age'), fill=guide_legend(title='Patient age'))               
p4 <- autoplot(gt2.km4) + guides(color=guide_legend(title='GT match'), fill=guide_legend(title='GT match'))               

# plot cox summary
p5 <- ggtexttable(summary(gt2.cox)$coefficients %>% round(., 4) %>% data.frame, theme=ttheme("minimal"))

pdf(height=7, width=8, file='./results/cox_survivals.pdf')
ggarrange(ggarrange(p1, p2, p3, p4, nrow=2, ncol=2, align='v'),
          p5, ncol=1, nrow=2, heights=c(3, 1))
dev.off()



