medFit <- glm(Precentral ~ treatment + NIHSS + age, data = ddd, family = gaussian)

expData <- neWeight(medFit)


neMod1 <- neModel(goodOutcome ~ factor(treatment0) + factor(treatment1) + NIHSS, family = binomial, expData = expData)

summary(neMod1)


impData <- neImpute(goodOutcome ~ treatment + 
                      Precentral * Postcentral + Insula + Parsopercularis + Parstriangularis + 
                      Superiortemporal + Transversetemporal + Supramarginal + Caudalmiddlefrontal +
                      age + NIHSS
                    , family = binomial("logit"), nMed = 9, data = ddd)
neMod6 <- neModel(goodOutcome ~ treatment0 + treatment1 + age + NIHSS, family = binomial("logit"), expData = impData, se = "robust")
summary(neMod6)
confint(neMod6)


plot(neMod6)

d.tmp <- ddd %>% filter(vol0>0 & !is.nan(goodOutcome))
d.tmp$llv0 <- log(d.tmp$vol0)
d.tmp %>% nrow()
impData <- neImpute(goodOutcome ~ NIHSS + 
                      Precentral * Postcentral + Insula + Parsopercularis + Parstriangularis + 
                      Superiortemporal + Transversetemporal + Supramarginal + Caudalmiddlefrontal
                    , family = binomial(), nMed = 9, data = d.tmp)
neMod6 <- neModel(goodOutcome ~ NIHSS0 + NIHSS1, family = binomial("logit"), expData = impData, se = "robust")
summary(neMod6)
confint(neMod6)
plot(neMod6)
