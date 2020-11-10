p.GGP.timecourse<-d.GGP %>% 
  filter(threshold == 1) %>% 
  ggplot(aes(x = tp-.3+as.numeric(relpos)/5, y = GGP
             , group = interaction(numID,relpos), color=relpos, shape=side
             , fill = relpos
             )
         )+
  #geom_line(size=.1, alpha=.2)+
  #geom_jitter(size=.5,alpha=.1, width = .2, shape=21)+
  sc_x_tp + 
  scale_y_continuous('')+
  geom_smooth(aes(group=relpos), method = "nls"
              , formula = y ~ a*exp(-b*x)+c, se = F
              , method.args = list(control=nls.control(warnOnly = F, maxiter=1e5, minFactor = 1e-10), start = list(a = 0.1, b=0.3, c=0.5))
              , size = .5)+
  stat_summary(aes(group = relpos, width=.75), size=.25, fun.data = mean_se, geom = "errorbar", linetype=1)+
  stat_summary(aes(group = relpos), geom = "point", fun.y = mean, shape = 21, size = .5, color='black', stroke=.2)+
  facet_wrap(. ~ lab, scales="free", labeller = labeller(lab = GGP_names)
             , ncol = length(unique(d.GGP$lab))#, strip.position = "left"
             )+
  guides(color=F, fill=F, shape=F, linetype=F)+
  #thm+
  theme(strip.background = element_blank(), strip.placement = "inside"
        , text = element_text(size = 8)
        , axis.text = element_text(size = 4)
        , axis.text.x = element_text(hjust = c(.75,.25,.5,.5))
        , axis.line = element_line(size = .25, colour = "black")
        , axis.ticks = element_line(size = .25, colour = 'black'))


p.GGP.timecourse
ggsave('./../Word/Floats/Figures/2-GGP.tiff', plot = p.GGP.timecourse, device = 'tiff'
       , width = 9, height = 4.5, units = 'cm', dpi = 1200)


stats.GGP<-d.GGP %>%
  group_by(lab, threshold) %>%
  arrange(lab, threshold) %>% 
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'GGP', preds=list(a='1', b='1',del='relpos'))$result)) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


AIC.GGP<-stats.GGP %>% 
  unnest(AIC,.drop = T) %>% 
  filter(threshold == 1)

tbl.GGP<-stats.GGP %>% 
  select(tidy) %>% unnest(tidy,.drop = T) %>%
  filter(threshold == 1) %>% 
  dplyr::select(-c(effect,group,threshold)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GGP %>%   
  write.table(file = './../Word/Floats/Tables/05-GGP_timecourse.tab', sep = ",", quote = FALSE, row.names = F)



stats.GGP.relpos<-d.GGP %>%
  group_by(lab, threshold, relpos) %>%
  arrange(lab, threshold) %>% 
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'GGP', preds=list(a='1', b='1',del='1'))$result)) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


AIC.GGP.relpos<-stats.GGP.relpos %>% 
  select(AIC) %>% unnest(AIC,.drop = T) %>% 
  filter(threshold == 1)

tbl.GGP.relpos<-stats.GGP.relpos %>% 
  select(tidy) %>% unnest(tidy,.drop = T) %>%
  filter(threshold == 1) %>% 
  dplyr::select(-c(effect,group,threshold)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GGP.relpos %>%   
  write.table(file = './../Word/Floats/Tables/05-GGP_timecourse_by_relpos.tab', sep = ",", quote = FALSE, row.names = F)


################################

stats.GGP.lin<-d.GGP %>%
  group_by(lab, threshold) %>%
  arrange(lab, threshold) %>% 
  nest() %>%
  mutate(mdl=map(data,~lme(GGP ~ tp*relpos, random = ~1|numID, data = ., na.action = na.exclude))) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


AIC.GGP.lin <- stats.GGP.lin %>% 
  select(AIC) %>% unnest(AIC,.drop = T) %>% 
  filter(threshold == 1)


stats.GGP.relpos.lin<-d.GGP %>%
  group_by(lab, threshold, relpos) %>%
  arrange(lab, threshold) %>% 
  nest() %>%
  mutate(mdl=map(data,~lme(GGP ~ tp, random = ~1|numID, data = ., na.action = na.exclude))) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


AIC.GGP.relpos.lin <- stats.GGP.relpos.lin %>% 
  select(AIC) %>% unnest(AIC,.drop = T) %>% 
  filter(threshold == 1)



#########################################

d <- d.GGP %>%
  merge(data.clinical_long) %>%  filter(clin.meas.name == 'NIHSS') %>% 
  ungroup() %>% 
  filter(threshold == 1 & relpos != '!ipsi') %>% dplyr::select(-GGP.diff) %>% 
  mutate(lab = as_labeller(GGP_names)(lab)) %>% 
  spread(lab,GGP)
idx.complete <- d %>% dplyr::select(c(numID,tp,Efficiency, q50)) %>% complete.cases()

fm <- nlme.exp(d, response = 'q50', preds=list(a='1', b='1',del='relpos'))
d$.fitted <- NA
d$.fitted[idx.complete] <- fitted(fm)

d %>% subset(tp>=0) %>% 
  ggplot(aes(x=.fitted, y = Efficiency, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)

lm(Efficiency ~ .fitted + relpos*exp(-.5*tp), data = d) %>% summary()
