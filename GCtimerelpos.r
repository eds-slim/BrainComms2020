UKEcols<-c(rgb(170,156,143, maxColorValue = 255),rgb(0,73,146, maxColorValue = 255),rgb(178,34,41, maxColorValue = 255))
UKE <- define_palette(
  swatch = UKEcols,
  gradient = c(lower = UKEcols[1], upper = UKEcols[3])
)

ggthemr(UKE)
ggthemr('light')
p.GC.timecourse<-d.conn %>%
  ggplot(aes(x=tp, y=conn.meas.value, group=interaction(numID,relpos),  color=relpos, fill=relpos))+
  #geom_line(size=.2, alpha=.2)+
  #geom_jitter(size=.5,alpha=.1, width = .2, shape=21)+
  stat_summary(aes(group = relpos), geom = "point", fun.y = mean, shape = 21, size = .75, color='grey')+
  stat_summary(aes(group = relpos, width=.75), size=.5, fun.data = mean_se, geom = "errorbar")+
  geom_smooth(aes(group=relpos), method = "nls"
              , formula = y ~ a+del*(1-exp(-b*x)), se = F, size=.5
              , method.args = list(control=nls.control(warnOnly = T, maxiter=1e5, minFactor = 1e-20), start = list(a = 0.5, b=0.5, del=-0.1)))+
  labs(y="Median connectivity")+
  sc_x_tp+
  guides(color=F, fill=F)+
  theme_pubr()+
  theme_update(panel.grid.major.x = element_line(colour = 'lightgrey', size = .2, linetype = 2)
        ,axis.text = element_text(size = 12,hjust = c(1,.25,0,.5)))
p.GC.timecourse

ggsave('./../Word/Floats/Figures/02-GC.tiff', plot = p.GC.timecourse, device = 'tiff', width = 8.5, height = 8.5, units = 'cm', dpi = 600)

stats.GC.exp <- d.conn %>%
  group_by(conn.meas.name) %>%
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'conn.meas.value', preds=list(a='1', b='1',del='relpos'))$result)) %>% 
  mutate(tidy = map(mdl, tidy), AIC = map(mdl,AIC))


AIC.GC.exp<-stats.GC.exp %>% 
  unnest(AIC,.drop = T)

tbl.GC.exp<-stats.GC.exp %>% 
  unnest(tidy,.drop = T) %>% 
  dplyr::select(-c(effect,group)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))

tbl.GC.exp %>% 
  write.table(file = './../Word/Floats/Tables/04-GC_timecourse.tab', sep = ",", quote = FALSE, row.names = F)


stats.GC.exp.relpos <- d.conn %>%
  group_by(conn.meas.name, relpos) %>%
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'conn.meas.value', preds=list(a='1', b='1',del='1'))$result)) %>% 
  mutate(tidy = map(mdl, tidy), AIC = map(mdl,AIC))

AIC.GC.exp.relpos<-stats.GC.exp.relpos %>% 
  unnest(AIC,.drop = T)

tbl.GC.exp.relpos<-stats.GC.exp.relpos %>% 
  unnest(tidy,.drop = T) %>% 
  dplyr::select(-c(effect,group)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))

tbl.GC.exp.relpos %>% 
  write.table(file = './../Word/Floats/Tables/05-GC_timecourse_by_relpos.tab', sep = ",", quote = FALSE, row.names = F)


stats.GC.lin <- d.conn %>%
  group_by(conn.meas.name) %>%
  nest() %>%
  mutate(mdl=map(data,~lme(conn.meas.value~tp*relpos, random=~1|numID, data=., na.action = na.omit, method = 'ML'))) %>% 
  mutate(tidy = map(mdl, tidy), AIC = map(mdl,AIC))

AIC.GC.lin<-stats.GC.lin %>% 
  unnest(AIC,.drop = T)

stats.GC.lin.by.relpos <- d.conn %>%
  group_by(conn.meas.name, relpos) %>%
  nest() %>%
  mutate(mdl=map(data,~lme(conn.meas.value~tp, random=~1|numID, data=., na.action = na.omit, method = 'ML'))) %>% 
  mutate(tidy = map(mdl, tidy), AIC = map(mdl,AIC))

AIC.GC.lin.relpos<-stats.GC.lin.by.relpos %>% 
  unnest(AIC,.drop = T)