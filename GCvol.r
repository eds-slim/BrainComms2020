stats.GC.vol<-d.conn %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  nlme(conn.meas.value~deriv_mod(tp,a,b,del)
       , fixed=list(a+b~1,del~logvol*relpos)
       , random=a~1|numID
       , data=.
       , start=c(1,1,1,1,1,1)
       , na.action = na.omit) %>% 
  tidy()

tbl.GC.vol<-stats.GC.vol %>% 
  dplyr::select(-c(effect,group)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GC.vol %>% 
  write.table(file = './../Word/Floats/Tables/07-GC_logvol.tab', sep = ",", quote = FALSE, row.names = F)

stats.GC.vol.relpos <- d.conn %>% 
  merge(data.clinical_long) %>% 
  dplyr::filter(clin.meas.name=='NIHSS') %>%
  group_by(conn.meas.name,relpos) %>%
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'conn.meas.value', preds=list(a='1', b='1',del='logvol'))$result)) %>% 
  mutate(tidy = map(mdl, tidy), AIC = map(mdl,AIC)) %>% 
  unnest(tidy)


tbl.GC.vol.relpos <- stats.GC.vol.relpos %>% 
  select(-c(effect,group)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GC.vol.relpos %>% 
  write.table(file = './../Word/Floats/Tables/07-GC_logvol_by_relpos.tab', sep = ",", quote = FALSE, row.names = F)




dd<-d.conn %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS')

labelInfo <-
  split(dd, list(dd$relpos, as.factor(dd$tp))) %>%
  lapply(function(t){
    data.frame(
      predAtMax = lm(conn.meas.value ~ logvol, data = t, na.action = na.omit) %>%
        predict(newdata = data.frame(logvol = max(t$logvol, na.rm = T)))
      , max = max(t$logvol, na.rm = T),
      relpos = t$relpos[1]
    )}) %>%
  bind_rows

labelInfo$label = levels(factor(interaction(dd$relpos, as.factor(dd$tp))))
labelInfo$label2<-unlist(lapply(as.character(labelInfo$label), reformatlabel))
labelInfo

p.GC.vol<-ggplot(data=dd, aes(x=log(vol), y=conn.meas.value, shape=as.factor(tp), color=relpos))+
  #geom_point(alpha=.16, size=1)+
  geom_smooth(aes(group=interaction(as.factor(tp),relpos),color=relpos
                  , linetype=as.factor(tp)), size=.5, alpha=.5, method="lm"
              , formula = 'y~x', se=F)+
  ylab("Median connectivity")+xlab("Lesion volume [ml] (log)")+
  geom_label_repel(data = labelInfo
                   , aes(x = max
                         , y = predAtMax
                         , label = label2
                         , color = relpos
                   )
                   , nudge_x = 3, inherit.aes = F, show.legend = F
                   , size=8/ggplot2:::.pt, segment.size = .1, label.size = .1)+
  coord_cartesian(xlim = c(min(dd$logvol, na.rm = T),10)-1,ylim=c(.425,.5))+
  scale_x_continuous(breaks=c(0,2,4,6))+
  guides(color=guide_legend(title="Hemisphere",order=1), linetype=F, shape=guide_legend(title='Time post stroke', order=0), text=F)+
  guides(shape=F, linetype=F, color=F)+
  theme_pubr()+
  theme(panel.grid.major.x = element_line(colour = 'lightgrey', size = .2, linetype = 2)
        ,axis.text = element_text(size = 12,hjust = c(1,.25,0,.5)))
#theme(legend.position = c(.15,.2), legend.direction = 'vertical')+
  #thm
p.GC.vol
print(p.GC.vol)

ggsave('./../Word/Floats/Figures/08-GC_logvol.tiff', device = 'tiff', width = 8.5, height = 8.5, units = 'cm', dpi = 600)



