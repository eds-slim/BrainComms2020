
stats.GGP.vol <- d.GGP %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  group_by(lab, threshold) %>% 
  arrange(lab) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(GGP~deriv_mod(tp,a,b,del)
                                    , fixed=list(a+b~1,del~logvol*relpos)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))



tbl.GGP.vol <- stats.GGP.vol %>% 
  select(tidy) %>% unnest(tidy,.drop = T) %>%
  filter(threshold == 1) %>% 
  dplyr::select(-c(effect,group,threshold)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GGP.vol %>% 
  write.table(file = './../Word/Floats/Tables/08-GGP_logvol.tab', sep = ",", quote = FALSE, row.names = F)

stats.GGP.vol.relpos <- d.GGP %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  group_by(lab, threshold,relpos) %>% 
  arrange(lab) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'GGP', preds=list(a='1', b='1',del='logvol'))$result)) %>% 
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


tbl.GGP.vol.relpos <- stats.GGP.vol.relpos %>% 
  select(tidy) %>% unnest(tidy,.drop = T) %>% 
  filter(threshold == 1) %>% 
  dplyr::select(-c(effect,group,threshold,df,statistic)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GGP.vol.relpos %>% 
  write.table(file = './../Word/Floats/Tables/08-GGP_logvol_by_relpos.tab', sep = ",", quote = FALSE, row.names = F)  




GGPvol_plt_fcn<-function(d,lab){
  
  labelInfo <-
    split(d, list(d$relpos, as.factor(d$tp))) %>%
    lapply(function(t){
      data.frame(
        predAtMax = lm(GGP ~ logvol, data = t) %>%
          predict(newdata = data.frame(logvol = max(t$logvol, na.rm = T)))
        , max = max(t$logvol,na.rm = T),
        relpos = t$relpos[1]
      )}) %>%
    bind_rows
  
  labelInfo$label = levels(factor(interaction(d$relpos, as.factor(d$tp))))
  labelInfo$label2<-unlist(lapply(as.character(labelInfo$label), reformatlabel))
  labelInfo
  
  scaleFUN <- function(x) sprintf("%.3f", x)
  
  p.GGP.vol<-ggplot(data=d, aes(x=log(vol), y=GGP, shape=as.factor(tp), color=relpos))+
    geom_point(alpha=.32, size=.5)+
    geom_smooth(aes(group=interaction(relpos,tp), color=relpos, linetype=as.factor(tp)), size=0.25, method="lm", se=F)+
    guides(shape=F, linetype=guide_legend(title="Time point"), color=guide_legend(title="Hemisphere", labels=c("Stroke","Intact")))+
    theme(legend.position = "top",legend.spacing.y = unit(-5, "cm"), legend.margin = margin(-0.5,0,0,0, unit="cm"))+
    theme(legend.position = c(0.15,.2), legend.direction = "vertical")+
    #geom_label_repel(data = labelInfo
    #                 , aes(x = max
    #                       , y = predAtMax
    #                       , label = label2
    #                       , color = relpos
    #                 )
    #                 , nudge_x = 5, inherit.aes = F, show.legend = F, max.iter = 1e6
    #                 , size=.5, segment.size = .1, label.size = .1
    #                 , label.padding = .1, point.padding = 0, force = 0.1, min.segment.length = 0)+
    #coord_cartesian(xlim = c(min(d$logvol, na.rm = T)-1,9))+
    scale_x_continuous(if_else(lab==1,"Lesion volume (log)",''), breaks=c(0,2,4,6))+
    scale_y_continuous('', breaks=pretty_breaks(3), labels=scaleFUN) +
    ggtitle(as_labeller(GGP_names)(lab)) +
    scale_shape_manual(values=c(15,16,17,18),labels=c('3-5 d','1 m','3 m','12 m'))+
    scale_linetype_manual(values = c('solid','31','3111','11'))+
    guides(color=guide_legend(title="Hemisphere",order=1), linetype=F, shape=guide_legend(title='Time post stroke', order=0), text=F)+
    guides(shape=F, linetype=F, color=F)+
    theme(strip.background = element_blank(), strip.placement = "inside"
          #, text = element_text(size = 8)
          , axis.text = element_text(size = 6)
          , axis.title = element_text(size = 8)
          , axis.text.x = element_text(hjust = c(.75,.25,.5,.5))
          , axis.line = element_line(size = .25, colour = "black")
          , axis.ticks = element_line(size = .25, colour = 'black')
          , title = element_text(size=6, hjust = .5)
          , plot.title = element_text(size=6, hjust = .5)
          #, plot.margin = unit(c(0,0,0,-.2), "cm")
          )
    #theme_pubr()+
    #theme(panel.grid.major.x = element_line(colour = 'lightgrey', size = .2, linetype = 2)
    #      , axis.text.x = element_text(hjust = c(.75,.25,.5,.5)))
    #thm
}


p.GGP.vol.list<-d.GGP %>% 
  filter(threshold==1) %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  group_by(lab) %>% 
  arrange(lab) %>% 
  nest() %>% 
  mutate(plt=pmap(list(data,lab),~GGPvol_plt_fcn(.x,.y)))

p.GGP.vol<-do.call(arrangeGrob,c(p.GGP.vol.list$plt
                                 , ncol=length(unique((d.GGP$lab)))
                                 , as.table=F
                                 , padding=0)
                   )
p<-p.GGP.vol %>% as_ggplot()

scaleFUN <- function(x) sprintf("%.3f", x)

p.GGP.vol <- d.GGP %>% 
  filter(threshold==1) %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  ggplot(aes(x=log(vol), y=GGP, shape=as.factor(tp), color=relpos))+
  geom_point(alpha=.5, size=.1)+
  geom_smooth(aes(group=interaction(relpos,tp), color=relpos, linetype=as.factor(tp)), size=0.25, method="lm", se=F)+
  #theme(legend.position = "top",legend.spacing.y = unit(-5, "cm"), legend.margin = margin(-0.5,0,0,0, unit="cm"))+
  #theme(legend.position = c(0.15,.2), legend.direction = "vertical")+
  scale_x_continuous("Lesion volume (log)", breaks=c(0,2,4,6))+
  scale_y_continuous('', breaks=pretty_breaks(3), labels=scaleFUN) +
  scale_shape_manual(values=c(15,16,17,18),labels=c('3-5 d','1 m','3 m','12 m'))+
  scale_linetype_manual(values = c('solid','31','3111','11')
                        , labels = c('3-5d', '1m', '3m', '12m'))+
  scale_colour_manual(values = swatch()[c(2,3)], labels = c('Contralesional', 'Ipsilesional'))+
  facet_wrap(. ~ lab, scales="free", labeller = labeller(lab = GGP_names)
             , ncol = length(unique(d.GGP$lab))#, strip.position = "left"
  )+
  guides(shape=F
         , linetype = guide_legend(title = 'Time point', override.aes = list(color='black'), direction = 'vertical', nrow = 2)
         , colour = guide_legend(title = 'Hemisphere condition', override.aes = list(shape=NA)))+
  theme(strip.background = element_blank(), strip.placement = "inside"
        , text = element_text(size = 8)
        , axis.text = element_text(size = 4)
        , axis.line = element_line(size = .25, colour = "black")
        , axis.ticks = element_line(size = .25, colour = 'black')
        , panel.grid.major.x = element_blank()
        )

p.GGP.vol.leg <- p.GGP.vol +
                  theme(legend.position = c(.51, -.26)
                  , legend.spacing = unit(3.4, 'cm')
                  #, legend.margin = unit(1.7, 'cm')
                  , legend.key.height = unit(0, 'cm')
                  , legend.direction = 'vertical'
                  , legend.box = 'horizontal'
                  , legend.text = element_text(size = 4)
                  , legend.title = element_text(size = 6)
                  ,legend.background = element_rect(fill = "white", colour = 'black', size = .1)
                  , legend.margin = margin(t = 1, r = 1, b = 1, l = 1)
                  )

p.GGP.vol.leg

ggsave('./../Word/Floats/Figures/3-GGPvol.tiff', plot = p.GGP.vol.leg, device = 'tiff', width = 9, height = 4.5, units = 'cm', dpi = 1200)



######################
##### sensi ana ######
######################ø


p.GGP.vol.sensiana.g1<-stats.GGP.vol %>% 
  select(tidy) %>% unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed'  & lab != 0) %>%
  mutate(term = dplyr::recode(term, `del.(Intercept)`='Δ', `del.logvol`='Δ~vol~', `del.logvol:relposipsi`='Δ~vol:ipsi~', `del.relposipsi`='Δ~ipsi~', )) %>% 
  mutate(term = forcats::fct_relevel(term, 'a', 'b', 'Δ', 'Δ~ipsi~', 'Δ~vol~')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = 'lightgrey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size = .5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_wrap(lab~term
             , scale = "free_x"
             , ncol = 6
             , labeller = labeller(lab=GGP_names)) +
  scale_x_continuous('Network density', breaks = seq(.1,1,.1))+
  scale_y_continuous('Estimate', breaks = pretty_breaks(3), labels = format.fcn)+
  scale_color_aaas()+
  scale_linetype_manual(values = c('solid','11'))+
  theme(
    strip.background = element_blank()
    , strip.text.x = element_blank()
    , panel.spacing.x = unit(4, "mm")
    , axis.text = element_text(size = 6)
    , panel.grid.major.x = element_blank()
    , axis.line = element_line(size = .25, colour = "black")
    , axis.ticks = element_line(size = .25, colour = 'black')
  )+
  guides(color=F,linetype=F)

p.GGP.vol.sensiana.g2<-stats.GGP.vol %>%
  select(tidy) %>% unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed'  & lab != 0) %>%
  mutate(term = dplyr::recode(term, `del.(Intercept)`='Δ', `del.logvol`='Δ~vol~', `del.logvol:relposipsi`='Δ~vol:ipsi~', `del.relposipsi`='Δ~ipsi~', )) %>% 
  mutate(term = forcats::fct_relevel(term, 'a', 'b', 'Δ', 'Δ~ipsi~', 'Δ~vol~')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = 'lightgrey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size = .5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_grid(lab~term
             , scales = "free_x"
             #, ncol = 6
             , labeller = labeller(lab=GGP_names)) +
  scale_x_continuous('Network density', breaks = seq(.1,1,.1))+
  scale_y_continuous('Estimate', breaks = pretty_breaks(3))+
  scale_color_aaas()+
  scale_linetype_manual(values = c('solid','11'))+
  theme(
    panel.spacing.x = unit(4, "mm")
  )+
  guides(color=F,linetype=F)

p.GGP.vol.sensiana <- change.facets.fcn(p.GGP.vol.sensiana.g1, p.GGP.vol.sensiana.g2)
ggsave('./../Word/Floats/Figures/S5-GGP-logvol-sensiana.tiff', plot = p.GGP.vol.sensiana, device = 'tiff'
       , width = 18.5, height = 9.25, units = 'cm', dpi = 1200)



p.GGP.vol.relpos.sensiana.g1 <- stats.GGP.vol.relpos %>% 
  select(tidy) %>% unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed' & !(threshold %in% c(2) & lab==30) & lab != 0) %>%
  ungroup() %>% 
  mutate(lab = as_labeller(GGP_names)(lab) %>% as.character()) %>% 
  mutate(term = dplyr::recode(term, `del.(Intercept)`='Δ', `del.logvol`='Δ~vol~')) %>% 
  mutate(term = forcats::fct_relevel(term, 'a', 'b', 'Δ', 'Δ~vol~')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = 'lightgrey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size =.5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_wrap(lab*relpos~term
             , scale = "free_x"
             , ncol = 8) +
  scale_x_continuous('Network density', breaks = seq(.1,1,.1))+
  scale_y_continuous('Estimate', breaks = pretty_breaks(3), labels = format.fcn)+
  scale_color_aaas()+
  scale_linetype_manual(values = c('solid','11'))+
  theme(
    strip.background = element_blank()
    , strip.text.x = element_blank()
    , panel.spacing.x = unit(4, "mm")
    , axis.text = element_text(size = 6)
    , panel.grid.major.x = element_blank()
    , axis.line = element_line(size = .25, colour = "black")
    , axis.ticks = element_line(size = .25, colour = 'black')
  )+
  guides(color=F,linetype=F)

p.GGP.vol.relpos.sensiana.g2 <- stats.GGP.vol.relpos %>% 
  select(tidy) %>% unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed' & !(threshold %in% c(2) & lab==30) & lab != 0) %>%
  ungroup() %>% 
  mutate(lab = as_labeller(GGP_names)(lab) %>% as.character()) %>% 
  mutate(term = dplyr::recode(term, `del.(Intercept)`='Δ', `del.logvol`='Δ~vol~')) %>% 
  mutate(term = forcats::fct_relevel(term, 'a', 'b', 'Δ', 'Δ~vol~')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = 'lightgrey') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size =.5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_grid(lab ~ relpos*term
             , scales = "free_x"
             #, ncol = 8
             ) +
  scale_x_continuous('Network density', breaks = seq(.1,1,.1))+
  scale_y_continuous('Estimate', breaks = pretty_breaks(2))+
  scale_color_aaas()+
  scale_linetype_manual(values = c('solid','11'))+
  theme(
    panel.spacing.x = unit(4, "mm")
  )+
  guides(color=F,linetype=F)


p.GGP.vol.relpos.sensiana <- change.facets.fcn(p.GGP.vol.relpos.sensiana.g1, p.GGP.vol.relpos.sensiana.g2)
ggsave('./../Word/Floats/Figures/S6-GGP-logvol-by_relpos-sensiana.tiff', plot = p.GGP.vol.relpos.sensiana
       , device = 'tiff', width = 18.5, height = 9.25, units = 'cm', dpi = 1200)
