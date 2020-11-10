p.GGP.timecourse.sensiana.g1<-stats.GGP %>% 
  unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed' & lab != 0) %>%
  mutate(term = dplyr::recode(term, `del.(Intercept)`='Δ', `del.relposipsi`='Δ~ipsi~')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color='lightgray') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size=.5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_wrap(lab~term
             , scale = "free_x"
             , ncol = 4
             , labeller = labeller(lab = GGP_names)) +
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
p.GGP.timecourse.sensiana.g1

p.GGP.timecourse.sensiana.g2<-stats.GGP %>% 
  unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed' & lab != 0) %>%
  mutate(term = dplyr::recode(term, `del.(Intercept)`='Δ', `del.relposipsi`='Δ~ipsi~')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color='lightgray') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size=.5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_grid(lab~term
             , scales = "free_x"
             , labeller = labeller(lab = GGP_names)) +
  scale_x_continuous('Network density', breaks = seq(.1,1,.1))+
  scale_y_continuous('Estimate', breaks = pretty_breaks(3))+
  scale_color_aaas()+
  scale_linetype_manual(values = c('solid','11'))+
  theme(panel.spacing.x = unit(4, "mm")
  )+
  guides(color=F,linetype=F)

p.GGP.timecourse.sensiana.g2

p.GGP.timecourse.sensiana <- change.facets.fcn(p.GGP.timecourse.sensiana.g1, p.GGP.timecourse.sensiana.g2)
p.GGP.timecourse.sensiana
ggsave('./../Word/Floats/Figures/S3-GGP-sensiana.tiff', plot = p.GGP.timecourse.sensiana
       ,device = 'tiff', width = 18.5, height = 9.25, units = 'cm', dpi = 1200)


p.GGP.timecourse.relpos.sensiana.g1<-stats.GGP.relpos %>%
  unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed' & lab != 0) %>%
  mutate(term = dplyr::recode(term, `del`='Δ')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = 'lightgray') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size = .5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_wrap(lab*relpos~term
             , scale = "free_x"
             , ncol = 6
             , labeller = labeller(lab = GGP_names)) +
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

p.GGP.timecourse.relpos.sensiana.g2<-stats.GGP.relpos %>% 
  unnest(tidy, .drop = T) %>% 
  filter(effect=='fixed' & lab != 0) %>%
  mutate(term = dplyr::recode(term, `del`='Δ')) %>% 
  ggplot(aes(x = threshold, y = estimate, linetype=conf.high*conf.low<0, color=as.factor(lab))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = 'lightgray') +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = .1, size = .5) +
  geom_point(size = .5) +
  coord_flip() +
  facet_grid(lab ~ relpos*term
             , scales = "free_x"
             #, ncol = 6
             , labeller = labeller(lab = GGP_names)) +
  scale_x_continuous('Network density', breaks = seq(.1,1,.1))+
  scale_y_continuous('Estimate', breaks = pretty_breaks(3))+
  scale_color_aaas()+
  scale_linetype_manual(values = c('solid','11'))+
  theme(
    panel.spacing.x = unit(4, "mm")
  )+
  guides(color=F,linetype=F)


p.GGP.timecourse.relpos.sensiana <- change.facets.fcn(p.GGP.timecourse.relpos.sensiana.g1, p.GGP.timecourse.relpos.sensiana.g2)

ggsave('./../Word/Floats/Figures/S4-GGP-relpos-sensiana.tiff', plot = p.GGP.timecourse.relpos.sensiana
       ,device = 'tiff', width = 18.5, height = 9.25, units = 'cm', dpi = 1200)
