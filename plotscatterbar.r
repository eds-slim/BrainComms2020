d.plot <- dd.delta %>% 
  mutate(ROI = forcats::fct_reorder(ROI,p)) %>% 
  unnest(data) %>% 
  mutate(ROI = forcats::fct_reorder(ROI,as.numeric(lobe))) %>% 
  filter(ROI %in% levels(d$ROI)[1:42]) %>% 
  dplyr::select(-c(lesionvolume,llv, delta.nemo)) 

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  #l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("^(.*)e", "e", l)
  # turn the 'e+' into plotmath format
  #l <- gsub("e", "%*%10^", l)
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

main.plot <- d.plot %>% 
  filter(nemoscore>0) %>% 
  spread(visit,nemoscore) %>% 
  ggplot(aes(x=V0, y=V3, color = as.factor(lobe), shape = treatment))+
  scale_x_continuous('ChaCo score at V0', trans = 'log10', labels=fancy_scientific) + # function(x)sprintf('%1.0e',x)
  scale_y_continuous('ChaCo score at V3', trans = 'log10', labels=fancy_scientific) +
  geom_point(alpha=.2) + 
  coord_fixed() + 
  geom_rug(alpha = .05, color = 'black') +
  stat_function(data=data.frame(x=c(1e-10,1)), fun = function(x)x, color='black', linetype = 3, inherit.aes = FALSE)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(.~ROI, scales = 'fixed', ncol = 6) + 
  theme(text = element_text(size=7)) + 
  guides(color=FALSE, shape =FALSE)


## A function to plot the inset 
get_inset <- function(df){
  p <- df %>% 
    group_by(ROI,visit, treatment) %>% 
    mutate(proppos = mean(nemoscore>0)) %>% 
    ggplot(aes(x = treatment, y = proppos, fill = visit, group = visit))+
    geom_col(position = position_dodge(width = 0.6), width = 0.4, color = 'darkgrey', size = .1, alpha = .2)+
    geom_text(data = . %>% ungroup()  %>% group_by(ROI, treatment) %>% mutate(max.proppos = max(proppos)) %>% filter(visit == 'V0') %>% ungroup() %>% mutate(label = forcats::fct_recode(treatment, A='rtPA', P='Placebo'))
                , aes(x = treatment, y = max.proppos, label = label), vjust = -0.5, size = 2, inherit.aes = FALSE) + 
    geom_segment(aes(color = lobe), x = .4, xend = .4, y = 0, yend = 1, arrow = arrow(length = unit(0.1, "cm"), type = 'closed', angle = 90, ends = 'both')) +
    annotate("text", x = -.4, y = 1, label = "1", hjust = -1, vjust = 1, size = 2) +
    annotate("text", x = -.4, y = 0, label = "0", hjust = -1, vjust = 0, size = 2) +
    theme_minimal() + 
    guides(fill=FALSE, color = FALSE) + xlab('') + 
    scale_y_continuous('', breaks = NULL, limits = c(0,1.3), expand = c(0,0)) + 
    scale_x_discrete(name = '', breaks = NULL) + scale_fill_grey(start = .8, end = .9) + scale_color_discrete(drop = FALSE) +
    facet_wrap(.~ROI, scales = 'free')+
    theme(text = element_text(size = 6)
          , strip.background = element_blank()
          , strip.text.x = element_blank()
          , panel.grid = element_blank()
          )
  return(p)
}
inset_plot <- get_inset(d.plot) 
inset_plot
main.plot + annotation_custom(grob = ggplotGrob(inset_plot)
                              , xmin = log(1e-5), ymin = log(1e-3), xmax = log(1), ymax = log(1e-1))

## This function allows us to specify which facet to annotate
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity
        , geom = ggplot2:::GeomCustomAnn,
        inherit.aes = FALSE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax
                                          ))
}


main.plot + annotation_custom2(grob=ggplotGrob(inset_plot), 
                               data = d.plot %>% subset(ROI=='Amygdala'),
                               ymin = -Inf, ymax=Inf, xmin=-Inf, xmax=Inf)


insets <- d.plot %>% droplevels() %>% 
  split(f = .$ROI) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)), 
    data = data.frame(.)
    , xmin = log(7e-2), ymin = log(3e-4), xmax = log(2e0), ymax = log(3e-2))
  )

p <- main.plot + insets
ggsave('./../../derivatives/figures/scatterplot.tiff', plot = p, device = 'tiff'
       , width = 21, height = 29.7, units = 'cm', dpi = 1200)
