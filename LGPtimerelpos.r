




plot_fcn<-function(d, p, node, LGP.name){
  if(p>0.05){return(zeroGrob())}
  if(LGP.name=='clustering'){
    title<-node
  }else{
    title<-node
  }
  if(node=='caudalanteriorcingulate'){
    ylab<-LGP.name
  }else{
    ylab<-''
  }
  # run nls_multstart
  fits <- d %>% 
    group_by(relpos) %>% 
    nest() %>% 
    mutate(fit=map(data,~nls_multstart(LGP.value ~ a+del*(1-exp(b*tp)),
                                       data = .,
                                       iter = 500,
                                       start_lower = c(a = 1, b = 0.2, del = -30),
                                       start_upper = c(a = 300, b = 1.5, del = -1),
                                       supp_errors = 'Y',
                                       na.action = na.omit)))
  
  print(fits)
  # new data frame of predictions
  new_preds <- d %>%
    do(., data.frame(tp = seq(0,12, length.out = 10), stringsAsFactors = FALSE))
  
  # create new predictions
  # preds2 <- fits %>%
  #   unnest(fit %>% map(augment, newdata = new_preds)) %>%
  #   rename(., LGP.value = .fitted) %>%
  #   ungroup()
  
  preds2 <- fits %>%
    mutate(fit = map(fit, ~augment(., newdata = new_preds))) %>% 
    select(fit) %>% unnest(fit) %>% 
    rename(., LGP.value = .fitted) %>%
    ungroup()
  
  # plot
  p<-ggplot(d,aes(x=tp-.15+as.numeric(relpos)/10, y=LGP.value, col = relpos, fill=relpos)) +
    #geom_point(size = 2)+
    geom_line(aes(group = relpos, y=LGP.value), data=preds2)+
    stat_summary(aes(group = relpos), geom = "point", fun.y = mean, shape = 21, size = 1.5, color='black')+
    stat_summary(aes(group = relpos, width=.75), size=.5, fun.data = mean_se, geom = "errorbar")+
    theme_pubr()+
    scale_x_continuous("Time post stroke", breaks = c(0,1,3,12), labels = c('3-5d','1m','3m','12m'))+
    scale_y_continuous(ylab)+
    thm+theme(text = element_text(size = 6))+
    ggtitle(title)+
    guides(color=F, fill=F)

  return(p)
}

b<-stats.LGP.ipsi %>% filter(relpos=='ipsi') %>% 
  filter(term=='del') %>% 
  dplyr::select(c(LGP.name,node,p.value))

a<-d.LGP %>% 
  filter(node %in% unique(nodes$node) & LGP.name %in% c('strength','efficiency','clustering')) %>% 
  merge(b) %>% 
  group_by(node,LGP.name) %>% 
  arrange(LGP.name,node) %>% 
  nest() %>% 
  mutate(p.value=map_dbl(data,~.$p.value[[1]])) %>% 
  mutate(plot=pmap(list(data,p.value,node, LGP.name),plot_fcn))

#tiff('./../Word/Floats/Figures/07-LGP-smallmults.tiff', width = 18.5, height = 37, units = 'cm', res = 600)
p.LGP.timecourse<-do.call(grid.arrange,c(a$plot,ncol=3, as.table=F)) %>% as_ggplot()
#p.LGP.timecourse
#dev.off()




