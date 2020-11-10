source('./nlme.exp.r')

relposfcn <- function(side,Stroke.side){
  if(Stroke.side=="control") return("healthy");
  if((side=="Left" & Stroke.side=="L") | (side=="Right" & Stroke.side=="R")) return("ipsi");
  if((side=="Left" & Stroke.side=="R") | (side=="Right" & Stroke.side=="L")) return("contra");
}


## lookup table for plots and tables.
GGP_names <- c(
  '0'='q50',
  '1'="Efficiency",
  '2'="Clustering",
  '3'="Modularity",
  '4'="Assortativity",
  '5'="DiffusionEfficiency"
)

GGP_names_short <- c(
  '0'='q50',
  '1'="E",
  '2'="C",
  '3'="M",
  '4'="A",
  '5'="DE"
)

ggthemr('fresh')
old_swatch <- swatch()
new_swatch <- old_swatch
#new_swatch[1:9] <- new_swatch[c(2,4,3)]
new_swatch[1:9] <- new_swatch[c(6,2,4,3,5,1,7,8,9)]
new_swatch %>% as.list() %>% unlist %>% set_swatch()
#darken_swatch(0.5)

thm<-theme_update(text = element_text(size = 10)
                  , axis.text = element_text(size = 6)
                  , panel.grid.major.y = element_blank()
                  , panel.grid.major.x = element_line(colour = 'lightgrey', size = .05)
                  , axis.text.x = element_text(hjust = c(.75,.25,.5,.5))
                  , panel.background = element_rect(fill = "transparent",colour = NA))
sc_x_tp<-scale_x_continuous("Time since stroke"
                            , breaks = c(0,1,3,12)
                            , labels = c('3-5d','1m','3m','12m'))



reformatlabel <- function(s){
  tp<-strsplit(s,"[.]")[[1]][2]
  if(tp==0){f<-"3-5 d"}
  else if(tp==1){f<-"1 m"}
  else if(tp==3){f<-"3 m"}
  else if(tp==12){f<-"12 m"}
  return(f)
}


add_preds<-function(d,mdl){
  d$pred<-predict(mdl, type = 'link', se.fit = T)$fit
  d$pred.se<-predict(mdl, type = 'link', se.fit = T)$se.fit
  return(d)
}
plot_preds<-function(response,d,pred,cmn,mdl){
    s<-scale_y_continuous(cmn)
    xlab<-paste0('Change in ', pred)
  
    d$y = mdl$family$linkinv(d$pred)
    d$lower = mdl$family$linkinv(d$pred - 1.96*d$pred.se)
    d$upper = mdl$family$linkinv(d$pred + 1.96*d$pred.se)
    if (cmn == 'FM'){
      d$clin.meas.value <- 66 - d$clin.meas.value
      d$y = 66 - d$y
      d$lower = 66 - d$lower
      d$upper = 66 - d$upper
    }  
  p<-d %>% 
    ggplot(aes(x=get(response), y=clin.meas.value, color=as.factor(tp), shape = as.factor(tp)))+
    geom_point(alpha = .5, size = 1)+
    geom_line(aes(y = y, linetype = as.factor(tp))) +
    geom_line(aes(y = upper, linetype = as.factor(tp)), alpha=.32) +
    geom_line(aes(y = lower, linetype = as.factor(tp)), alpha=.32) +
    scale_y_continuous(if_else(pred=='q50',cmn,''))+
    scale_x_continuous(if_else(cmn == 'FM',paste0('Change in ', pred),''))+
        #scale_alpha_manual(values=c(.5,1))+
    scale_linetype_manual(values = c('solid','33','11'))+
    guides(alpha=F, color=F, shape=F, linetype=F) +
    theme(strip.background = element_blank(), strip.placement = "inside"
          , text = element_text(size = 8)
          , axis.text = element_text(size = 6)
          , axis.title = element_text(size = 8)
          , plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm')
          , axis.line = element_line(size = .25, colour = "black")
          , axis.ticks = element_line(size = .25, colour = 'black')
    ) +
    coord_cartesian(ylim = c(0,max(d$clin.meas.value, na.rm = T))) %>% 
    return()
}


joint_fit_CV<-function(response,d,pred,cmn,fam = 'quasipoisson'){
  
  #d<-subset(d, !(numID==40))
  
  
  fm <- lme(as.formula(paste0(response, '~ tp')), data = d, random = ~ 1 | numID, na.action = na.exclude)
  ffm <- fitted(fm)
  mdl<-glm(clin.meas.value ~ ffm, data = d, family = fam, na.action = na.exclude)
  
  d$ffm<-fitted(fm)
  
  newd <- data.frame(ffm=seq(min(ffm, na.rm = T),max(ffm, na.rm = T),length.out = 100))
  preds<-predict(mdl, type = 'link', se.fit = T, newdata = newd)
  newd$.fitted<-preds$fit
  newd$.se<-preds$se.fit
  
  newd$y <- mdl$family$linkinv(newd$.fitted)
  newd$lower <-  mdl$family$linkinv(newd$.fitted - 1.96*newd$.se)
  newd$upper <-mdl$family$linkinv(newd$.fitted + 1.96*newd$.se)
  
  if (cmn == 'FM'){
    d$clin.meas.value <- 66 - d$clin.meas.value
    newd$y = 66 - newd$y
    newd$lower = 66 - newd$lower
    newd$upper = 66 - newd$upper
  }
  plt<-d %>% 
    ggplot(aes(x=ffm,y=clin.meas.value, color=as.factor(tp), shape = as.factor(tp)))+
    geom_point(alpha=.5, size=1)+
    geom_line(aes(x = ffm, y = y), data = newd, color='black', inherit.aes = F)+
    geom_line(aes(x = ffm, y = lower), data = newd, color='black', linetype=2, inherit.aes = F)+
    geom_line(aes(x = ffm, y = upper), data = newd, color='black', linetype=2, inherit.aes = F)+
    scale_y_continuous(if_else(pred=='q50',cmn,''))+
    scale_x_continuous(if_else(cmn == 'FM',paste0('Change in ', pred),''))+
    coord_cartesian(ylim = c(0,max(d$clin.meas.value, na.rm = T)))+ 
    guides(color=F, shape=F)+
    theme(strip.background = element_blank(), strip.placement = "inside"
          , text = element_text(size = 8)
          , axis.text = element_text(size = 6)
          , axis.title = element_text(size = 8)
          , plot.margin = unit(c(0.1,0.1,0.1,0.1),'cm')
          , axis.line = element_line(size = .25, colour = "black")
          , axis.ticks = element_line(size = .25, colour = 'black')
    )
  
  return(plt)
}

joint_fit_CV_mdl<-function(response,d,pred,cmn,fam = 'quasipoisson', covs = ''){
  
  #d<-subset(d, !(numID==40))
  
  fm <- lme(as.formula(paste0(response, '~ tp')), data = d, random = ~ 1 | numID, na.action = na.exclude)
  mdl<-glm(as.formula(paste0('clin.meas.value ~ fitted(fm)', covs)), data = d, family = fam, na.action = na.exclude)
  
  
  return(mdl)
}


tbl.flex <- function(tbl,caption, fontsize=9, fontsizeheader=9){
  tbl %>% 
    flextable() %>% 
    flextable::theme_booktabs() %>%
    autofit() %>% 
    theme_zebra() %>% 
    align(align = 'left', part = 'footer')
}


tbl.flex.post <- function(tbl, caption, fontsize=9, fontsizeheader=9){
  tbl %>% 
    fontsize(size = fontsize, part = 'body') %>% 
    fontsize(size = fontsizeheader, part = 'header') %>% 
    padding(padding = 1) %>% 
    autofit() %>% 
    #add_footer_lines(value=caption) %>%  
    theme_zebra() %>% 
    align(align = 'center', part = 'header') %>% 
    align(align = 'left', part = 'footer')#%>% 
    #align(align = 'left', j = 1, part = 'body')
} 

tbl.kable <- function(tbl, caption){
  tbl %>% 
    kable(format = 'html', escape = F, caption = caption) %>% 
    kable_styling(bootstrap_options = c('striped','bordered','condensed')) %>% 
    collapse_rows(columns = 1)
}

pvalformatter <- Vectorize(
  function(x){
    if(!is.finite(x)){
      return('')
    }
    if(x==0){
      return(0)
    }
  if(x>0.0001){
    sprintf('%1.4f',x) %>% return()
  }else if(x<0.0001){
  sprintf('%1.2e',x) %>% return()
  }
  })

table.ref <- function(x) {
  stringr::str_extract(table_nums(x), "[^:]*")
}
figure.ref <- function(x) {
  stringr::str_extract(figure_nums(x), "[^:]*")
}
  



change.facets.fcn <- function(g1, g2){
  gt1 = ggplot_gtable(ggplot_build(g1))
  gt2 = ggplot_gtable(ggplot_build(g2))
  gt1 <- gtable_add_rows(gt1, heights = unit(0.5, 'cm'), pos = 2)
  gt1 <- gtable_add_grob(gt1, grobs = gt2$grobs[grep('strip-t', gt2$layout$name)], t = 2, l = gt1$layout[grep('strip-t.+1$', gt1$layout$name),]$l)
  
  
  #gt1$grobs[grep('strip-t.+1$', gt1$layout$name)] = gt2$grobs[grep('strip-t', gt2$layout$name)]
  grid.draw(gt1)
  
  gt1 = gtable_add_cols(gt1, widths=gt1$widths[1], pos = -1)
  
  panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
  gt.side1 = gtable_filter(gt2, 'strip-r-1')
  gt.side2 = gtable_filter(gt2, 'strip-r-2')
  gt.side3 = gtable_filter(gt2, 'strip-r-3')
  gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))
  gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))
  
  return(gt1 %>% as_ggplot())
}

format.fcn <- Vectorize(
  function(x){
    if(!is.numeric(x)) return('bla')
    if(is.na(x)) return('na')
    if(x==0)return(0)
    if (abs(x)<0.0000001){
      return(sprintf('%.0e',x))
    }
    else{
      return(x)
    }
  }
)
