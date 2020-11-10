###### experimental ###############

d.conn %>% 
  merge(data.clinical_long) %>% 
  group_by(conn.meas.name, clin.meas.name, relpos) %>% 
  nest() %>% 
  mutate(mdl=map(data,~lme(clin.meas.value~conn.meas.value.diff, random=~1|tp, data=., na.action = na.omit))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(effect=="fixed" & term == "conn.meas.value.diff")

d %>% 
  lme(clin.meas.value~conn.meas.value.diff, random=~1|tp, data=., na.action = na.omit) %>% 
  tidy()

require(MuMIn)
d %>%  
  ggplot(aes(x=clin.meas.value, y=conn.meas.value.diff,  color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=1), method = 'lm', se = F, color='green')+
  facet_wrap(.~clin.meas.name, scales = 'free')



stats<-d %>% 
  merge(data.clinical_long) %>% 
  group_by(tp,clin.meas.name,relpos) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(lm)(clin.meas.value~conn.meas.value
                                  , na.action = na.omit
                                  , data=.
  )$result
  )
  ) %>% 
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))

stats %>% 
  unnest(tidy, .drop = T) %>% 
  filter(p.value<0.05)



df1<-tibble(t1=c(0,1,3,12),LGP=c("x","y","z","w"),NIHSS=c("a","b","c","d"))
df2<-tibble(t2=c(0,1,3,12),NIHSS=c("a","b","c","d"))


cbind(df1,df2) %>% 
  group_by(t1,t2) %>% 
  nest()
dd<-d %>% 
  merge(data.clinical_long)
for(t1 in c(0,1,3,12)){
  for(t2 in c(0,1,3,12)){
    x<-dd$clin.meas.value[dd$relpos=='ipsi' & dd$tp==t1 & dd$clin.meas.name=='NIHSS']
    y<-dd$conn.meas.value[dd$relpos=='ipsi' & dd$tp==t2 & dd$clin.meas.name=='NIHSS']
    corr.test(x,y) %>% print()
  }
}
dd<-d %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS' & relpos=='ipsi') %>% 
  as_tibble()

expand.grid(c(0,1,3,12),c(0,1,3,12)) %>% as_tibble() %>% group_by(Var1,Var2) %>% 
  nest() %>% 
  mutate(d1=map(Var1,~subset(dd,tp==.)),d2=map(Var2,~subset(dd,tp==.))) %>% 
  mutate(mdl=pmap(list(d1,d2),~lm(.x$clin.meas.value~.y$conn.meas.value, na.action = na.omit))) %>% 
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe())) %>% 
  unnest(tidy) %>% 
  filter(term==".y$conn.meas.value")

dd<-d %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='rGS' & relpos=='ipsi')

dd %>% 
  lme(clin.meas.value~exp(-tp)*conn.meas.value, random=~1|numID, data=., na.action = na.omit) %>% 
  tidy()
