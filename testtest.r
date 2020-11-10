dd<-d.LGP %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS' 
         & relpos=='ipsi' 
         & node == 'postcentral'
         &LGP.name == 'efficiency'
         & tp>=0) %>% as_tibble()


d3.lag1<-dd %>% 
  group_by(numID,clin.meas.name,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(LGP.value.delta = LGP.value - dplyr::lag(LGP.value, n = 1, default = NA)
               # ,clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 1, default = NA)
                ,tp.mean = (tp + dplyr::lag(tp, n = 1, default = NA))/2) %>% 
  dplyr::select(-tp)

d3.lag2<-dd %>% 
  group_by(numID,clin.meas.name,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(LGP.value = LGP.value - dplyr::lag(LGP.value, n = 2, default = NA)
                ,clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 2, default = NA)
                ,tp = (tp + dplyr::lag(tp, n = 2, default = NA))/2)
d3.lag3<-dd %>% 
  group_by(numID,clin.meas.name,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(LGP.value = LGP.value - dplyr::lag(LGP.value, n = 3, default = NA)
                ,clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 3, default = NA)
                ,tp.mean = (tp + dplyr::lag(tp, n = 3, default = NA))/2)

d3<-bind_rows(d3.lag1) %>% 
  filter(!is.na(tp)) %>% 
  bind_rows(dd) %>% 
  arrange(numID,tp)

dplyr::full_join(dd,d3.lag1) %>% 
  filter(numID==24) %>% 
  dplyr::select(c(tp,tp.mean))

d3.lag1 %>% 
  filter(tp==12) %>% 
  ggplot(aes(x=LGP.value.diff, y=clin.meas.value, color=as.factor(tp), shape=as.factor(tp.mean)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)

d3 %>% 
  filter(tp%%1==0.5 | tp==2) %>% 
  ggplot(aes(x=LGP.value,y=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(aes(group=tp), method = 'lm', se=F)+
  geom_smooth(aes(group=1), method = 'lm', se=F, color='green')


d3 %>% 
  filter(tp%%1==0.5 | tp==2) %>% 
  lme(clin.meas.value~LGP.value, random = ~1|numID, data = ., na.action = na.omit, correlation = corCAR1(form=~tp)) %>% 
  tidy()

d3 %>% 
  filter(tp%%1==0.5 | tp==2) %>% 
  rmcorr(participant = tp, measure1 = clin.meas.value, measure2 = LGP.value, dataset = .)

p.list<-c()
beta.list<-c()
for( t1 in unique(d3$tp)){ #q50
  for( t2 in unique(d3$tp)){ # clinical
    if(t1>t2){
      p.list<-c(p.list,NA)
      beta.list<-c(beta.list,NA)
      next()
      }
    x<-d3 %>% filter(tp==t1) %>% select(-c(tp,clin.meas.value))
    y<-d3 %>% filter(tp==t2) %>% select(-c(tp,conn.meas.value))
    df<-merge(x,y) %>% as_tibble()
    m<-df %>% 
      group_by(clin.meas.name,relpos) %>% 
      nest() %>% 
      mutate(mdl=map(data,~lm(clin.meas.value~conn.meas.value,data=., na.action = na.omit))) %>% 
      mutate(tidy=map(mdl,tidy))
    p<-m %>% 
      unnest(tidy) %>% 
      filter(term=='conn.meas.value') %>% 
      pull('p.value')
    p.list<-c(p.list,p)
    beta<-m %>% 
      unnest(tidy) %>% 
      filter(term=='conn.meas.value') %>% 
      pull('estimate')
    beta.list<-c(beta.list,beta)
    
    next()
    m<-lm(y~x,df,na.action = na.omit) %>% tidy
    p<-m%>% filter(term=='x') %>% pull('p.value')
    beta<-m%>% filter(term=='x') %>% pull('estimate')
    if(p<0.05){sprintf('t1=%f, t2=%f; p=%f, beta=%f',t1,t2,p,beta) %>% print()}
  } 
}
require(gplots)
matrix(p.list, ncol = length(unique(d3$tp))) %T>% print() %>% 
  apply(2, function(y) as.numeric(y > 0.05)) %>% 
  heatmap.2(Rowv = NA, Colv = NA, dendrogram = 'none', trace='none', key=F
            , labRow = unique(d3$tp), labCol = unique(d3$tp)
            , symm=F,symkey=F,symbreaks=T, scale="none")

matrix(beta.list, ncol = length(unique(d3$tp))) %T>% print() %>% 
 # apply(2, function(y) as.numeric(y > 0.05)) %>% 
  heatmap.2(Rowv = NA, Colv = NA, dendrogram = 'none', trace='none', key=T
            , labRow = unique(d3$tp), labCol = unique(d3$tp)
            , symm=F,symkey=F,symbreaks=T, scale="none")
