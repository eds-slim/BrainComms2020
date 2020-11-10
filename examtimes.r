source('./loaddata.r')

data.time <- read.csv('./../SFB_structural_connectomes/SFB_patient_list_31082017_mod.csv', header = T)

d.time.acute <- data.time %>% 
  select(c(ID,T0,T1.exam, T2.exam,T3.exam, T4.exam)) %>% 
  mutate(T0 = as.Date(T0, format = '%d/%m/%Y')) %>% 
  mutate(T1.exam = as.Date(T1.exam, format = '%d/%m/%Y')) %>% 
  mutate(T1.delta = difftime(T1.exam, T0, units = 'days'))%>% 
  summarise(mean = mean(T1.delta, na.rm=T)
            , sd = sd(T1.delta, na.rm=T)
            , median = median(T1.delta, na.rm=T)
            , lower = quantile(T1.delta, na.rm=T, probs = 0.25)
            , upper = quantile(T1.delta, na.rm=T, probs = 0.75))

d.time <- data.time %>% 
  select(c(ID,T0,T1.exam, T2.exam,T3.exam, T4.exam)) %>% 
  mutate(T0 = as.Date(T0, format = '%d/%m/%Y')) %>% 
  mutate(T1.exam = as.Date(T1.exam, format = '%d/%m/%Y')) %>% 
  mutate(T2.exam = as.Date(T2.exam, format = '%d/%m/%Y')) %>% 
  mutate(T3.exam = as.Date(T3.exam, format = '%d/%m/%Y')) %>% 
  mutate(T4.exam = as.Date(T4.exam, format = '%d/%m/%Y')) %>% 
  filter(ID %in% data.clinical$ID) %>% 
  #mutate(T1.delta = difftime(T1.exam, T0, units = 'days'))%>% 
  mutate(T2.delta = difftime(T2.exam, T0, units = 'weeks'))%>% 
  mutate(T3.delta = difftime(T3.exam, T0, units = 'weeks'))%>% 
  mutate(T4.delta = difftime(T4.exam, T0, units = 'weeks')) %>% 
  gather(time, val, T2.delta:T4.delta) %>% 
  group_by(time) %>% 
  summarise(mean = mean(val, na.rm=T)
            , sd = sd(val, na.rm=T)
            , median = median(val, na.rm=T)
            , lower = quantile(val, na.rm=T, probs = 0.25)
            , upper = quantile(val, na.rm=T, probs = 0.75))

d.time[1,'median']
