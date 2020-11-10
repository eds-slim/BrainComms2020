require(gdata)
require(tibble)
require(plyr)
require(tidyr)
require(magrittr)
require(dplyr)
require(purrr)
require(broom)
require(broom.mixed)


require(foreach)
require(effsize)

require(ggplot2)
require(ggsignif)
library(ggsci)
require(ggpubr)
require(ggthemr)
require(ggrepel)


require(lmerTest)
require(nlme)
require(nls.multstart)
require(stringr)

require(gridExtra)
require(kableExtra)
require(scales)
require(flextable)


library(grid)
library(gtable) 

source('./auxfcns.r')

data.overview<-read.xls('../SFB_structural_connectomes/SFB_structural_connectomes_data_overview.xlsx', sheet = 1, header = TRUE, pattern = "ID", na.strings=c("","NA"))
str(data.overview)
data.overview<-within(data.overview,
  lesion.volume..ml.<-as.numeric(as.character(lesion.volume..ml.)))

data.overview$numID <- seq.int(nrow(data.overview))

data.patients<-data.overview[!is.na(data.overview$Stroke.side),c('numID','ID', 'Stroke.side', 'lesion.volume..ml.')]

str(data.patients)


data.clinical<-read.xls("../SFB_structural_connectomes/clinical_data/clinical_data_structural_connectomes.xlsx")
str(data.clinical)

data.clinical<-within(data.clinical,{
  gender<-as.factor(gender)
  levels(gender)<-c("female","male")
})


data.clinical<-merge(data.patients, data.clinical)
data.clinical$vol<-coalesce(data.clinical$lesion.volume..ml.,data.clinical$volume)
data.clinical<-within(data.clinical,rm(lesion.volume..ml.,volume,X,X.1,X.2,X.3))
data.clinical$logvol <- log(data.clinical$vol)

## long format
data.clinical$initialNIHSS <- data.clinical$NIHSS_T1
data.clinical$initialrGS <- data.clinical$Gripstrength_ratio_T1
data.clinical$initialFM <- data.clinical$FM_T1

data.clinical_long<-reshape(data.clinical, direction="long",varying = 7:38, timevar="tp", idvar="numID", sep="_T")
data.clinical_long$tp<-as.numeric(revalue(as.character(data.clinical_long$tp),c("1"="0","2"="1","3"="3","4"="12")))

data.clinical_long$NHP_ratio<-1/data.clinical_long$NHP_ratio
#data.clinical_long$FM<-66 - data.clinical_long$FM


source("patientsubgroups.r")

data.clinical_long<-data.clinical_long %>% 
  gather('clin.meas.name','clin.meas.value',c('NIHSS','FM','NHP_ratio','NHP.affected_hand','NHP.not_affected_hand','Gripstrength_ratio','Gripstrenght_affected_hand','Gripstrenght_not_affected_hand')) %>%
  filter(!stringr::str_detect(clin.meas.name,'NHP') & !str_detect(clin.meas.name,'not_affected')) %>% 
  mutate(clin.meas.name = dplyr::recode(clin.meas.name, Gripstrength_ratio = 'rGS'))



## hemisphere data

data.hemisphere<-merge(data.overview[,"numID", drop=FALSE], data.patients[,c("numID","Stroke.side")], all=TRUE)
levels(data.hemisphere$Stroke.side)<-c("L","R","control")
data.hemisphere$Stroke.side[which(is.na(data.hemisphere$Stroke.side))]<-"control"
data.hemisphere<-data.hemisphere[rep(row.names(data.hemisphere),each=2),]
data.hemisphere$side<-c("Left","Right")
data.hemisphere$relpos <- as.factor(mapply(relposfcn,data.hemisphere$side, data.hemisphere$Stroke.side))

data.hemisphere$relpos <- relevel(data.hemisphere$relpos,"healthy")
data.hemisphere$relpos<-factor(data.hemisphere$relpos,levels=c("healthy","contra","ipsi"))


ID.small <- c(3, 5, 7,13,14,15,20,22,23,24,26,25,30,32,33,35,44,45,46,47,51,55,53) # definitely include
ID.large <- c(9,12,16,36,43,52,54)

ID.missing <- c(13,23,24,53,54,55)

ID.odd <- c(34)
ID.bs <- c(31,41,50) # definitely exclude


ID.keep <-union(ID.small,ID.large) 
#ID.keep <-union(ID.keep, 34) 
#ID.keep <- union(ID.keep, ID.bs)
#ID.keep <- setdiff(ID.keep, ID.missing)

numID.keep <- data.overview$numID[data.overview$ID %in% ID.keep & data.overview$X=='patient (p)']

numID.exclude<-c(37, 40, 41,43, 49) #p31 (bs), ?p34 (odd shape), p35 (lesion outside brain), p41 (bs), p50 (bs)

data.patients<-subset(data.patients, numID %in% numID.keep)
data.clinical<-subset(data.clinical, numID %in% numID.keep)
data.clinical_long<-subset(data.clinical_long, numID %in% numID.keep)
data.hemisphere<-subset(data.hemisphere, numID %in% numID.keep)

idx.complete <- data.clinical_long %>% filter(clin.meas.name != 'Gripstrenght_affected_hand') %>% ungroup() %>%
  dplyr::select(c(ID,tp,clin.meas.name,clin.meas.value)) %>% 
  unite(temp,tp,clin.meas.name) %>% 
  spread(temp,clin.meas.value) %>% 
  dplyr::select(-ID) %>% 
  complete.cases()

#####
# time course relative to initial level
data.clinical_long %<>% 
  group_by(numID,clin.meas.name) %>% 
  arrange(tp) %>% 
  mutate(clin.meas.value.diff=clin.meas.value-clin.meas.value[tp==0]) %>% 
  arrange(clin.meas.name)

###############################
### load connectivity data ####
###############################

data.conn<-read.csv('./../q50-full.dat', header = TRUE)
data.conn$side<-c("Left","Right")[data.conn$side]

conn.measures<-c("q50","lambda")
data.conn<-merge(data.conn, data.hemisphere)

d.conn <- data.conn %>% 
  gather(key = 'lab', value = 'GGP', c('q50','lambda')) %>%
  filter(relpos != 'healthy', lab == 'q50') %>% droplevels() %>% 
  mutate(lab=0)
d.conn$threshold <- 1




###############################
### load GGP data ####
###############################
data.GGP<-read.csv('./../structural-proportional-intra-full-sum/intermediate/GGP-raw.csv', header = TRUE)
data.GGP$side<-c("Left","Right")[data.GGP$side]

d.GGP<-merge(data.GGP, data.hemisphere)
d.GGP <- rbind(d.GGP,d.conn)

d.GGP %<>% subset(lab %in% c(0,1,3))
d.GGP %<>% subset(relpos !='healthy') %>% droplevels()


## GGP ipsi-contra
# d.GGP %<>% 
#   group_by(numID,tp,lab, threshold) %>% 
#   arrange(relpos) %>% 
#   mutate(GGP=GGP[relpos=='ipsi']-GGP[relpos=='contra']) %>% 
#   filter(relpos=='ipsi') %>% 
#   mutate(relpos='delta') %>% 
#   bind_rows(d.GGP)

## GGP relative to t0
d.GGP %<>% 
  group_by(numID,relpos,lab, threshold) %>% 
  arrange(tp) %>% 
  mutate(GGP.diff=GGP-GGP[tp==0])


###############################
### load LGP data ####
###############################
data.LGP<-read.csv('./../structural-proportional-intra-full-sum/intermediate/LGP.csv', header = TRUE)
data.LGP$side<-c("Left","Right")[data.LGP$side]

data.LGP<-merge(data.LGP, data.hemisphere)
data.LGP$relpos<-factor(data.LGP$relpos,levels=c("healthy","ipsi","contra"))

labels <- read.table("../ROIs.dat", header = FALSE)
data.LGP$node <-labels$V1[data.LGP$node]

d.LGP <- data.LGP %>% 
  gather(key = 'LGP.name', value = 'LGP.value', c('strength','efficiency','clustering','nodelambda','betweenness')) %>%
  filter(relpos != 'healthy') %>% droplevels()


## LGP ipsi-contra
# d.LGP %<>% 
#   group_by(numID,tp,LGP.name,node) %>% 
#   arrange(relpos) %>% 
#   mutate(LGP.value=LGP.value[relpos=='ipsi']-LGP.value[relpos=='contra']) %>% 
#   filter(relpos=='ipsi') %>% 
#   mutate(relpos='delta') %>% 
#   bind_rows(d.LGP)

## LGP relative to t0
d.LGP %<>% 
  group_by(numID,relpos,LGP.name, node) %>% 
  arrange(tp) %>% 
  mutate(LGP.value.diff=LGP.value-LGP.value[tp==0]) 

######################
deriv_mod <- deriv(~a+del*(1-exp(-b*tp)), c('tp','a','b','del'),function(tp,a,b,del){})

