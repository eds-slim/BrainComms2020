data<-data.clinical_long

data$FM_norm<-data$FM/66

data$comp_motor<-rowMeans(data[c("Gripstrength_ratio","NHP_ratio","FM_norm")], na.rm = T)


idx.minimal<-subset(data, comp_motor>.6 & tp==0)$numID

idx.good<-subset(data, !numID%in%idx.minimal & comp_motor>.6 & tp==3)$numID

idx.poor<-subset(data, !numID%in%idx.minimal & comp_motor<.6 & tp==3)$numID


data.clinical_long$recovery<-NA
data.clinical_long$recovery[data.clinical_long$numID%in%idx.minimal]<-"minimal"
data.clinical_long$recovery[data.clinical_long$numID%in%idx.good]<-"good"
data.clinical_long$recovery[data.clinical_long$numID%in%idx.poor]<-"poor"
data.clinical_long$recovery<-as.factor(data.clinical_long$recovery)


