mean(data.clinical$age)
sd(data.clinical$age)

xtabs(~gender, data=data.clinical) # 1=male
prop.test(count(data.clinical, vars="gender")[1,2],nrow(data.clinical), p=.5)

plyr::count(data.patients$Stroke.side)
prop.test(plyr::count(data.patients, vars="Stroke.side")[1,2],nrow(data.clinical), p=.5)       

quantile(data.clinical$vol, na.rm = TRUE)


LM<-lm(log(vol)~Stroke.side, data=data.clinical)
summary(LM)

t.test(data.clinical$age[data.clinical$Stroke.side=="L"],data.clinical$age[data.clinical$Stroke.side=="R"], var.equal = T)

xtabs(~Stroke.side+gender, data=data.clinical) %>% chisq.test()

quantile(data.clinical$NIHSS_T1,na.rm = TRUE)
wilcox.test(data.clinical$NIHSS_T1[data.clinical$Stroke.side=="L"],data.clinical$NIHSS_T1[data.clinical$Stroke.side=="R"], var.equal = T)

LM<-lm(NIHSS_T1~gender, data=data.clinical)
summary(LM)
t.test(data.clinical$NIHSS_T1[data.clinical$gender=="male"],data.clinical$NIHSS_T1[data.clinical$gender=="female"], var.equal = T)


quantile(data.clinical$rGS_T4,na.rm = TRUE)
wilcox.test(data.clinical$rGS_T1[data.clinical$Stroke.side=="L"],data.clinical$rGS_T1[data.clinical$Stroke.side=="R"], var.equal = T)

LM<-lm(NIHSS_T4~log(vol)+Stroke.side, data=data.clinical)
summary(LM)
t.test(data.clinical$rGS_T1[data.clinical$gender=="male"],data.clinical$rGS_T1[data.clinical$gender=="female"], var.equal = T)



quantile(data.clinical$NHP_ratio_T4,na.rm = TRUE)




## two-way
dd<-merge(data.patients[,c("ID","Stroke.side")],data.clinical[,c("ID","vol","age","gender","NIHSS_T1","FM_T1","NHP_ratio_T1","rGS_T1","dominant_affected")])

dd<-data.clinical[,c("numID","Stroke.side","vol","age","gender","NIHSS_T1","FM_T1","NHP_ratio_T1","rGS_T1","dominant_affected")]
dd<-unique(dd)
dd$NHP_ratio_T1<-1/dd$NHP_ratio_T1
dd$logvol<-log(dd$vol)
str(dd)


#tiff("contvars.tiff",res=600, width=7,height=7, units="in")
pal<-pal_jama()(4)
pairs.panels(dd[,c(4,11,6,7,8,9)], 
             method = "pearson", # correlation method
             ci=TRUE,
             hist.col = pal[3],
             density = TRUE,  # show density plots
             ellipses = FALSE, # show correlation ellipses
             lm = TRUE,
             stars = TRUE,
             bg = pal[dd$gender],
             pch=21,
             alpha=.05,
             col=pal[4],
             cex.labels=1.2,
             labels=c("Age","log volume","NIHSS","Fugl-Meyer","NHP ratio","Grip strength ratio")
             )
#dev.off()

dd$gender<-as.factor(dd$gender)

dd$dominant_affected<-as.factor(dd$dominant_affected)
levels(dd$dominant_affected)<-c("no","yes")

dd$NHP_ratio_T1<-1/dd$NHP_ratio_T1
dd.long<-gather(dd,key="cont.key", value="cont.value",age, logvol,NIHSS_T1,FM_T1,NHP_ratio_T1,rGS_T1)
dd.long<-gather(dd.long,key="disc.key", value="disc.value", Stroke.side,gender, dominant_affected)

dd.long$cont.key = factor(dd.long$cont.key,levels(factor(dd.long$cont.key))[c(1,4,6,2,3,5)])

#tiff("disc-con.tiff",res=600, width=10,height=7, units="in")
p<-ggplot(dd.long, aes(x=disc.value, y=cont.value, fill=disc.value))+
  geom_violin(alpha=.1)+
  geom_boxplot(alpha=.1, width=.1)+
  geom_point(position=position_jitter(width=.1), alpha=1, shape=23, color="black")+
  scale_fill_jama()+scale_color_jama()+
  facet_grid(cont.key~disc.key, scales="free", switch="y", labeller = label_value)+
  scale_x_discrete("")+scale_y_continuous("", expand=c(0.5,0))+
  theme_bw()+
  guides(fill=FALSE)
print(p)
#dev.off()

xtabs(~dominant_affected+Stroke.side, data=dd) %>% chisq.test()

