gc_data<-read.csv("/Users/jongewirtzman/Downloads/Gewirtzman_Oak_Samples.csv")

str(gc_data)
plot(gc_data$X.CH4.[which(gc_data$Sample.Type=="Lab Air Blank")])
plot(log(gc_data$X.CH4.))


library(ggplot2)


##

gc_data.1<-gc_data[which(is.na(gc_data$Tree.Tissue)==F),]
gc_data.1<-gc_data.1[-which(gc_data.1$Sample.Type=="Lab Air Blank"),]
gc_data.1<-gc_data.1[-which(gc_data.1$Sample.Type=="Aerobic (UZA)"),]
gc_data.1<-gc_data.1[-which(gc_data.1$Sample.Type=="Check Standard"),]
gc_data.1<-gc_data.1[-which(is.na(gc_data.1$Tree.Height)==T),]

ggplot(data=gc_data.1,
       aes(x=as.numeric(Tree.Height), y=X.CH4., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()+ggtitle("Anaerobic Incubation CH4")

ggplot(data=gc_data.1,
       aes(x=as.numeric(Tree.Height), y=X.CO2., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()

ggplot(data=gc_data.1,
       aes(x=as.numeric(Tree.Height), y=X.N2O., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()

ggplot(data=gc_data.1,
       aes(x=as.numeric(Tree.Height), y=X.O2., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()

##

gc_data.2<-gc_data[which(is.na(gc_data$Tree.Tissue)==F),]
gc_data.2<-gc_data.2[-which(gc_data.2$Sample.Type=="Lab Air Blank"),]
gc_data.2<-gc_data.2[-which(gc_data.2$Sample.Type=="Anaerobic (N2)"),]
gc_data.2<-gc_data.2[-which(gc_data.2$Sample.Type=="Check Standard"),]
gc_data.2<-gc_data.2[-which(is.na(gc_data.2$Tree.Height)==T),]

ggplot(data=gc_data.2,
       aes(x=as.numeric(Tree.Height), y=X.CH4., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()+ggtitle("Aerobic Incubation CH4")

ggplot(data=gc_data.2,
       aes(x=as.numeric(Tree.Height), y=X.CO2., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()

ggplot(data=gc_data.2,
       aes(x=as.numeric(Tree.Height), y=X.N2O., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()

ggplot(data=gc_data.2,
       aes(x=as.numeric(Tree.Height), y=X.O2., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()


heights<-c(rev(unique(gc_data.2$Tree.Height)[1:7]))
ch4flux<-
c(0.0232,
  0.0658,
  0.0243,
  0.2758,
  0.0302,
  0.0312,
  0.0023)
fluxdf<-as.data.frame(cbind(heights, ch4flux))


pgg1<-
ggplot(data=fluxdf,
      aes(x=as.numeric(heights), y=as.numeric(ch4flux)))+
  geom_point()+geom_smooth(se=F)+coord_flip()+ggtitle("Tree Stem Flux")

pgg2<-
ggplot(data=gc_data.1,
       aes(x=as.numeric(Tree.Height), y=X.CH4., color=Tree.Tissue))+
  geom_point()+geom_smooth(se=F)+coord_flip()+ggtitle("Anaerobic Incubation CH4 & Internal Gas")

library(ggpubr)
ggarrange(pgg1, pgg2)
