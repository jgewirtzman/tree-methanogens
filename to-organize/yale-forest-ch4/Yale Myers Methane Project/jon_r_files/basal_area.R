library(tidyverse)
trees<-read.csv("/Users/jongewirtzman/Dropbox (Yale_FES)/Yale Myers Methane Project/Environmental variables/ForestGEO_Data_Final_2019.csv")

str(trees)
unique(trees$Species_Code)

trees$DBH<-as.numeric(trees$DBH)
trees<-trees[-which(trees$DBH>100), ]
trees<-trees[-which(trees$Species_Code==""), ]

sort(trees$DBH, decreasing=T)
plot(trees$DBH~trees$Species_Code)


trees$basal_area<-pi*(trees$DBH/2)^2


basal_area<-
  trees %>% 
  group_by(Species_Code) %>% 
  summarise(basal = sum(basal_area, na.rm=T))

basal_area<-basal_area[order(-basal_area$basal),]

basal_area %>% 
  ggplot(aes(reorder(Species_Code,
                     -basal),
             basal))+
  geom_col()+
  labs(x="Species", y="Basal Area",
       title="YMF Forest-GEO Basal Area by Species")



species.target<-as.matrix(basal_area[c(1:8, 12),1])


par(mfrow=c(3,3))
for(i in 1:length(species.target)){
  hist(trees$DBH[which(trees$Species_Code==species.target[i] & trees$DBH>20)],
       main=species.target[i],
       xlab="DBH (Stems >20cm)")
  max.dbh<-max(trees$DBH[which(trees$Species_Code==species.target[i] & trees$DBH>20)], na.rm=T)
  min.dbh<-min(trees$DBH[which(trees$Species_Code==species.target[i] & trees$DBH>20)], na.rm=T)
  range<-max.dbh-min.dbh
  abline(v=1/3*range+min.dbh, col="red", lty=2)
  abline(v=2/3*range+min.dbh, col="red", lty=2)
}

trees_sample<-trees
trees_sample<-trees[which(trees$Species_Code %in% species.target), ]
trees_sample<-trees[which(trees$DBH>20), ]

par(mfrow=c(1,1))
plot(trees_sample$PX, trees_sample$PY)


max<-c()
min<-c()
for(i in 1:length(species.target)){
  print(nrow(trees_sample[which(trees_sample$Species_Code==as.character(species.target[i])),]))
  max[i]<-max(trees_sample$DBH[which(trees_sample$Species_Code==as.character(species.target[i]))], na.rm=T)
  min[i]<-min(trees_sample$DBH[which(trees_sample$Species_Code==as.character(species.target[i]))], na.rm=T)
}
tree_crit<-data.frame(species.target, max, min)
tree_crit$range<-tree_crit$max-tree_crit$min
tree_crit$q1<-tree_crit$min+1/4*tree_crit$range
tree_crit$q2<-tree_crit$min+2/4*tree_crit$range
tree_crit$q3<-tree_crit$min+3/4*tree_crit$range


selected<-c()

for(i in 1:length(tree_crit$Species_Code)){

  q1tf<-  which(trees_sample$Species_Code==as.character(tree_crit[i,1])  &
                  trees_sample$DBH>=as.numeric(tree_crit[i, 'min']) &
                  trees_sample$DBH<=as.numeric(tree_crit[i, 'q1']) )
  
  q2tf<-  which(trees_sample$Species_Code==as.character(tree_crit[i,1])  &
                  trees_sample$DBH>=as.numeric(tree_crit[i, 'q1']) &
                  trees_sample$DBH<=as.numeric(tree_crit[i, 'q2']) )
  
  q3tf<-  which(trees_sample$Species_Code==as.character(tree_crit[i,1])  &
                  trees_sample$DBH>=as.numeric(tree_crit[i, 'q2']) &
                  trees_sample$DBH<=as.numeric(tree_crit[i, 'q3']) )
  
  q4tf<-  which(trees_sample$Species_Code==as.character(tree_crit[i,1])  &
                  trees_sample$DBH>=as.numeric(tree_crit[i, 'q3']) &
                  trees_sample$DBH<=as.numeric(tree_crit[i, 'max']) )
  
  selected<-c(selected, 
              if(length(q1tf)==1){q1tf}else{sample(q1tf, 
                     size=min(c(8, length(q1tf)))
              )},
              
              if(length(q2tf)==1){q2tf}else{sample(q2tf, 
                     min(c(7, length(q2tf)))
              )},
              
              if(length(q3tf)==1){q3tf}else{sample(q3tf, 
                     min(c(7, length(q3tf)))
              )},
              
              if(length(q4tf)==1){q4tf}else{sample(q4tf, 
                     min(c(8, length(q4tf)))
              )}
              )
  
}

trees_sample<-trees_sample[selected,]

ggplot(data=trees_sample, aes(x=PX, y=PY, color=Species_Code, size=DBH))+
  geom_point()+
  xlim(0,200)+ylim(0,200)

ggplot(data=trees_sample, aes(x=Species_Code, y=DBH))+
  geom_violin()+
  geom_point()
  xlim(0,200)+ylim(0,200)
  
  write.csv(trees_sample, "trees_to_sample.csv")
  
  
  
  trees2<-read.csv("/Users/jongewirtzman/Google Drive/Methane Research/Trees_DBH.csv")
  ggplot(data=trees2, aes(x=Species, y=as.numeric(as.character(DBH))))+
    geom_violin()+
    geom_point()
  xlim(0,200)+ylim(0,200)