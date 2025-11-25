#Compare clusters found in heatmaps with actual clusters (genetic pedigree or social pedigree)

#code from Lewis et. al 2021, used for getting kinship matrix
library(kinship2)

#import manual note classifications
manual.note.classes=read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_UnitTable2_Tempo.csv",stringsAsFactors=FALSE,fill=TRUE)
#import spectral data on individual notes
spectral.note.data=read.csv("Lewisetal2021_AcousticParameters.csv",header=TRUE)
#import metadata on individual birds
meta.data=read.csv("~/Dropbox (The University of Manchester)/Java_Sparrow_Temporal/Lewisetal2021_metadata.csv")
meta.data$X1 <- NULL
meta.data=meta.data[order(meta.data$Bird.ID),]

#convert data to characters
meta.data$Bird.ID=as.character(meta.data$Bird.ID)
meta.data$Genetic.Father=as.character(meta.data$Genetic.Father)
meta.data$Genetic.Mother=as.character(meta.data$Genetic.Mother)
meta.data$GGF_G=as.character(meta.data$GGF_F)
meta.data$GGM_F=as.character(meta.data$GGM_F)
meta.data$GGF_M=as.character(meta.data$GGF_M)
meta.data$GGM_M=as.character(meta.data$GGM_M)

#create pedigree to use for regressions
#clean the meta.data to replace unknown parents with NAs
cols.w.ors=grep(" or ",meta.data)
for (i in 1:length(cols.w.ors)){
  meta.data[grep(" or ",meta.data[,cols.w.ors[i]]),cols.w.ors[i]]=NA
}
cols.w.ors=grep("nknown",meta.data)
for (i in 1:length(cols.w.ors)){
  meta.data[grep("nknown",meta.data[,cols.w.ors[i]]),cols.w.ors[i]]=NA
}
recover.meta.data=meta.data

#use information from the grandparents columns to extend the pedigree
male.count=sum(!is.na(unique(c(meta.data$Bird.ID,meta.data$Genetic.Father))))
ped.id=c(meta.data$Bird.ID,unique(meta.data$Genetic.Father[which(!(meta.data$Genetic.Father %in% meta.data$Bird.ID) & !is.na(meta.data$Genetic.Father))]),unique(meta.data$Genetic.Mother[which(!(meta.data$Genetic.Mother %in% meta.data$Bird.ID) & !is.na(meta.data$Genetic.Mother))]))
ped.gf=c(meta.data$Genetic.Father,rep(NA,length(ped.id)-length(meta.data$Genetic.Father)))
ped.gm=c(meta.data$Genetic.Mother,rep(NA,length(ped.id)-length(meta.data$Genetic.Mother)))
#use information about grandparents
for (i in 1:length(ped.id)){
  if (ped.id[i] %in% ped.gf){
    if (is.na(ped.gf[i])){
      ped.gf[i]=as.character(meta.data$GGF_F[match(ped.id[i],meta.data$Genetic.Father)])
    }
    if (is.na(ped.gm[i])){
      ped.gm[i]=as.character(meta.data$GGM_F[match(ped.id[i],meta.data$Genetic.Father)])
    }
  }
  if (ped.id[i] %in% ped.gm){
    if (is.na(ped.gf[i])){
      ped.gf[i]=as.character(meta.data$GGF_M[match(ped.id[i],meta.data$Genetic.Mother)])
    }
    if (is.na(ped.gm[i])){
      ped.gm[i]=as.character(meta.data$GGM_M[match(ped.id[i],meta.data$Genetic.Mother)])
    }
  }
}
#include new additions in ped.id
par.list=unique(c(ped.gf,ped.gm))
par.list=par.list[!is.na(par.list)]
new.ids=setdiff(par.list,ped.id)
new.ids=new.ids[which(new.ids!="Unknown")]
new.sex=rep(2,length(new.ids))
new.sex[which(new.ids %in% ped.gf)]=1
ped.id=c(ped.id,new.ids)
ped.gf=c(ped.gf,rep(NA,length(new.ids)))
ped.gm=c(ped.gm,rep(NA,length(new.ids)))

#create the pedigree
full.ped=pedigree(id=ped.id,dadid=ped.gf,momid=ped.gm,sex=c(rep(1,male.count),rep(2,length(ped.id)-male.count-length(new.sex)),new.sex),affected=0+(ped.id %in% manual.note.classes$song_individual),missid="Unknown")
plot.pedigree(full.ped,id=substr(ped.id,4,6),affected=full.ped$affected,cex=.7,align=c(1.5,2))
#build the kinship matrix
kin=kinship(full.ped)
#trim kinship matrix, removing mothers and grandparents we have no recordings of
ID=as.character(meta.data$Bird.ID)
kin.trim=kin[match(ID,rownames(kin)),match(ID,colnames(kin))]

#Lewis et all code ends here

#from alignment script

k <- heatmap.2(t(note.prop.mat), 
               col=terrain.colors(256), trace="none",
               main="Proportion of notes recorded (per bird)",
               key.title="Key", key.xlab="Proportion of notes",
               margins=c(4,7))

plot(k$colDendrogram)


M <- 2*kin.trim 
#-1 because no relationship is considered further apart
as.dist(1-M)
gen_tree <- hclust(as.dist(1-M))
gen_tree<- as.dendrogram(gen_tree)

library(dendextend)
comp_trees = dendlist(gen_tree, k$colDendrogram)
cor.dendlist(comp_trees)
#cophenetic correlation of the trees is ~0.2

hclust(kin.trim)
as.dist(kin.trim)
?hclust()

