# --------------------------------------------------------
# --------------------------------------------------------
#    Unraveling the mortality, growth and recruitment processes 
#    underlying the population dynamics of a small pelagic system 
#    using a Bayesian life cycle model 
#
# Alice Bordes, M2 Internship, Febrary - July 2023  
# A two-species (sardine and anchovy) life-cycle model - Golfe de Gascogne - 2000 à 2022
#
# --------------------------------------------------------
# --------------------------------------------------------

# Données biologique (Abondance et longueur corporelle) proviennent des campagnes PELGAS
# Données de capture proviennent du groupe de travail du CIEM
# Données "environnementales" proviennent du modèle SEAPODYM
# cas d'étude : anchois 

rm(list=ls())

#importation librairie
library(dclone)
library(rjags)
library(MASS)
library(coda)
library(reshape2)
library(tidyverse)
library(rlang)
load.module("dic") # module permettant de stocker les paramètres (les déviances et les termes qui permettent de calculer la pénalisation par la complexité du modèle) pour le calcul du WAIC lors de l'exécution du modèle


##nb_annees
n_annee=23

####################################

#anhovy
#complet sur M
covar.set.M1.anch<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")
covar.set.M2.anch<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")
covar.set.M3.anch<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")

#complet sur g
covar.set.g1.anch<-c("nsprat","AMO","NAO","Temp.jul_oct","Zoo.jul_oct")
covar.set.g2.anch<-c("nsprat","AMO","NAO","Temp.jul_oct","Zoo.jul_oct") 
covar.set.g3.anch<-c("nsprat","AMO","NAO","Temp.jul_oct","Zoo.jul_oct")

#complet sur Z.SSN
covar.set.ZSSN.anch<-c("AMO","NAO","Temp.jul_oct","Zoo.jul_oct","Temp.nov_feb","Zoo.nov_feb","Temp.apr_jun","Zoo.apr_jun")

#sardine
#complet sur M1
covar.set.M1.sard<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")
covar.set.M2.sard<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")
covar.set.M3.sard<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")
covar.set.M4.sard<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")
covar.set.M5.sard<-c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb")

#complet sur g
covar.set.g1.sard<-c("nsprat","AMO","NAO","Temp.jun_oct","Zoo.jun_oct")
covar.set.g2.sard<-c("nsprat","AMO","NAO","Temp.jun_oct","Zoo.jun_oct")
covar.set.g3.sard<-c("nsprat","AMO","NAO","Temp.jun_oct","Zoo.jun_oct")
covar.set.g4.sard<-c("nsprat","AMO","NAO","Temp.jun_oct","Zoo.jun_oct")
covar.set.g5.sard<-c("nsprat","AMO","NAO","Temp.jun_oct","Zoo.jun_oct")

#complet sur Z.SSN
covar.set.ZSSN.sard<-c("AMO","NAO","Temp.jun_oct","Zoo.jun_oct","Temp.nov_feb","Zoo.nov_feb","Temp.mar_may","Zoo.mar_may","Temp.jun_oct_ante","Zoo.jun_oct_ante","Temp.nov_feb_ante","Zoo.nov_feb_ante")

####################################################

# --------------------------------------------------------
#              Importation données 
# --------------------------------------------------------

#par convention : 1 = anchois, 2 = sardine

##loading data
base<-("C:/Users/alice/Documents/stage_M2_bayesien/work2023/Stage_M2_Alice_PELGAS/Github/stageM2_2023_spf/")
setwd(paste(base,"02_TRANSFORMED_DATA",sep=""))
data_anchovy <- read.table("data_anchovy.txt", sep="\t", h=T) # Données biologiques
data_sardine <- read.table("data_sardine.txt", sep="\t", h=T) # Données biologiques

##matrice de maturite
maturite_anch<-matrix(nrow=23,ncol=6,rep(1))
maturite_anch[,5]<-rep(NA)
maturite_anch[,6]<-rep(NA)
maturite.sard_obs<-read.csv2(paste(base,"01_RAW_DATA/Age_At_Mat_obs_sardines.csv",sep=""), sep = ";", h=T,dec=",")
maturite.sard_obs<-maturite.sard_obs[,-c(1:2)]
maturite_sard<-matrix(nrow=23,ncol=6,rep(1))
for(i in 1:nrow(maturite.sard_obs)){
  for(j in 1:ncol(maturite.sard_obs)){
    maturite_sard[i,j]<-maturite.sard_obs[i,j]
  }
}

mat.list=as.vector(c(maturite_anch,maturite_sard))
mat=array(mat.list,dim=c(23,6,2),dimnames=list(1:23,c("age1","age2","age3","age4","age5","age6plus"),c("1","2")))


##abondance sprats
setwd(paste(base,"01_RAW_DATA",sep=""))
data_sprat<- read.table("biomCVallSpPELGASseriesLong.csv", sep=";",dec=",", h=T) # Données bio
data_sprat<-data_sprat %>% filter(sp=="SPRA-SPR")  %>% group_by(year,wabun) %>% summarise(year,wabun)
data_sprat<-as.data.frame(data_sprat)
data_sprat<-rbind(data_sprat,c(2020,mean(data_sprat[16:20,2])))
data_sprat<-data_sprat%>%arrange(year)

### Densité ###
n1.obs.anch.cr <- (data_anchovy$N1-mean(data_anchovy$N1[1:n_annee]))/sd(data_anchovy$N1[1:n_annee]) # Abondance observee à l'âge 1
n2.obs.anch.cr <- (data_anchovy$N2-mean(data_anchovy$N2[1:n_annee]))/sd(data_anchovy$N2[1:n_annee])
n3.obs.anch.cr <- (data_anchovy$N3-mean(data_anchovy$N3[1:n_annee]))/sd(data_anchovy$N3[1:n_annee])
n4plus.obs.anch.cr <- (data_anchovy$N4plus-mean(data_anchovy$N4plus[1:n_annee]))/sd(data_anchovy$N4plus[1:n_annee])
nlarge.obs.anch.cr <- ((data_anchovy$N3+data_anchovy$N4plus)-mean((data_anchovy$N3+data_anchovy$N4plus)[1:n_annee]))/sd((data_anchovy$N3+data_anchovy$N4plus)[1:n_annee])
n.ssn.obs.anch.cr<-  ((data_anchovy$N1*mat[,1,1]+data_anchovy$N2*mat[,2,1]+data_anchovy$N3*mat[,3,1]+data_anchovy$N4plus*mat[,4,1])-mean((data_anchovy$N1*mat[,1,1]+data_anchovy$N2*mat[,2,1]+data_anchovy$N3*mat[,3,1]+data_anchovy$N4plus*mat[,4,1]))[1:n_annee][1])/sd((data_anchovy$N1*mat[,1,1]+data_anchovy$N2*mat[,2,1]+data_anchovy$N3*mat[,3,1]+data_anchovy$N4plus*mat[,4,1]))[1:n_annee][1]
n.ssn.obs.anch<-  (data_anchovy$N1*mat[,1,1]+data_anchovy$N2*mat[,2,1]+data_anchovy$N3*mat[,3,1]+data_anchovy$N4plus*mat[,4,1])

n1.obs.sard.cr <- (data_sardine$N1-mean(data_sardine$N1[1:n_annee]))/sd(data_sardine$N1[1:n_annee]) # Abondance observee à l'âge 1
n2.obs.sard.cr <- (data_sardine$N2-mean(data_sardine$N2[1:n_annee]))/sd(data_sardine$N2[1:n_annee])
n3.obs.sard.cr <- (data_sardine$N3-mean(data_sardine$N3[1:n_annee]))/sd(data_sardine$N3[1:n_annee])
n4.obs.sard.cr <- (data_sardine$N4-mean(data_sardine$N4[1:n_annee]))/sd(data_sardine$N4[1:n_annee])
n5.obs.sard.cr <- (data_sardine$N5-mean(data_sardine$N5[1:n_annee]))/sd(data_sardine$N5[1:n_annee])
n6plus.obs.sard.cr <- (data_sardine$N6plus-mean(data_sardine$N6plus[1:n_annee]))/sd(data_sardine$N6plus[1:n_annee])
nlarge.obs.sard.cr <- ((data_sardine$N3+data_sardine$N4+data_sardine$N5+data_sardine$N6plus)-mean((data_sardine$N3+data_sardine$N4+data_sardine$N5+data_sardine$N6plus)[1:n_annee]))/sd((data_sardine$N3+data_sardine$N4+data_sardine$N5+data_sardine$N6plus)[1:n_annee])
n.ssn.obs.sard.cr<-  ((data_sardine$N1*mat[,1,2]+data_sardine$N2*mat[,2,2]+data_sardine$N3*mat[,3,2]+data_sardine$N4*mat[,4,2]+data_sardine$N5*mat[,5,2]+data_sardine$N6plus*mat[,6,2])-mean((data_sardine$N1*mat[,1,2]+data_sardine$N2*mat[,2,2]+data_sardine$N3*mat[,3,2]+data_sardine$N4*mat[,4,2]+data_sardine$N5*mat[,5,2]+data_sardine$N6plus*mat[,6,2]))[1:n_annee][1])/sd((data_sardine$N1*mat[,1,2]+data_sardine$N2*mat[,2,2]+data_sardine$N3*mat[,3,2]+data_sardine$N4*mat[,4,2]+data_sardine$N5*mat[,5,2]+data_sardine$N6plus*mat[,6,2]))[1:n_annee][1]
n.ssn.obs.sard<-(data_sardine$N1*mat[,1,2]+data_sardine$N2*mat[,2,2]+data_sardine$N3*mat[,3,2]+data_sardine$N4*mat[,4,2]+data_sardine$N5*mat[,5,2]+data_sardine$N6plus*mat[,6,2])

nsprat.obs.cr <- (data_sprat$wabun-mean(data_sprat$wabun[1:n_annee]))/sd(data_sprat$wabun[1:n_annee])

l1.obs.anch.cr <- (data_anchovy$L1-mean(data_anchovy$L1[1:n_annee]))/sd(data_anchovy$L1[1:n_annee]) # Abondance observee à l'âge 1
l2.obs.anch.cr <- (data_anchovy$L2-mean(data_anchovy$L2[1:n_annee]))/sd(data_anchovy$L2[1:n_annee])
l3.obs.anch.cr <- (data_anchovy$L3-mean(data_anchovy$L3[1:n_annee]))/sd(data_anchovy$L3[1:n_annee])
l4.obs.anch.cr <- (data_anchovy$L4-mean(data_anchovy$L4[1:n_annee]))/sd(data_anchovy$L4[1:n_annee])
#l.ssn.obs.anch.cr<-(((data_anchovy$L1*(data_anchovy$N1*mat[,1,1]/data_anchovy$N_total)+
#                     data_anchovy$L2*(data_anchovy$N2*mat[,2,1]/data_anchovy$N_total)+
#                     data_anchovy$L3*(data_anchovy$N3*mat[,3,1]/data_anchovy$N_total)+
#                     data_anchovy$L4*(data_anchovy$N4plus*mat[,4,1]/data_anchovy$N_total))
#                    -mean(  data_anchovy$L1[1:n_annee]*(data_anchovy$N1[1:n_annee]*mat[,1,1]/data_anchovy$N_total[1:n_annee])+
#                            data_anchovy$L2[1:n_annee]*(data_anchovy$N2[1:n_annee]*mat[,2,1]/data_anchovy$N_total[1:n_annee])+
#                            data_anchovy$L3[1:n_annee]*(data_anchovy$N3[1:n_annee]*mat[,3,1]/data_anchovy$N_total[1:n_annee])+
#                            data_anchovy$L4[1:n_annee]*(data_anchovy$N4plus[1:n_annee]*mat[,4,1]/data_anchovy$N_total[1:n_annee])))
#                    /sd( data_anchovy$L1[1:n_annee]*(data_anchovy$N1[1:n_annee]*mat[,1,1]/data_anchovy$N_total[1:n_annee])+
#                         data_anchovy$L2[1:n_annee]*(data_anchovy$N2[1:n_annee]*mat[,2,1]/data_anchovy$N_total[1:n_annee])+
#                         data_anchovy$L3[1:n_annee]*(data_anchovy$N3[1:n_annee]*mat[,3,1]/data_anchovy$N_total[1:n_annee])+
#                         data_anchovy$L4[1:n_annee]*(data_anchovy$N4plus[1:n_annee]*mat[,4,1]/data_anchovy$N_total[1:n_annee])))
l.ssn.obs.anch.cr<-((((data_anchovy$L1*(data_anchovy$N1*mat[,1,1])+
                        data_anchovy$L2*(data_anchovy$N2*mat[,2,1])+
                        data_anchovy$L3*(data_anchovy$N3*mat[,3,1])+
                        data_anchovy$L4*(data_anchovy$N4plus*mat[,4,1]))/n.ssn.obs.anch)
                     -mean(  (data_anchovy$L1[1:n_annee]*(data_anchovy$N1[1:n_annee]*mat[,1,1])+
                               data_anchovy$L2[1:n_annee]*(data_anchovy$N2[1:n_annee]*mat[,2,1])+
                               data_anchovy$L3[1:n_annee]*(data_anchovy$N3[1:n_annee]*mat[,3,1])+
                               data_anchovy$L4[1:n_annee]*(data_anchovy$N4plus[1:n_annee]*mat[,4,1]))/n.ssn.obs.anch[1:n_annee]))
                    /sd( (data_anchovy$L1[1:n_annee]*(data_anchovy$N1[1:n_annee]*mat[,1,1])+
                           data_anchovy$L2[1:n_annee]*(data_anchovy$N2[1:n_annee]*mat[,2,1])+
                           data_anchovy$L3[1:n_annee]*(data_anchovy$N3[1:n_annee]*mat[,3,1])+
                           data_anchovy$L4[1:n_annee]*(data_anchovy$N4plus[1:n_annee]*mat[,4,1]))/n.ssn.obs.anch[1:n_annee]))

l1.obs.sard.cr <- (data_sardine$L1-mean(data_sardine$L1[1:n_annee]))/sd(data_sardine$L1[1:n_annee]) # Abondance observee à l'âge 1
l2.obs.sard.cr <- (data_sardine$L2-mean(data_sardine$L2[1:n_annee]))/sd(data_sardine$L2[1:n_annee])
l3.obs.sard.cr <- (data_sardine$L3-mean(data_sardine$L3[1:n_annee]))/sd(data_sardine$L3[1:n_annee])
l4.obs.sard.cr <- (data_sardine$L4-mean(data_sardine$L4[1:n_annee]))/sd(data_sardine$L4[1:n_annee])
l5.obs.sard.cr <- (data_sardine$L5-mean(data_sardine$L5[1:n_annee]))/sd(data_sardine$L5[1:n_annee])
l6.obs.sard.cr <- (data_sardine$L6-mean(data_sardine$L6[1:n_annee]))/sd(data_sardine$L6[1:n_annee])
#l.ssn.obs.sard.cr<-(((  data_sardine$L1*(data_sardine$N1*mat[,1,2]/data_sardine$N_total)+
#                        data_sardine$L2*(data_sardine$N2*mat[,2,2]/data_sardine$N_total)+
#                       data_sardine$L3*(data_sardine$N3*mat[,3,2]/data_sardine$N_total)+
#                       data_sardine$L4*(data_sardine$N4*mat[,4,2]/data_sardine$N_total)+
#                       data_sardine$L5*(data_sardine$N5*mat[,5,2]/data_sardine$N_total)+
#                       data_sardine$L6*(data_sardine$N6plus*mat[,6,2]/data_sardine$N_total))
#                    -mean(  data_sardine$L1[1:n_annee]*(data_sardine$N1[1:n_annee]*mat[,1,2]/data_sardine$N_total[1:n_annee])+
#                            data_sardine$L2[1:n_annee]*(data_sardine$N2[1:n_annee]*mat[,2,2]/data_sardine$N_total[1:n_annee])+
#                            data_sardine$L3[1:n_annee]*(data_sardine$N3[1:n_annee]*mat[,3,2]/data_sardine$N_total[1:n_annee])+
#                            data_sardine$L4[1:n_annee]*(data_sardine$N4[1:n_annee]*mat[,4,2]/data_sardine$N_total[1:n_annee])+
#                            data_sardine$L5[1:n_annee]*(data_sardine$N5[1:n_annee]*mat[,5,2]/data_sardine$N_total[1:n_annee])+
#                            data_sardine$L6[1:n_annee]*(data_sardine$N6plus[1:n_annee]*mat[,6,2]/data_sardine$N_total[1:n_annee])))
#                   /sd( data_sardine$L1[1:n_annee]*(data_sardine$N1[1:n_annee]*mat[,1,2]/data_sardine$N_total[1:n_annee])+
#                        data_sardine$L2[1:n_annee]*(data_sardine$N2[1:n_annee]*mat[,2,2]/data_sardine$N_total[1:n_annee])+
#                        data_sardine$L3[1:n_annee]*(data_sardine$N3[1:n_annee]*mat[,3,2]/data_sardine$N_total[1:n_annee])+
#                        data_sardine$L4[1:n_annee]*(data_sardine$N4[1:n_annee]*mat[,4,2]/data_sardine$N_total[1:n_annee])+
#                        data_sardine$L5[1:n_annee]*(data_sardine$N5[1:n_annee]*mat[,5,2]/data_sardine$N_total[1:n_annee])+
#                        data_sardine$L6[1:n_annee]*(data_sardine$N6plus[1:n_annee]*mat[,6,2]/data_sardine$N_total[1:n_annee])))
l.ssn.obs.sard.cr<-((((data_sardine$L1*(data_sardine$N1*mat[,1,2])+
                         data_sardine$L2*(data_sardine$N2*mat[,2,2])+
                         data_sardine$L3*(data_sardine$N3*mat[,3,2])+
                         data_sardine$L4*(data_sardine$N4*mat[,4,2])+
                         data_sardine$L5*(data_sardine$N5*mat[,5,2])+
                         data_sardine$L6*(data_sardine$N6plus*mat[,6,2]))/n.ssn.obs.sard)
                     -mean(  (data_sardine$L1[1:n_annee]*(data_sardine$N1[1:n_annee]*mat[,1,2])+
                                data_sardine$L2[1:n_annee]*(data_sardine$N2[1:n_annee]*mat[,2,2])+
                                data_sardine$L3[1:n_annee]*(data_sardine$N3[1:n_annee]*mat[,3,2])+
                                data_sardine$L4[1:n_annee]*(data_sardine$N4[1:n_annee]*mat[,4,2])+
                                data_sardine$L5[1:n_annee]*(data_sardine$N5[1:n_annee]*mat[,5,2])+
                                data_sardine$L6[1:n_annee]*(data_sardine$N6plus[1:n_annee]*mat[,6,2]))/n.ssn.obs.sard[1:n_annee]))
                    /sd( (data_sardine$L1[1:n_annee]*(data_sardine$N1[1:n_annee]*mat[,1,2])+
                            data_sardine$L2[1:n_annee]*(data_sardine$N2[1:n_annee]*mat[,2,2])+
                            data_sardine$L3[1:n_annee]*(data_sardine$N3[1:n_annee]*mat[,3,2])+
                            data_sardine$L4[1:n_annee]*(data_sardine$N4[1:n_annee]*mat[,4,2])+
                            data_sardine$L5[1:n_annee]*(data_sardine$N5[1:n_annee]*mat[,5,2])+
                            data_sardine$L6[1:n_annee]*(data_sardine$N6plus[1:n_annee]*mat[,6,2]))/n.ssn.obs.sard[1:n_annee]))



## AMO and NAO
setwd(paste(base,"01_RAW_DATA",sep=""))
data_AMO <- read.table("amon.us.data.txt", sep="",dec=".",fill = TRUE) # Atlantic Multi-decadal Oscillation (AMO) Index
names(data_AMO)<-c("Year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
data_AMO<-data_AMO[-c(77:81),]
data_AMO<-data_AMO %>% filter(Year>=1999)
data_AMO$Year<-as.numeric(data_AMO$Year)
data_AMO$Jan<-as.numeric(data_AMO$Jan)
data_AMO$Feb<-as.numeric(data_AMO$Feb)
data_AMO$Mar<-as.numeric(data_AMO$Mar)
data_AMO$Apr<-as.numeric(data_AMO$Apr)
data_AMO$May<-as.numeric(data_AMO$May)
data_AMO$Jun<-as.numeric(data_AMO$Jun)
data_AMO$Jul<-as.numeric(data_AMO$Jul)
#data_AMO$Annuel<-rowMeans(data_AMO[,2:13])

AMO<-data.frame("Year"=2000:2022,"AMO_index"=rep(NA))
for(i in 1:23){
  selecA<-data_AMO[data_AMO$Year== 1999+i,c("May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")]
  selecB<-data_AMO[data_AMO$Year== 2000+i,c("Jan","Feb","Mar","Apr")]
  vect<-c(selecA,selecB)
  vect<-unlist(vect)
  AMO$AMO_index[i]<-mean(vect)
}


##NAO
# = anomalies de T°C basé sur le calcul de la moy glissante sur 10 années, soustrait l'effet du CC considéré comme un effet linéaire
data_NAO <- read.table("nao_monthly.txt", sep="",dec=".",fill = TRUE) # Hurrell North Atlantic Oscillation (NAO) Index (station-based)
data_NAO<-cbind("Year"=rownames(data_NAO),data_NAO)
data_NAO<-data_NAO %>% filter(Year>=1999)
data_NAO$Year<-as.numeric(data_NAO$Year)
#data_NAO$Annuel<-rowMeans(data_NAO[,2:13])
#data_NAO<-data_NAO[-nrow(data_NAO),]
NAO<-data.frame("Year"=2000:(max(data_NAO$Year)-1),"NAO_index"=rep(NA))
for(i in 1:22){
  selecA<-data_NAO[data_NAO$Year== 1999+i,c("May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")]
  selecB<-data_NAO[data_NAO$Year== 2000+i,c("Jan","Feb","Mar","Apr")]
  vect<-c(selecA,selecB)
  vect<-unlist(vect) 
  NAO$NAO_index[i]<-mean(vect)
  
  #pour t=1, on aura les mois de mai à dec de l'an 2000 + janv à avril de 2001 = conditions rencontrés par les N1 de l'an 1 (2000) durant l'an qui les aura menés à N2
}
NAO[23,2]<-mean(NAO[18:22,2])


##env data
##env anchois
setwd(paste(base,"02_TRANSFORMED_DATA",sep=""))
env_anchovy_stand<- read.table("Env_stand_anchovy.txt", sep="\t", h=T) # Données env
dtplus<-matrix(ncol=ncol(env_anchovy_stand),nrow=2)
dtplus[,1]<-c(2021,2022)
colnames(dtplus)<-colnames(env_anchovy_stand)
env_anchovy_stand<-rbind(env_anchovy_stand,dtplus)
for(i in 22:nrow(env_anchovy_stand))
{
  for(j in 2:ncol(env_anchovy_stand))
  {
    env_anchovy_stand[i,j]<-mean(env_anchovy_stand[17:21,j]) 
  }
}

#env sardine
env_sardine_stand<- read.table("Env_stand_sardine.txt", sep="\t", h=T) # Données env
dtplus<-matrix(ncol=ncol(env_sardine_stand),nrow=2)
dtplus[,1]<-c(2021,2022)
colnames(dtplus)<-colnames(env_sardine_stand)
env_sardine_stand<-rbind(env_sardine_stand,dtplus)
for(i in 22:nrow(env_sardine_stand))
{
  for(j in 2:ncol(env_sardine_stand))
  {
    env_sardine_stand[i,j]<-mean(env_sardine_stand[17:21,j]) 
  }
}



##------- covar -------##

data_covar.all<-data.frame( "Year"=2000:(1999+(n_annee-1)),
                            "n1.anch"=n1.obs.anch.cr[1:(n_annee-1)],"n2.anch"=n2.obs.anch.cr[1:(n_annee-1)],"n3.anch"=n3.obs.anch.cr[1:(n_annee-1)],"n4plus.anch"=n4plus.obs.anch.cr[1:(n_annee-1)],
                            "n1.sard"=n1.obs.sard.cr[1:(n_annee-1)],"n2.sard"=n2.obs.sard.cr[1:(n_annee-1)],"n3.sard"=n3.obs.sard.cr[1:(n_annee-1)],"n4.sard"=n4.obs.sard.cr[1:(n_annee-1)],"n5.sard"=n5.obs.sard.cr[1:(n_annee-1)],"n6.sard"=n6plus.obs.sard.cr[1:(n_annee-1)],
                            "nlarge.anch"=nlarge.obs.anch.cr[1:(n_annee-1)],"nlarge.sard"=nlarge.obs.sard.cr[1:(n_annee-1)],"nsprat"=nsprat.obs.cr[1:(n_annee-1)],
                            "l1.anch"=l1.obs.anch.cr[1:(n_annee-1)],"l2.anch"=l2.obs.anch.cr[1:(n_annee-1)],"l3.anch"=l3.obs.anch.cr[1:(n_annee-1)],"l4.anch"=l4.obs.anch.cr[1:(n_annee-1)],
                            "l1.sard"=l1.obs.sard.cr[1:(n_annee-1)],"l2.sard"=l2.obs.sard.cr[1:(n_annee-1)],"l3.sard"=l3.obs.sard.cr[1:(n_annee-1)],"l4.sard"=l4.obs.sard.cr[1:(n_annee-1)],"l5.sard"=l5.obs.sard.cr[1:(n_annee-1)],"l6.sard"=l6.obs.sard.cr[1:(n_annee-1)],
                            "NAO"=NAO$NAO_index[1:(n_annee-1)],"AMO"=AMO$AMO_index[1:(n_annee-1)],
                            "Zoo.nov_feb"=env_anchovy_stand[1:((n_annee-1)),2],"Temp.nov_feb"=env_anchovy_stand[1:((n_annee-1)),3],"Zoo.jul_oct"=env_anchovy_stand[1:(n_annee-1),4],"Temp.jul_oct"=env_anchovy_stand[1:(n_annee-1),5],
                            "Zoo.jun_oct"=env_sardine_stand[1:(n_annee-1),4],"Temp.jun_oct"=env_sardine_stand[1:(n_annee-1),5],
                            "Zoo.apr_jun"=env_anchovy_stand[1:((n_annee-1)),6],"Temp.apr_jun"=env_anchovy_stand[1:((n_annee-1)),7],
                            "Zoo.mar_may"=env_sardine_stand[1:((n_annee-1)),6],"Temp.mar_may"=env_sardine_stand[1:((n_annee-1)),7],
                            "Zoo.nov_feb_ante"=env_sardine_stand[1:((n_annee-1)),8],"Temp.nov_feb_ante"=env_sardine_stand[1:((n_annee-1)),9],
                            "Zoo.jun_oct_ante"=env_sardine_stand[1:((n_annee-1)),10],"Temp.jun_oct_ante"=env_sardine_stand[1:((n_annee-1)),11],
                            "n.ssn.anch"=n.ssn.obs.anch.cr[1:(n_annee-1)],"n.ssn.sard"=n.ssn.obs.sard.cr[1:(n_annee-1)],
                            "l.ssn.anch"=l.ssn.obs.anch.cr[1:(n_annee-1)],"l.ssn.sard"=l.ssn.obs.sard.cr[1:(n_annee-1)])

#anchovy
data_covar.M1.anch <- data_covar.all %>% select(covar.set.M1.anch)
data_covar.M2.anch <- data_covar.all %>% select(covar.set.M2.anch)
data_covar.M3.anch <- data_covar.all %>% select(covar.set.M3.anch)

data_covar.g1.anch <- data_covar.all %>% select(covar.set.g1.anch)
data_covar.g2.anch <- data_covar.all %>% select(covar.set.g2.anch)
data_covar.g3.anch <- data_covar.all %>% select(covar.set.g3.anch)

data_covar.Z.SSN.anch <- data_covar.all %>% select(covar.set.ZSSN.anch)


ncov.M.anch <- ncol(data_covar.M1.anch)   #nb de covar : NAO, n1.obs... = ncol()
ncov.g.anch <- ncol(data_covar.g1.anch)
ncov.zssn.anch <- ncol(data_covar.Z.SSN.anch)


X.M1.anch=as.matrix(data_covar.M1.anch[,1:ncol(data_covar.M1.anch)])
X.M2.anch=as.matrix(data_covar.M2.anch[,1:ncol(data_covar.M2.anch)])
X.M3.anch=as.matrix(data_covar.M3.anch[,1:ncol(data_covar.M3.anch)])

X.g1.anch=as.matrix(data_covar.g1.anch[,1:ncol(data_covar.g1.anch)])
X.g2.anch=as.matrix(data_covar.g2.anch[,1:ncol(data_covar.g2.anch)])
X.g3.anch=as.matrix(data_covar.g3.anch[,1:ncol(data_covar.g3.anch)])

X.ZSSN.anch=as.matrix(data_covar.Z.SSN.anch[,1:ncol(data_covar.Z.SSN.anch)])

#sardine
data_covar.M1.sard <- data_covar.all %>% select(covar.set.M1.sard)
data_covar.M2.sard <- data_covar.all %>% select(covar.set.M2.sard)
data_covar.M3.sard <- data_covar.all %>% select(covar.set.M3.sard)
data_covar.M4.sard <- data_covar.all %>% select(covar.set.M4.sard)
data_covar.M5.sard <- data_covar.all %>% select(covar.set.M5.sard)

data_covar.g1.sard <- data_covar.all %>% select(covar.set.g1.sard)
data_covar.g2.sard <- data_covar.all %>% select(covar.set.g2.sard)
data_covar.g3.sard <- data_covar.all %>% select(covar.set.g3.sard)
data_covar.g4.sard <- data_covar.all %>% select(covar.set.g4.sard)
data_covar.g5.sard <- data_covar.all %>% select(covar.set.g5.sard)

data_covar.Z.SSN.sard <- data_covar.all %>% select(covar.set.ZSSN.sard)


ncov.M.sard <- ncol(data_covar.M1.sard)   #nb de covar : NAO, n1.obs... = ncol()
ncov.g.sard <- ncol(data_covar.g1.sard)
ncov.zssn.sard <- ncol(data_covar.Z.SSN.sard)


X.M1.sard=as.matrix(data_covar.M1.sard[,1:ncol(data_covar.M1.sard)])
X.M2.sard=as.matrix(data_covar.M2.sard[,1:ncol(data_covar.M2.sard)])
X.M3.sard=as.matrix(data_covar.M3.sard[,1:ncol(data_covar.M3.sard)])
X.M4.sard=as.matrix(data_covar.M4.sard[,1:ncol(data_covar.M4.sard)])
X.M5.sard=as.matrix(data_covar.M5.sard[,1:ncol(data_covar.M5.sard)])

X.g1.sard=as.matrix(data_covar.g1.sard[,1:ncol(data_covar.g1.sard)])
X.g2.sard=as.matrix(data_covar.g2.sard[,1:ncol(data_covar.g2.sard)])
X.g3.sard=as.matrix(data_covar.g3.sard[,1:ncol(data_covar.g3.sard)])
X.g4.sard=as.matrix(data_covar.g4.sard[,1:ncol(data_covar.g4.sard)])
X.g5.sard=as.matrix(data_covar.g5.sard[,1:ncol(data_covar.g5.sard)])

X.ZSSN.sard=as.matrix(data_covar.Z.SSN.sard[,1:ncol(data_covar.Z.SSN.sard)])

#M
X.M1.list=as.vector(c(X.M1.anch,X.M1.sard))
X.M1=array(X.M1.list,dim=c((n_annee-1),ncov.M.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb"),c("1","2")))

X.M2.list=as.vector(c(X.M2.anch,X.M2.sard))
X.M2=array(X.M2.list,dim=c((n_annee-1),ncov.M.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb"),c("1","2")))

X.M3.list=as.vector(c(X.M3.anch,X.M3.sard))
X.M3=array(X.M3.list,dim=c((n_annee-1),ncov.M.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb"),c("1","2")))

X.M4.list=as.vector(c(X.M1.anch,X.M4.sard))
X.M4=array(X.M4.list,dim=c((n_annee-1),ncov.M.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb"),c("1","2")))

X.M5.list=as.vector(c(X.M1.anch,X.M5.sard))
X.M5=array(X.M5.list,dim=c((n_annee-1),ncov.M.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AMO","NAO","Temp.nov_feb","Zoo.nov_feb"),c("1","2")))

#g
X.g1.list=as.vector(c(X.g1.anch,X.g1.sard))
X.g1=array(X.g1.list,dim=c((n_annee-1),ncov.g.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AgO","NAO","Temp.jun/jul_oct","Zoo.jun/jul_oct"),c("1","2")))

X.g2.list=as.vector(c(X.g2.anch,X.g2.sard))
X.g2=array(X.g2.list,dim=c((n_annee-1),ncov.g.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AgO","NAO","Temp.jun/jul_oct","Zoo.jun/jul_oct"),c("1","2")))

X.g3.list=as.vector(c(X.g3.anch,X.g3.sard))
X.g3=array(X.g3.list,dim=c((n_annee-1),ncov.g.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AgO","NAO","Temp.jun/jul_oct","Zoo.jun/jul_oct"),c("1","2")))

X.g4.list=as.vector(c(X.g1.anch,X.g4.sard))
X.g4=array(X.g4.list,dim=c((n_annee-1),ncov.g.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AgO","NAO","Temp.jun/jul_oct","Zoo.jun/jul_oct"),c("1","2")))

X.g5.list=as.vector(c(X.g1.anch,X.g5.sard))
X.g5=array(X.g5.list,dim=c((n_annee-1),ncov.g.anch,2),dimnames=list(1:(n_annee-1),c("nsprat","AgO","NAO","Temp.jun/jul_oct","Zoo.jun/jul_oct"),c("1","2")))


X.ZSSN.anch<-cbind(X.ZSSN.anch,v9 = rep(NA,(n_annee-1)), V10 = rep(NA,(n_annee-1)),v11 = rep(NA,(n_annee-1)),v12 = rep(NA,(n_annee-1)))

X.ZSSN.list=as.vector(c(X.ZSSN.anch,X.ZSSN.sard))
X.ZSSN=array(X.ZSSN.list,dim=c((n_annee-1),ncov.zssn.sard,2),dimnames=list(1:(n_annee-1),c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12"),c("1","2")))


##------------observation data-----------##

data <- list(
  "X.M1"=X.M1,"X.M2"=X.M2,"X.M3"=X.M3,"X.M4"=X.M4,"X.M5"=X.M5,"X.g1"=X.g1,"X.g2"=X.g2,"X.g3"=X.g3,"X.g4"=X.g4,"X.g5"=X.g5,"X.ZSSN"=X.ZSSN,
  
  "ncov.M"=c(ncov.M.anch,ncov.M.sard),"ncov.g"=c(ncov.g.anch,ncov.g.sard),"ncov.zssn"=c(ncov.zssn.anch,ncov.zssn.sard),
  
  "n_annee" = length(data_anchovy$Year), # durée de l'analyse (en années)
  # Variables observées  
  "n1.obs" = matrix(c(data_anchovy$N1+1,data_sardine$N1+1),ncol=2), # Abondance age 1
  "n2.obs" = matrix(c(data_anchovy$N2+1,data_sardine$N2+1),ncol=2), # Abondance age 2
  "n3.obs" = matrix(c(data_anchovy$N3+1,data_sardine$N3+1),ncol=2), # Abondance age 3
  "n4.obs" = matrix(c(data_anchovy$N4plus+1,data_sardine$N4+1),ncol=2), # Abondance age 4,
  "n5.obs" = matrix(c(rep(NA,23),data_sardine$N5+1),ncol=2), # Abondance age 5,
  "n6.obs" = matrix(c(rep(NA,23),data_sardine$N6plus+1),ncol=2), # Abondance age 6,
  
  "nlarge.obs"= matrix(c(data_anchovy$N3+data_anchovy$N4plus,data_sardine$N2+data_sardine$N3+data_sardine$N4+data_sardine$N5+data_sardine$N6plus),ncol=2), # Abondance older fishes,
  "n.ssn.obs"= matrix(c(data_anchovy$N1+data_anchovy$N2+data_anchovy$N3+data_anchovy$N4plus,data_sardine$N1*maturite_sard[,1]+data_sardine$N2*maturite_sard[,2]+data_sardine$N3*maturite_sard[,3]+data_sardine$N4*maturite_sard[,4]+data_sardine$N5*maturite_sard[,5]+data_sardine$N6plus*maturite_sard[,6]),ncol=2), # Abondance mature fishes,
  
  "l1.obs" = matrix(c(data_anchovy$L1+1,data_sardine$L1+1),ncol=2), # Taille age 1
  "l2.obs" = matrix(c(data_anchovy$L2+1,data_sardine$L2+1),ncol=2), # Taille age 2
  "l3.obs" = matrix(c(data_anchovy$L3+1,data_sardine$L3+1),ncol=2), # Taille age 3
  "l4.obs" = matrix(c(data_anchovy$L4+1,data_sardine$L4+1),ncol=2), # Taille age 4
  "l5.obs" = matrix(c(rep(NA,23),data_sardine$L5+1),ncol=2), # Taille age 5
  "l6.obs" = matrix(c(rep(NA,23),data_sardine$L6+1),ncol=2), # Taille age 6
  
  "c1.obs" = matrix(c(data_anchovy$C1+1,data_sardine$C1+1),ncol=2), # Capture age 1
  "c2.obs" = matrix(c(data_anchovy$C2+1,data_sardine$C2+1),ncol=2), # Capture age 2
  "c3.obs" = matrix(c(data_anchovy$C3+1,data_sardine$C3+1),ncol=2), # Capture age 3
  "c4.obs" = matrix(c(rep(NA,23),data_sardine$C4+1),ncol=2), # Capture age 4
  "c5.obs" = matrix(c(rep(NA,23)+1,data_sardine$C5+1),ncol=2), # Capture age 5
  
  "mat" =mat,  # matrice de maturite aux ages 1 a 4 anchois + matrice de maturite aux ages 1 a 6 sardines
  
  # erreur observation dans les données 
  "sigma.n.obs" = matrix(c(data_anchovy$sigma_N,data_sardine$sigma_N),ncol=2), # erreur annuelle des abondances = ecartype tot (de ts les ages) pr chaq an
  "sigma.l1.obs" = matrix(c(data_anchovy$sigma_L1,data_sardine$sigma_L1),ncol=2), # erreur  des Longueur corporelle age 1
  "sigma.l2.obs" = matrix(c(data_anchovy$sigma_L2,data_sardine$sigma_L2),ncol=2), # erreur des Longueur corporelle age 2
  "sigma.l3.obs" = matrix(c(data_anchovy$sigma_L3,data_sardine$sigma_L3),ncol=2), # erreur des Longueur corporelle age 3
  "sigma.l4.obs" = matrix(c(data_anchovy$sigma_L4,data_sardine$sigma_L4),ncol=2), # erreur des Longueur corporelle age 4
  "sigma.l5.obs" = matrix(c(rep(NA,23),data_sardine$sigma_L5),ncol=2), # erreur des Longueur corporelle age 4
  "sigma.l6.obs" = matrix(c(rep(NA,23),data_sardine$sigma_L6),ncol=2) # erreur des Longueur corporelle age 4
  
  # Covariables
  
  
)
# "+1" sur chaque ligne pour éviter les données de 0 pour pouvoir tirer les valerus dans une lognormale


# --------------------------------------------------------
#                         Modèle
# --------------------------------------------------------

# lecture du modele a implementer 
setwd(paste(base,"04_MODELS/M3_conjoint",sep=""))
source("M3_joint_SSVS_rjags.txt") # nom du modele a aller chercher, ecrit ds un fichier txt 

model <- textConnection(model_simple) # on specifie a R qu'on utilise un language JAGS


# --------------------------------------------------------
#                  Variables stockées
# --------------------------------------------------------

variables_set <- c(
  # Variables d'états (variables estimées)
  "N1","N2","N3","N4","N5","N6","N.SSN", # Abondance
  "I1","I2","I3","I4","I5","I6","L.SSN", # Indice d'abondance absolue
  "L1","L2","L3","L4","L5","L6", # Longueur corporelle
  "C1","C2","C3","C4","C5", # Capture
  "M1","M2","M3","M4","M5", # mortalité naturelle
  "F1","F2","F3","F4","F5", # mortalité par pêche
  "Z1","Z2","Z3","Z4","Z5","Z.SSN", # mortalité totale
  "g1","g2","g3","g4","g5", # taux de croissance
  
  # Variabilité interannuelle
  "sigma.M1","sigma.M2","sigma.M3","sigma.M4","sigma.M5", # mortalité naturelle 
  "sigma.F1","sigma.F2","sigma.F3","sigma.F4","sigma.F5", # mortalité par pêche     
  "sigma.g1","sigma.g2","sigma.g3","sigma.g4","sigma.g5",# croissance
  "sigma.Z.SSN",# mortalité des ages matures
  
  "n1.pp","n2.pp","n3.pp","n4.pp","n5.pp","n6.pp",     # abondance predite /esperances de la prediction des abondances
  "c1.pp","c2.pp","c3.pp","c4.pp","c5.pp",     # captures predites /esperances de la prediction des captures
  "l1.pp","l2.pp","l3.pp","l4.pp","l5.pp","l6.pp",     # taille predite /esperances de la prediction des tailles
  
  
  # Coefficient directeur des régressions - /!\ SI MODELE CROISE 
  "alpha_M1",
  "alpha_M2","alpha.bis_M2",
  "alpha_M3","alpha.bis_M3","alpha.ter_M3", 
  "alpha_M4","alpha.bis_M4",
  "alpha_M5",
  
  "alpha_Z.SSN",
  
  "beta_M1","beta_M5",
  
  #var set pour SSVS
  "tau_in",
  
  "tau_beta.M1","tau_beta.M2","tau_beta.M3","tau_beta.M4","tau_beta.M5",
  "p_inclusion.M1","p_inclusion.M2","p_inclusion.M3","p_inclusion.M4","p_inclusion.M5",
  "beta.M1","beta.M2","beta.M3","beta.M4","beta.M5",
  
  "tau_beta.g1","tau_beta.g2","tau_beta.g3","tau_beta.g4","tau_beta.g5",
  "p_inclusion.g1","p_inclusion.g2","p_inclusion.g3","p_inclusion.g4","p_inclusion.g5",
  "beta.g1","beta.g2","beta.g3","beta.g4","beta.g5",
  
  "tau_beta.zssn",
  "p_inclusion.zssn",
  "beta.zssn"
)

# --------------------------------------------------------
#               Exécution du modèle
# --------------------------------------------------------

n_chains = 3 # nombre de chaines MCMC
n_adapt = 5000 # nombre d'itérations qui servent de burn-in = poubelle
# modèle croisé => n_adapt = 190000

time_compile <- system.time(
  model_jags <- jags.model(model, data=data, n.chain=n_chains, n.adapt=n_adapt)
)
time_compile

# --------------------------------------------------------
# Exécution des chaines de markov et stockage de l'échantillon
# --------------------------------------------------------

thin = 100 # on récupère une valeur toutes les 100 valeurs #100 par défaut
# modèle croisé => thin = 600
n_samples = 3000 # nombre d'itérations que l'on récupère par chaînes #3000 par défaut
n_iter = n_samples*thin # nombre d'itérations totales réalisées par le modèle
n.burnin = 5000


# pour runner le modèle en divisant mes 3 chaines MCMC sur 3 cores de mon ordi 
source("M3_joint_SSVS_rjags.txt") # nom du modele a aller chercher, ecrit ds un fichier txt 

PARA = T

# RUN ITERATIONS

if(PARA==T)
{
  
  cluster <- makePSOCKcluster(names=c(n_chains)) 
  time_compile <- system.time(
    res <- jags.parfit(cluster,data=data,params=variables_set,model="M3jointSSVS.txt",n.chains=n_chains,n.adapt=n_adapt,n.update=n.burnin,n.iter=n_iter,thin=thin,inits=NULL)
    #warnings()
  )
  time_compile
}


# --------------------------------------------------------
# sauvegarde des outputs
# --------------------------------------------------------
# creation d un dossier pour y stocker les outputs du modele

#tau_log_mu="tau_log_mu=0.1_"
nbiter=paste("iter=",n_iter,"_")
effet="_demo_SSVS.covar_sigma.g=0et5"

# _N1.anch.stand_N2.anch.stand_L1_sur_M1.anch__and__N1.sard.stand_L1_sur_M1.sard_and_N2.anch.stand_sur_M2.anch_and_N2.sard.stand_sur_M2.sard

setwd(paste(base,"05_OUTPUTS",sep=""))
st_wd <- paste(base,"05_OUTPUTS/M3_joint/",sep="")

DateFile = paste0(
  st_wd,
  nbiter,
  format(Sys.time(), "%d-%b-%Y %H.%M"),
  effet
)
dir.create(DateFile)
setwd(DateFile)
save(res, file = "M3_joint.Rdata")
write(model_simple,file("M3joint.SSVS.txt"))
dir.create("figures")
setwd(paste(DateFile, "/figures", sep = ""))
dir.create("posterior")
dir.create("time_series_variables")
dir.create("coefs")
setwd(paste(DateFile, "/figures/coefs", sep = ""))
dir.create("SSVS_anchovy")
dir.create("SSVS_sardine")
setwd(paste(DateFile, "/figures/posterior", sep = ""))
dir.create("anchovy")
dir.create("sardine")
setwd(paste(DateFile, "/figures/time_series_variables", sep = ""))
dir.create("anchovy")
dir.create("sardine")

