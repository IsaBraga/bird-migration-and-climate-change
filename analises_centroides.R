--------------------------------------------
  ##ANÁLISE DESLOCAMENTO DOS CENTRÓIDES##
--------------------------------------------
library(maptools)

nomes_index<-read.table("nomes_index.txt", h=T) #para a indexaçao 
resultados <- aves_migratorias_model_2$resultados_spp #RData originário do Modelo Bioclim


#listas das spp separadas por sistemas migratórios
spp_austrais<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_austrais.CSV", h=F)
spp_nearticos<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_nearticas.CSV", h=F)
spp_intra<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_intratropicais.CSV", h=F)


#reproducao
Cent0k.rep <- t(resultados[nomes_index$epoca=="0k" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",5:6,])
Cent26.rep <- t(resultados[nomes_index$epoca=="26" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",5:6,])
Cent85.rep <- t(resultados[nomes_index$epoca=="85" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",5:6,])


 
#invernada
Cent0k.inv <- t(resultados[nomes_index$epoca=="0k" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",5:6,])
Cent26.inv <- t(resultados[nomes_index$epoca=="26" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",5:6,])
Cent85.inv<- t(resultados[nomes_index$epoca=="85" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",5:6,])


#residencia
Cent0k.res <- t(resultados[nomes_index$epoca=="0k" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",5:6,])
Cent26.res <- t(resultados[nomes_index$epoca=="26" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",5:6,])
Cent85.res <- t(resultados[nomes_index$epoca=="85" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",5:6,])



## Diferenças entre as Centróides do presente e o futuro e entre os hemisferios
#(ângulo em graus da direção da mudança dos centroides)
library(rgeos)
library(maptools)
library(circular)

####REPRODUCAO
#calcular o azimute
Desloc_centRep_0k26<-matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.rep)[[1]],"azimute"))
for(j in 1: nrow(Cent0k.rep)){
  Desloc_centRep_0k26[j,1]<-gzAzimuth(from=Cent0k.rep[j,1:2],to=Cent26.rep[j,1:2])
}

Desloc_centRep_0k85<-matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.rep)[[1]],"azimute"))
for(j in 1: nrow(Cent0k.rep)){
  Desloc_centRep_0k85[j,1]<-gzAzimuth(from=Cent0k.rep[j,1:2],to=Cent85.rep[j,1:2])
}

##tornar os dados circulares e separar os sistemas migratórios
#cenário de baixa emissão
Desloc_centRep_0k26_Neart <- as.circular(Desloc_centRep_0k26[dimnames(Desloc_centRep_0k26)[[1]] %in% spp_nearticos$V1,], type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


Desloc_centRep_0k26_Aust<-as.circular(Desloc_centRep_0k26[dimnames(Desloc_centRep_0k26)[[1]] %in% spp_austrais$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))

Desloc_centRep_0k26_Intrat<-as.circular(Desloc_centRep_0k26[dimnames(Desloc_centRep_0k26)[[1]] %in% spp_intra$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


#cenário de alta emissão
Desloc_centRep_0k85_Neart <- as.circular(Desloc_centRep_0k85[dimnames(Desloc_centRep_0k85)[[1]] %in% spp_nearticos$V1,], type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


Desloc_centRep_0k85_Aust<-as.circular(Desloc_centRep_0k85[dimnames(Desloc_centRep_0k85)[[1]] %in% spp_austrais$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))

Desloc_centRep_0k85_Intrat<-as.circular(Desloc_centRep_0k85[dimnames(Desloc_centRep_0k85)[[1]] %in% spp_intra$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))




####INVERNADA
#calcular o azimute 
Desloc_centInv_0k26<-matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.inv)[[1]],"azimute"))
for(j in 1: nrow(Cent0k.inv)){
  Desloc_centInv_0k26[j,1]<-gzAzimuth(from=Cent0k.inv[j,1:2],to=Cent26.inv[j,1:2])
}

Desloc_centInv_0k85<-matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.inv)[[1]],"azimute"))
for(j in 1: nrow(Cent0k.inv)){
  Desloc_centInv_0k85[j,1]<-gzAzimuth(from=Cent0k.inv[j,1:2],to=Cent85.inv[j,1:2])
}



#retirar spp com centroides outliers: S. saya e E.minimus
Desloc_centInv_0k85[c(22,63),]<-NA


##tornar os dados circulares e separar os sistemas migratórios
#cenário de baixa emissão
Desloc_centInv_0k26_Neart <- as.circular(Desloc_centInv_0k26[dimnames(Desloc_centInv_0k26)[[1]] %in% spp_nearticos$V1,], type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


Desloc_centInv_0k26_Aust<-as.circular(Desloc_centInv_0k26[dimnames(Desloc_centInv_0k26)[[1]] %in% spp_austrais$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))

Desloc_centInv_0k26_Intrat<-as.circular(Desloc_centInv_0k26[dimnames(Desloc_centInv_0k26)[[1]] %in% spp_intra$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


#cenário de alta emissão
Desloc_centInv_0k85_Neart <- as.circular(Desloc_centInv_0k85[dimnames(Desloc_centInv_0k85)[[1]] %in% spp_nearticos$V1,], type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


Desloc_centInv_0k85_Aust<-as.circular(Desloc_centInv_0k85[dimnames(Desloc_centInv_0k85)[[1]] %in% spp_austrais$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))

Desloc_centInv_0k85_Intrat<-as.circular(Desloc_centInv_0k85[dimnames(Desloc_centInv_0k85)[[1]] %in% spp_intra$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))



####RESIDENCIA
#calcular o azimute
Desloc_centRes_0k26<-matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.res)[[1]],"azimute"))
for(j in 1: nrow(Cent0k.res)){
  Desloc_centRes_0k26[j,1]<-gzAzimuth(from=Cent0k.res[j,1:2],to=Cent26.res[j,1:2])
}

Desloc_centRes_0k85<-matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.res)[[1]],"azimute"))
for(j in 1: nrow(Cent0k.res)){
  Desloc_centRes_0k85[j,1]<-gzAzimuth(from=Cent0k.res[j,1:2],to=Cent85.res[j,1:2])
}

#retirar spp com centroides outliers: S. saya 
Desloc_centRes_0k85[63,]<-NA


##tornar os dados circulares e separar os sistemas migratórios
#cenário de baixa emissão
Desloc_centRes_0k26_Neart <- as.circular(Desloc_centRes_0k26[dimnames(Desloc_centRes_0k26)[[1]] %in% spp_nearticos$V1,], type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


Desloc_centRes_0k26_Aust<-as.circular(Desloc_centRes_0k26[dimnames(Desloc_centRes_0k26)[[1]] %in% spp_austrais$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))

Desloc_centRes_0k26_Intrat<-as.circular(Desloc_centRes_0k26[dimnames(Desloc_centRes_0k26)[[1]] %in% spp_intra$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


#cenário de alta emissão
Desloc_centRes_0k85_Neart <- as.circular(Desloc_centRes_0k85[dimnames(Desloc_centRes_0k85)[[1]] %in% spp_nearticos$V1,], type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))


Desloc_centRes_0k85_Aust<-as.circular(Desloc_centRes_0k85[dimnames(Desloc_centRes_0k85)[[1]] %in% spp_austrais$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))

Desloc_centRes_0k85_Intrat<-as.circular(Desloc_centRes_0k85[dimnames(Desloc_centRes_0k85)[[1]] %in% spp_intra$V1,],type=c("angles"),units=c("degrees"),template=c("none"),zero=1.5708,rotation=c("clock"))





#### análises com o teste de rayleigh juntando ângulo médio por sistema e pvalor
#spp neárticas
angulos_0k2685_Neart <- rbind(
  cbind(mean.circular(Desloc_centRep_0k26_Neart,na.rm=T)[[1]],  rayleigh.test(Desloc_centRep_0k26_Neart)$p.value),
  cbind(mean.circular(Desloc_centInv_0k26_Neart,na.rm=T)[[1]],  rayleigh.test(Desloc_centInv_0k26_Neart)$p.value),
  cbind(mean.circular(Desloc_centRes_0k26_Neart,na.rm=T)[[1]],  rayleigh.test(Desloc_centRes_0k26_Neart)$p.value),
  cbind(mean.circular(Desloc_centRep_0k85_Neart,na.rm=T)[[1]],  rayleigh.test(Desloc_centRep_0k85_Neart)$p.value),
  cbind(mean.circular(Desloc_centInv_0k85_Neart,na.rm=T)[[1]],  rayleigh.test(Desloc_centInv_0k85_Neart)$p.value),
  cbind(mean.circular(Desloc_centRes_0k85_Neart,na.rm=T)[[1]],  rayleigh.test(Desloc_centRes_0k85_Neart)$p.value))
rownames(angulos_0k2685_Neart)<-rep(c("Reprodutivas", "Não-reprodutivas", "Residência"), times=2)
colnames(angulos_0k2685_Neart)<-c( "Ângulo médio", "pvalor")

angulos_0k2685_Neart  <- ifelse(angulos_0k2685_Neart<0,360+angulos_0k2685_Neart, angulos_0k2685_Neart) #transformar os ângulos negativos em positivos


#spp austrais
angulos_0k2685_Aust <- rbind(
  cbind(mean.circular(Desloc_centRep_0k26_Aust,na.rm=T)[[1]],  rayleigh.test(Desloc_centRep_0k26_Aust)$p.value),
  cbind(mean.circular(Desloc_centInv_0k26_Aust,na.rm=T)[[1]],  rayleigh.test(Desloc_centInv_0k26_Aust)$p.value),
  cbind(mean.circular(Desloc_centRes_0k26_Aust,na.rm=T)[[1]],  rayleigh.test(Desloc_centRes_0k26_Aust)$p.value),
  cbind(mean.circular(Desloc_centRep_0k85_Aust,na.rm=T)[[1]],  rayleigh.test(Desloc_centRep_0k85_Aust)$p.value),
  cbind(mean.circular(Desloc_centInv_0k85_Aust,na.rm=T)[[1]],  rayleigh.test(Desloc_centInv_0k85_Aust)$p.value),
  cbind(mean.circular(Desloc_centRes_0k85_Aust,na.rm=T)[[1]],  rayleigh.test(Desloc_centRes_0k85_Aust)$p.value))
rownames(angulos_0k2685_Aust)<-rep(c("Reprodutivas", "Não-reprodutivas", "Residência"), times=2)
colnames(angulos_0k2685_Aust)<-c( "Ângulo médio", "pvalor")

angulos_0k2685_Aust<-ifelse(angulos_0k2685_Aust<0,360+angulos_0k2685_Aust,angulos_0k2685_Aust) #transformar os ângulos negativos em positivos


#spp Intratropicais
angulos_0k2685_Intrat <- rbind(
  cbind(mean.circular(Desloc_centRep_0k26_Intrat,na.rm=T)[[1]],  rayleigh.test(Desloc_centRep_0k26_Intrat)$p.value),
  cbind(mean.circular(Desloc_centInv_0k26_Intrat,na.rm=T)[[1]],  rayleigh.test(Desloc_centInv_0k26_Intrat)$p.value),
  cbind(mean.circular(Desloc_centRes_0k26_Intrat,na.rm=T)[[1]],  rayleigh.test(Desloc_centRes_0k26_Intrat)$p.value),
  cbind(mean.circular(Desloc_centRep_0k85_Intrat,na.rm=T)[[1]],  rayleigh.test(Desloc_centRep_0k85_Intrat)$p.value),
  cbind(mean.circular(Desloc_centInv_0k85_Intrat,na.rm=T)[[1]],  rayleigh.test(Desloc_centInv_0k85_Intrat)$p.value),
  cbind(mean.circular(Desloc_centRes_0k85_Intrat,na.rm=T)[[1]],  rayleigh.test(Desloc_centRes_0k85_Intrat)$p.value))
rownames(angulos_0k2685_Intrat)<-rep(c("Reprodutivas", "Não-reprodutivas", "Residência"), times=2)
colnames(angulos_0k2685_Intrat)<-c( "Ângulo médio", "pvalor")

angulos_0k2685_Intrat<-ifelse(angulos_0k2685_Intrat<0,360+angulos_0k2685_Intrat,angulos_0k2685_Intrat) #transformar os ângulos negativos em positivos

angulos_0k2685<-cbind(angulos_0k2685_Neart, angulos_0k2685_Intrat, angulos_0k2685_Aust)

write.csv2(angulos_0k2685 , "angulos_0k2685.csv")





#####DISTANCIAS ENTRE OS CENTROIDES PRESENTE E FUTURU PARA AS MESMAS ÁREAS

#calculo das distancias (km) utilizando a funcao distVincentyEllipsoid
library(geosphere)
library(bbmle)

###REPRODUÇÃO
DistRep0k_26 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.rep)[[1]],"distância (km)"))
for(i in 1:nrow(DistRep0k_26)){
  DistRep0k_26[i,] <- distVincentyEllipsoid(p1=Cent0k.rep[i,], p2=Cent26.rep[i,])/1000
}

DistRep0k_85 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.rep)[[1]],"distância (km)"))
for(i in 1:nrow(DistRep0k_85)){
  DistRep0k_85[i,] <- distVincentyEllipsoid(p1=Cent0k.rep[i,], p2=Cent85.rep[i,])/1000
}

###INVERNADA
DistInv0k_26 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.inv)[[1]],"distância (km)"))
for(i in 1:nrow(DistInv0k_26)){
  DistInv0k_26[i,] <- distVincentyEllipsoid(p1=Cent0k.inv[i,], p2=Cent26.inv[i,])/1000
}

DistInv0k_85 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.inv)[[1]],"distância (km)"))
for(i in 1:nrow(DistInv0k_85)){
  DistInv0k_85[i,] <- distVincentyEllipsoid(p1=Cent0k.inv[i,], p2=Cent85.inv[i,])/1000
}

##retirar outliers do rcp85
library(outliers)
outlier(DistInv0k_85, logical=T) #Sayornis saya
DistInv0k_85[63,]<-NA

##?? quando aplico a funcao outlier no objeto DistInv0k_85, aparce um outlier, a spp S. saya, retiro ela e aplico a funcao novamente, então aparece mais um outlier: a spp E. minimus.
DistInv0k_85[22,]<-NA


###RESIDÊNCIA
DistRes0k_26 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.res)[[1]],"distância (km)"))
for(i in 1:nrow(DistRes0k_26)){
  DistRes0k_26[i,] <- distVincentyEllipsoid(p1=Cent0k.res[i,], p2=Cent26.res[i,])/1000
}

DistRes0k_85 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(Cent0k.res)[[1]],"distância (km)"))
for(i in 1:nrow(DistRes0k_85)){
  DistRes0k_85[i,] <- distVincentyEllipsoid(p1=Cent0k.res[i,], p2=Cent85.res[i,])/1000
}

##retirar outliers do rcp85
library(outliers)
outlier(DistRes0k_85, logical=T) #Sayornis saya
DistRes0k_85[63,]<-NA

###separar spp entre os sistemas

##cenario de baixa emissão
DistRep0k_26Aust <-DistRep0k_26[dimnames(DistRep0k_26)[[1]] %in% spp_austrais$V1,]
DistRep0k_26Neart <-DistRep0k_26[dimnames(DistRep0k_26)[[1]] %in% spp_nearticos$V1,]
DistRep0k_26Intrat <-DistRep0k_26[dimnames(DistRep0k_26)[[1]] %in% spp_intra$V1,]

DistInv0k_26Aust <-DistInv0k_26[dimnames(DistInv0k_26)[[1]] %in% spp_austrais$V1,]
DistInv0k_26Neart <-DistInv0k_26[dimnames(DistInv0k_26)[[1]] %in% spp_nearticos$V1,]
DistInv0k_26Intrat <-DistInv0k_26[dimnames(DistInv0k_26)[[1]] %in% spp_intra$V1,]

DistRes0k_26Aust <-DistRes0k_26[dimnames(DistRes0k_26)[[1]] %in% spp_austrais$V1,]
DistRes0k_26Neart <-DistRes0k_26[dimnames(DistRes0k_26)[[1]] %in% spp_nearticos$V1,]
DistRes0k_26Intrat <-DistRes0k_26[dimnames(DistRes0k_26)[[1]] %in% spp_intra$V1,]


##CENÁRIO DE ALTA EMISSÃO

DistRep0k_85Aust <-DistRep0k_85[dimnames(DistRep0k_85)[[1]] %in% spp_austrais$V1,]
DistRep0k_85Neart <-DistRep0k_85[dimnames(DistRep0k_85)[[1]] %in% spp_nearticos$V1,]
DistRep0k_85Intrat <-DistRep0k_85[dimnames(DistRep0k_85)[[1]] %in% spp_intra$V1,]

DistInv0k_85Aust <-DistInv0k_85[dimnames(DistInv0k_85)[[1]] %in% spp_austrais$V1,]
DistInv0k_85Neart <-DistInv0k_85[dimnames(DistInv0k_85)[[1]] %in% spp_nearticos$V1,]
DistInv0k_85Intrat <-DistInv0k_85[dimnames(DistInv0k_85)[[1]] %in% spp_intra$V1,]

DistRes0k_85Aust <-DistRes0k_85[dimnames(DistRes0k_85)[[1]] %in% spp_austrais$V1,]
DistRes0k_85Neart <-DistRes0k_85[dimnames(DistRes0k_85)[[1]] %in% spp_nearticos$V1,]
DistRes0k_85Intrat <-DistRes0k_85[dimnames(DistRes0k_85)[[1]] %in% spp_intra$V1,]




##ANÁLISE DO DESLOCAMENTOS DOS CENTRÓIDES (KM) modelo linear

##CENÁRIO DE BAIXA EMISSÃO
#ÁREAS DE REPRODUÇÃO
reprodutivos26<-data.frame(especies = names(c(DistRep0k_26Neart, DistRep0k_26Intrat,DistRep0k_26Aust)),
              distancia = c(DistRep0k_26Neart, DistRep0k_26Intrat,DistRep0k_26Aust), 
              sistema = c(rep("Neártico-Neotropical",each=29), rep("Intratropical",each=11), rep("Austral-Neotropical",40)), row.names = NULL)


invernada26<-data.frame(especies = names(c(DistInv0k_26Neart, DistInv0k_26Intrat,DistInv0k_26Aust)),
            distancia = c(DistInv0k_26Neart, DistInv0k_26Intrat,DistInv0k_26Aust), 
            sistema = c(rep("Neártico-Neotropical",each=29), rep("Intratropical",each=11), rep("Austral-Neotropical",40)), row.names = NULL)

residencia26<-data.frame(especies = names(c(DistRes0k_26Neart, DistRes0k_26Intrat,DistRes0k_26Aust)),
            distancia = c(DistRes0k_26Neart, DistRes0k_26Intrat,DistRes0k_26Aust), 
            sistema = c(rep("Neártico-Neotropical",each=29), rep("Intratropical",each=11), rep("Austral-Neotropical",40)), row.names = NULL)



##CENÁRIO DE ALTA EMISSÃO
reprodutivos85<-data.frame(especies = names(c(DistRes0k_26Neart, DistRes0k_26Intrat,DistRes0k_26Aust)),
              distancia = c(DistRep0k_85Neart, DistRep0k_85Intrat,DistRep0k_85Aust), 
              sistema = c(rep("Neártico-Neotropical",each=29), rep("Intratropical",each=11), rep("Austral-Neotropical",40)), row.names = NULL)

invernada85<-data.frame(especies = names(c(DistInv0k_85Neart, DistInv0k_85Intrat,DistInv0k_85Aust)),
              distancia = c(DistInv0k_85Neart, DistInv0k_85Intrat,DistInv0k_85Aust), 
              sistema = c(rep("Neártico-Neotropical",each=29), rep("Intratropical",each=11), rep("Austral-Neotropical",40)), row.names = NULL)


residencia85<-data.frame(especies = names(c(DistRes0k_85Neart, DistRes0k_85Intrat,DistRes0k_85Aust)),
              distancia = c(DistRes0k_85Neart, DistRes0k_85Intrat,DistRes0k_85Aust), 
              sistema = c(rep("Neártico-Neotropical",each=29), rep("Intratropical",each=11), rep("Austral-Neotropical",40)), row.names = NULL)

         
         
##análise linear model do deslocamentos dos centroides (km)
##REPRODUÇÃO
##RCP 2.6
distCent_lm0_Rep26<-lm(distancia ~ 0, data = reprodutivos26)
distCent_lm1_Rep26<-lm(distancia ~ 1, data = reprodutivos26)
distCent_lmSist_Rep26<-lm(distancia ~ sistema, data = reprodutivos26)
library(bbmle)
AICtab(distCent_lm0_Rep26, distCent_lm1_Rep26, distCent_lmSist_Rep26)

dAIC df
distCent_lm1_Rep26     0.0 2 
distCent_lmSist_Rep26  2.8 4 
distCent_lm0_Rep26    63.8 1 

coefficients(distCent_lm1_Rep26)
(Intercept) 
227.1459 #média do deslocamento (km) dos centróides para todos os sistemas

##RCP 8.5
library(lme4)
library(bbmle)
distCent_lm0_Rep85<-lm(distancia ~ 0, data = reprodutivos85)
distCent_lm1_Rep85<-lm(distancia ~ 1, data = reprodutivos85)
distCent_lmSist_Rep85<-lm(distancia ~ sistema, data = reprodutivos85)
AICtab(distCent_lm0_Rep85, distCent_lm1_Rep85, distCent_lmSist_Rep85)

dAIC df
distCent_lmSist_Rep85  0.0 4 
distCent_lm1_Rep85     9.4 2 
distCent_lm0_Rep85    61.8 1 

plot(distCent_lmSist_Rep85)
coefficients(distCent_lmSist_Rep85)
(Intercept)        sistemaIntratropical sistemaNeártico-Neotropical 
404.0179                   -133.7514                    332.5178 

#média do deslocamento (km) dos centróides para sistema austral = 404.0179 
#média do deslocamento (km) dos centróides para sistema intratropical 404.0179 - 133.7514 =  270.2665
#média do deslocamento (km) dos centróides para sistema neártico  404.0179 + 332.5178 = 736.5357


TukeyHSD(aov(distCent_lmSist_Rep85))
$sistema
                                            diff        lwr      upr     p adj
Intratropical-Austral-Neotropical        -133.7514 -655.85878 388.3560 0.8267845
Neártico-Neotropical-Austral-Neotropical  332.5178  102.11704 562.9186 0.0032636
Neártico-Neotropical-Intratropical        466.2692  -51.69991 984.2383 0.0847261

IC neartico = > 1.96*94.26
[1] 184.7496

IC austral = > 69.30*1.96
[1] 135.828

summary(distCent_lmSist_Rep85)

Call:
  lm(formula = distancia ~ sistema, data = reprodutivos85)

Residuals:
  Min      1Q  Median      3Q     Max 
-597.52 -154.38  -56.96  166.04  929.95 

Coefficients:
                              Estimate Std. Error t value Pr(>|t|)    
(Intercept)                   404.02      69.30   5.830 1.17e-06 ***
  sistemaIntratropical       -133.75     213.60  -0.626  0.53515    
sistemaNeártico-Neotropical   332.52      94.26   3.528  0.00117 ** 
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 285.7 on 36 degrees of freedom
(41 observations deleted due to missingness)
Multiple R-squared:  0.2903,	Adjusted R-squared:  0.2509 
F-statistic: 7.364 on 2 and 36 DF,  p-value: 0.002084



##INVERNADA 
#RCP2.6
distCent_lm0_Inv26<-lm(distancia ~ 0, data = invernada26)
distCent_lm1_Inv26<-lm(distancia ~ 1, data = invernada26)
distCent_lmSist_Inv26<-lm(distancia ~ sistema, data = invernada26)

AICtab(distCent_lm0_Inv26, distCent_lm1_Inv26, distCent_lmSist_Inv26)
dAIC df
distCent_lm1_Inv26     0.0 2 
distCent_lmSist_Inv26  1.3 4 
distCent_lm0_Inv26    19.5 1 

coefficients(distCent_lm1_Inv26)
(Intercept) 
391.9026 #média do deslocamento (km) dos centróides para todos os sistemas

#TukeyHSD(aov(distCent_lmSist_Inv26))
library(bbmle)
#RCP8.5
distCent_lm0_Inv85<-lm(distancia ~ 0, data = invernada85)
distCent_lm1_Inv85<-lm(distancia ~ 1, data = invernada85)
distCent_lmSist_Inv85<-lm(distancia ~ sistema, data = invernada85)

AICtab(distCent_lm0_Inv85, distCent_lm1_Inv85, distCent_lmSist_Inv85)
dAIC df
distCent_lm1_Inv85     0.0 2 
distCent_lmSist_Inv85  3.6 4 
distCent_lm0_Inv85    61.6 1 

plot(distCent_lm1_Inv85)
coefficients(distCent_lm1_Inv85)
(Intercept) 
835.4236  #média do deslocamento (km) dos centróides para todos os sistemas


##RESIDENCIA
#RCP2.6
distCent_lm0_Res26<-lm(distancia ~ 0, data = residencia26)
distCent_lm1_Res26<-lm(distancia ~ 1, data = residencia26)
distCent_lmSist_Res26<-lm(distancia ~ sistema, data = residencia26)

AICtab(distCent_lm0_Res26, distCent_lm1_Res26, distCent_lmSist_Res26)
dAIC df
distCent_lm1_Res26     0.0 2 
distCent_lmSist_Res26  1.4 4 
distCent_lm0_Res26    45.9 1

coefficients(distCent_lm1_Res26)
(Intercept) 
298.2125
         
  

#RCP8.5
distCent_lm0_Res85<-lm(distancia ~ 0, data = residencia85)
distCent_lm1_Res85<-lm(distancia ~ 1, data = residencia85)
distCent_lmSist_Res85<-lm(distancia ~ sistema, data = residencia85)

AICtab(distCent_lm0_Res85, distCent_lm1_Res85, distCent_lmSist_Res85)
dAIC df
distCent_lm1_Res85     0.0 2 
distCent_lmSist_Res85  3.8 4 
distCent_lm0_Res85    47.7 1 

plot(distCent_lm1_Res85)
hat values (leverages) are all = 0.02631579
and there are no factor predictors; no plot no. 5

coefficients(distCent_lm1_Res85)
(Intercept) 
777.4691  

      
-------------------
##Maraísa Braga
##Leo Ré Jorge  
-------------------  

