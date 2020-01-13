-------------------------------------
  ##ANÁLISE DISTANCIA DE MIGRAÇÃO##
-------------------------------------

#distancias de migracao!! **a cereja do bolo**

library(geosphere)
library("bbmle")


nomes_index<-read.table("nomes_index.txt", h=T) #para a indexaçao 
resultados <- aves_migratorias_model_2$resultados_spp #RData originário do Modelo Bioclim


#listas das spp separadas por sistemas migratórios
spp_austrais<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_austrais.CSV", h=F)
spp_nearticos<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_nearticas.CSV", h=F)
spp_intra<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_intratropicais.csv", h=F)



#reproducao
midRep_0k <- t(resultados[nomes_index$tempo=="0k" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",5:6,])
midRep_26 <- t(resultados[nomes_index$tempo=="26" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",5:6,])
midRep_85 <- t(resultados[nomes_index$tempo=="85" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",5:6,])

#invernada
midInv_0k <- t(resultados[nomes_index$tempo=="0k" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",5:6,])
midInv_26 <- t(resultados[nomes_index$tempo=="26" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",5:6,])
midInv_85 <- t(resultados[nomes_index$tempo=="85" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",5:6,])

#residencia
midRes_0k <- t(resultados[nomes_index$tempo=="0k" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",5:6,])
midRes_26 <- t(resultados[nomes_index$tempo=="26" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",5:6,])
midRes_85 <- t(resultados[nomes_index$tempo=="85" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",5:6,])




library(geosphere)

#distancias ente reproducao e invernada (km)
DistInvRep0k <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midInv_0k)[[1]], "distancia mig"))
for(j in 1:nrow(DistInvRep0k)){
  DistInvRep0k[j,1] <- distVincentyEllipsoid(p1=midInv_0k[j,1:2], p2=midRep_0k[j,1:2])/1000
} 


DistInvRep26 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midInv_26)[[1]], "distancia mig"))
for(j in 1:nrow(DistInvRep26)){
  DistInvRep26[j,1] <- distVincentyEllipsoid(p1=midInv_26[j,1:2], p2=midRep_26[j,1:2])/1000
} 


DistInvRep85 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midInv_85)[[1]], "distancia mig"))
for(j in 1:nrow(DistInvRep85)){
  DistInvRep85[j,1] <- distVincentyEllipsoid(p1=midInv_85[j,1:2], p2=midRep_85[j,1:2])/1000
}

##**retirar a spp Sayornis saya pois ela tá zuando o barraco. tem outliers nas áreas de residência e de invernada no deslocamento de seus centroides.Retirar a spp Empidonax minimus, tem outlier no deslocamento do centroide de invernada.

DistInvRep85[63,]<-NA
DistInvRes0k[63,]<-NA
DistInvRep85[22,]<-NA
DistInvRes0k[22,]<-NA

#distancia de migracao entre invernada e residencia
DistInvRes0k <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midInv_0k)[[1]], "distancia mig"))
for(j in 1:nrow(DistInvRes0k)){
  DistInvRes0k[j,1] <- distVincentyEllipsoid(p1=midInv_0k[j,1:2], p2=midRes_0k[j,1:2])/1000
} 


DistInvRes26 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midInv_26)[[1]], "distancia mig"))
for(j in 1:nrow(DistInvRes26)){
  DistInvRes26[j,1] <- distVincentyEllipsoid(p1=midInv_26[j,1:2], p2=midRes_26[j,1:2])/1000
} 


DistInvRes85 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midInv_85)[[1]], "distancia mig"))
for(j in 1:nrow(DistInvRes85)){
  DistInvRes85[j,1] <- distVincentyEllipsoid(p1=midInv_85[j,1:2], p2=midRes_85[j,1:2])/1000
}

##**retirar a spp Sayornis saya pois ela tá zuando o barraco. tem outliers nas áreas de residência e de invernada no deslocamento de seus centroides.

DistInvRes85[63,]<-NA
DistInvRes0k[63,]<-NA

#distancia de migracao entre reproducao e residencia
DistRepRes0k <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midRep_0k)[[1]], "distancia mig"))
for(j in 1:nrow(DistRepRes0k)){
  DistRepRes0k[j,1] <- distVincentyEllipsoid(p1=midRep_0k[j,1:2], p2=midRes_0k[j,1:2])/1000
} 


DistRepRes26 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midRep_26)[[1]], "distancia mig"))
for(j in 1:nrow(DistRepRes26)){
  DistRepRes26[j,1] <- distVincentyEllipsoid(p1=midRep_26[j,1:2], p2=midRes_26[j,1:2])/1000
} 


DistRepRes85 <- matrix(nrow=82, ncol=1,dimnames=list(dimnames(midRep_85)[[1]], "distancia mig"))
for(j in 1:nrow(DistRepRes85)){
  DistRepRes85[j,1] <- distVincentyEllipsoid(p1=midRep_85[j,1:2], p2=midRes_85[j,1:2])/1000
}

##**retirar a spp Sayornis saya pois ela tá zuando o barraco. tem outliers nas áreas de residência e de invernada no deslocamento de seus centroides.

DistRepRes85[63,]<-NA
DistRepRes0k[63,]<-NA

##testar se as distancias de migracao mudaram do presente para o futuro
library(lme4)
library(bbmle)



#REPRODUCAO - INVERNADA
#RCP26 - separar entre os sistemas migratorios
DistInvRep_0k26_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                            tempo= rep(c("0k", "26"), each=29), 
                            especie = spp_nearticos$V1,
                            distancia=c(DistInvRep0k[dimnames(DistInvRep0k)[[1]] %in% spp_nearticos$V1,],
                                        DistInvRep26[dimnames(DistInvRep26)[[1]] %in% spp_nearticos$V1,]))

DistInvRep_0k26_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                                  tempo= rep(c("0k", "26"), each=40), 
                                  especie = spp_austrais$V1,
                                  distancia=c(DistInvRep0k[dimnames(DistInvRep0k)[[1]] %in% spp_austrais$V1,],
                                              DistInvRep26[dimnames(DistInvRep26)[[1]] %in% spp_austrais$V1,]))

DistInvRep_0k26_Intrat<-data.frame(sistema=rep("intrat", each=22),
                                 tempo= rep(c("0k", "26"), each=11), 
                                 especie = spp_intra$V1,
                                 distancia=c(DistInvRep0k[dimnames(DistInvRep0k)[[1]] %in% spp_intra$V1,],
                                             DistInvRep26[dimnames(DistInvRep26)[[1]] %in% spp_intra$V1,]))


DistInvRep_0k26_lmer<-rbind(DistInvRep_0k26_Aust,DistInvRep_0k26_Neart, DistInvRep_0k26_Intrat)


#modelo misto
DistInvRep26_tempo<-lmer(distancia ~ tempo + (1|especie) ,data=DistInvRep_0k26_lmer,REML=F)
DistInvRep26_nulo<-lmer(distancia ~ 1 + (1|especie) ,data=DistInvRep_0k26_lmer,REML=F)
DistInvRep26_tempoSistema<-lmer(distancia ~ tempo + sistema + (1|especie) ,data=DistInvRep_0k26_lmer,REML=F)
DistInvRep26_sist<-lmer(distancia ~ sistema + (1|especie) ,data=DistInvRep_0k26_lmer,REML=F)
DistInvRep26_interacao<-lmer(distancia ~ tempo * sistema + (1|especie) ,data=DistInvRep_0k26_lmer,REML=F)

AICtab(DistInvRep26_tempo,DistInvRep26_nulo,DistInvRep26_tempoSistema,DistInvRep26_sist,DistInvRep26_interacao) 

dAIC df
DistInvRep26_tempoSistema  0.0 6 
DistInvRep26_interacao     2.8 8 
DistInvRep26_sist          3.6 5 
DistInvRep26_tempo         6.9 4 
DistInvRep26_nulo         10.2 3 

plot(DistInvRep26_tempoSistema)
fixef(DistInvRep26_tempoSistema)
(Intercept)           tempo26 sistemaneart_neot     sistemaintrat 
1943.9577          294.3598         1749.0526          589.8184 


##REFAZER....
#distância de migraçao austral-neotropical no presente = 1944
#distância de migraçao austral-neotropical no futuro 1943.9577 + 294.3598  = 2238.318
(2238.318*100)/1943.9577 - 100 #aumento de 15% nas distâncias mig

#distância de migraçao neártico-neotropical no presente 1943.9577 + 1749.0526   = 3693.01
#distância de migraçao neártico-neotropical no futuro 1943.9577 + 294.3598 + 1749.0526 =  3987.37
(3987.37*100)/3693.01 - 100 #aumento de 7% nas distâncias mig

#distância de migraçao intratropical no presente 1943.9577 + 589.8184    = 2533.776
#distância de migraçao intratropical no futuro 1943.9577 + 294.3598 + 589.8184  =  2828.136
(2828.136*100)/2533.776 - 100 #aumento de 11% nas distâncias mig


library(emmeans)
emmeans(DistInvRep26_tempoSistema, pairwise ~ sistema, adjust = "Tukey")
$emmeans
sistema    emmean  SE  df asymp.LCL asymp.UCL
aust_neot    2091 354 Inf      1396      2786
neart_neot   3978 339 Inf      3314      4643
intrat       2681 934 Inf       850      4512

Results are averaged over the levels of: tempo 
Degrees-of-freedom method: asymptotic 
Confidence level used: 0.95 

$contrasts
contrast               estimate  SE  df z.ratio p.value
aust_neot - neart_neot    -1887 490 Inf -3.849  0.0003 
aust_neot - intrat         -590 999 Inf -0.591  0.8251 
neart_neot - intrat        1297 994 Inf  1.306  0.3920 

Results are averaged over the levels of: tempo 
P value adjustment: tukey method for comparing a family of 3 estimates 


#rcp85 - separar entre os sistemas migratorios
library(lme4)
library(bbmle)
DistInvRep_0k85_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                                  tempo= rep(c("0k", "85"), each=29), 
                                  especie = spp_nearticos$V1,
                                  distancia=c(DistInvRep0k[dimnames(DistInvRep0k)[[1]] %in% spp_nearticos$V1,],
                                              DistInvRep85[dimnames(DistInvRep85)[[1]] %in% spp_nearticos$V1,]))

DistInvRep_0k85_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                                 tempo= rep(c("0k", "85"), each=40), 
                                 especie = spp_austrais$V1,
                                 distancia=c(DistInvRep0k[dimnames(DistInvRep0k)[[1]] %in% spp_austrais$V1,],
                                             DistInvRep85[dimnames(DistInvRep85)[[1]] %in% spp_austrais$V1,]))

DistInvRep_0k85_Intrat<-data.frame(sistema=rep("intrat", each=22),
                                   tempo= rep(c("0k", "85"), each=11), 
                                   especie = spp_intra$V1,
                                   distancia=c(DistInvRep0k[dimnames(DistInvRep0k)[[1]] %in% spp_intra$V1,],
                                               DistInvRep85[dimnames(DistInvRep85)[[1]] %in% spp_intra$V1,]))





DistInvRep_0k85_lmer<-rbind(DistInvRep_0k85_Aust,DistInvRep_0k85_Neart, DistInvRep_0k85_Intrat)
teste<-DistInvRep_0k85_lmer[complete.cases(DistInvRep_0k85_lmer),]
table(teste$tempo, teste$sistema)

    aust_neot neart_neot intrat
0k        21         24      3
85        13         13      2


DistInvRep85_tempo<-lmer(distancia ~ tempo + (1|especie) ,data=DistInvRep_0k85_lmer,REML=F)
DistInvRep85_nulo<-lmer(distancia ~ 1 + (1|especie) ,data=DistInvRep_0k85_lmer,REML=F)
DistInvRep85_tempoSistema<-lmer(distancia ~ tempo + sistema + (1|especie) ,data=DistInvRep_0k85_lmer,REML=F)
DistInvRep85_sist<-lmer(distancia ~ sistema + (1|especie) ,data=DistInvRep_0k85_lmer,REML=F)
DistInvRep85_interacao<-lmer(distancia ~ tempo * sistema + (1|especie) ,data=DistInvRep_0k85_lmer,REML=F)

AICtab(DistInvRep85_tempo,DistInvRep85_nulo,DistInvRep85_tempoSistema,DistInvRep85_sist,DistInvRep85_interacao) 
dAIC df
DistInvRep85_interacao     0.0 8 
DistInvRep85_tempoSistema  6.2 6 
DistInvRep85_sist         13.6 5 
DistInvRep85_tempo        14.7 4 
DistInvRep85_nulo         21.3 3 

plot(DistInvRep85_interacao)
fixef(DistInvRep85_interacao)
(Intercept)                   tempo85 
2006.3308                   15.0660 
sistemaneart_neot             sistemaintrat 
1638.2381                  478.3634 
tempo85:sistemaneart_neot     tempo85:sistemaintrat 
1025.5589                  433.2820 

#distância de migraçao austral-neotropical no presente = 2006.3308 
#distância de migraçao austral-neotropical no futuro 2006.33083  + 15.0660  =  2021.397
(2021.397*100)/2006.33083 - 100 #aumento de 0.75% nas distâncias mig
2021.397 - 2006.3308 = 15 km

#distância de migraçao neártico-neotropical no presente 2006.33083  + 1638.2381   =   3644.569
#distância de migraçao neártico-neotropical no futuro 2006.33083  + 15.0660 + 1638.2381 + 1025.5589    =  4685.194
(4685.194*100)/ 3644.569 - 100 #aumento de 28,6% nas distâncias mig = 1040 km
4685.194 - 3644.569  = 1040.625 km

#distância de migraçao intratropical no presente 2006.3308 + 478.3634     = 2484.694
#distância de migraçao intratropical no futuro 2006.33083 +15.0660 + 478.3634 + 433.2820   =  2933.042
(2933.042*100)/2484.694 - 100 #aumento de 18% nas distâncias mig
2933-2485 =  448.348 km

emmeans(DistInvRep85_interacao, pairwise ~ sistema|tempo, adjust = "Tukey")
$emmeans
tempo = 0k:
  sistema    emmean   SE   df lower.CL upper.CL
aust_neot    2006  388 54.5     1229     2784
neart_neot   3645  363 54.5     2917     4372
intrat       2485 1026 54.4      427     4542

tempo = 85:
  sistema    emmean   SE   df lower.CL upper.CL
aust_neot    2021  412 65.2     1199     2844
neart_neot   4685  393 68.7     3901     5470
intrat       2933 1078 63.3      779     5087

Degrees-of-freedom method: kenward-roger 
Confidence level used: 0.95 

$contrasts
tempo = 0k:
  contrast               estimate   SE   df t.ratio p.value
aust_neot - neart_neot    -1638  531 54.5 -3.084  0.0089 
aust_neot - intrat         -478 1097 54.4 -0.436  0.9007 
neart_neot - intrat        1160 1089 54.4  1.065  0.5394 

tempo = 85:
  contrast               estimate   SE   df t.ratio p.value
aust_neot - neart_neot    -2664  569 66.9 -4.678  <.0001 
aust_neot - intrat         -912 1154 63.6 -0.790  0.7105 
neart_neot - intrat        1752 1148 64.0  1.527  0.2853 

Degrees-of-freedom method: kenward-roger 
P value adjustment: tukey method for comparing a family of 3 estimates 


#REPRODUCAO - RESIDENCIA
library(lme4)
library(bbmle)
#RCP26
DistRepRes_0k26_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                                  tempo= rep(c("0k", "26"), each=29), 
                                  especie = spp_nearticos$V1,
                                  distancia=c(DistRepRes0k[dimnames(DistRepRes0k)[[1]] %in% spp_nearticos$V1,],
                                              DistRepRes26[dimnames(DistRepRes26)[[1]] %in% spp_nearticos$V1,]))

DistRepRes_0k26_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                                 tempo= rep(c("0k", "26"), each=40), 
                                 especie = spp_austrais$V1,
                                 distancia=c(DistRepRes0k[dimnames(DistRepRes0k)[[1]] %in% spp_austrais$V1,],
                                             DistRepRes26[dimnames(DistRepRes26)[[1]] %in% spp_austrais$V1,]))

DistRepRes_0k26_Intrat<-data.frame(sistema=rep("intrat", each=22),
                                   tempo= rep(c("0k", "26"), each=11), 
                                   especie = spp_intra$V1,
                                   distancia=c(DistRepRes0k[dimnames(DistRepRes0k)[[1]] %in% spp_intra$V1,],
                                               DistRepRes26[dimnames(DistRepRes26)[[1]] %in% spp_intra$V1,]))


DistRepRes_0k26_lmer<-rbind(DistRepRes_0k26_Aust,DistRepRes_0k26_Neart, DistRepRes_0k26_Intrat)


#modelo misto
DistRepRes26_tempo<-lmer(distancia ~ tempo + (1|especie) ,data=DistRepRes_0k26_lmer,REML=F)
DistRepRes26_nulo<-lmer(distancia ~ 1 + (1|especie) ,data=DistRepRes_0k26_lmer,REML=F)
DistRepRes26_tempoSistema<-lmer(distancia ~ tempo + sistema + (1|especie) ,data=DistRepRes_0k26_lmer,REML=F)
DistRepRes26_sist<-lmer(distancia ~ sistema + (1|especie) ,data=DistRepRes_0k26_lmer,REML=F)
DistRepRes26_interacao<-lmer(distancia ~ tempo * sistema + (1|especie) ,data=DistRepRes_0k26_lmer,REML=F)

AICtab(DistRepRes26_tempo,DistRepRes26_nulo,DistRepRes26_tempoSistema,DistRepRes26_sist,DistRepRes26_interacao) 
dAIC df
DistRepRes26_nulo         0.0  3 
DistRepRes26_tempo        2.0  4 
DistRepRes26_sist         2.9  5 
DistRepRes26_interacao    4.6  8 
DistRepRes26_tempoSistema 4.9  6 

#rcp85
DistRepRes_0k85_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                                  tempo= rep(c("0k", "85"), each=29), 
                                  especie = spp_nearticos$V1,
                                  distancia=c(DistRepRes0k[dimnames(DistRepRes0k)[[1]] %in% spp_nearticos$V1,],
                                              DistRepRes85[dimnames(DistRepRes85)[[1]] %in% spp_nearticos$V1,]))

DistRepRes_0k85_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                                 tempo= rep(c("0k", "85"), each=40), 
                                 especie = spp_austrais$V1,
                                 distancia=c(DistRepRes0k[dimnames(DistRepRes0k)[[1]] %in% spp_austrais$V1,],
                                             DistRepRes85[dimnames(DistRepRes85)[[1]] %in% spp_austrais$V1,]))

DistRepRes_0k85_Intrat<-data.frame(sistema=rep("intrat", each=22),
                                   tempo= rep(c("0k", "85"), each=11), 
                                   especie = spp_intra$V1,
                                   distancia=c(DistRepRes0k[dimnames(DistRepRes0k)[[1]] %in% spp_intra$V1,],
                                               DistRepRes85[dimnames(DistRepRes85)[[1]] %in% spp_intra$V1,]))


DistRepRes_0k85_lmer<-rbind(DistRepRes_0k85_Aust,DistRepRes_0k85_Neart, DistRepRes_0k85_Intrat)
teste<-DistRepRes_0k85_lmer[complete.cases(DistRepRes_0k85_lmer),]
table(teste$tempo, teste$sistema)

aust_neot neart_neot intrat
0k        10         13      4
85         6          1      0

#DistRepRes_0k85_Intrat não possui áreas no futuro rcp85

library(lme4)
DistRepRes85_tempo<-lmer(distancia ~ tempo + (1|especie) ,data=DistRepRes_0k85_lmer,REML=F)
DistRepRes85_nulo<-lmer(distancia ~ 1 + (1|especie) ,data=DistRepRes_0k85_lmer,REML=F)
DistRepRes85_tempoSistema<-lmer(distancia ~ tempo + sistema + (1|especie) ,data=DistRepRes_0k85_lmer,REML=F)
DistRepRes85_sist<-lmer(distancia ~ sistema + (1|especie) ,data=DistRepRes_0k85_lmer,REML=F)
DistRepRes85_interacao<-lmer(distancia ~ tempo * sistema + (1|especie) ,data=DistRepRes_0k85_lmer,REML=F)
#fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

#AICtab(DistRepRes85_tempo,DistRepRes85_nulo,DistRepRes85_tempoSistema,DistRepRes85_sist,DistRepRes85_interacao) 


#INVERNADA - RESIDENCIA
library(lme4)
library(bbmle)
#RCP26 ********
DistInvRes_0k26_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                                  tempo= rep(c("0k", "26"), each=29), 
                                  especie = spp_nearticos$V1,
                                  distancia=c(DistInvRes0k[dimnames(DistInvRes0k)[[1]] %in% spp_nearticos$V1,],
                                              DistInvRes26[dimnames(DistInvRes26)[[1]] %in% spp_nearticos$V1,]))

DistInvRes_0k26_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                                 tempo= rep(c("0k", "26"), each=40), 
                                 especie = spp_austrais$V1,
                                 distancia=c(DistInvRes0k[dimnames(DistInvRes0k)[[1]] %in% spp_austrais$V1,],
                                             DistInvRes26[dimnames(DistInvRes26)[[1]] %in% spp_austrais$V1,]))

DistInvRes_0k26_Intrat<-data.frame(sistema=rep("intrat", each=22),
                                   tempo= rep(c("0k", "26"), each=11), 
                                   especie = spp_intra$V1,
                                   distancia=c(DistInvRes0k[dimnames(DistInvRes0k)[[1]] %in% spp_intra$V1,],
                                               DistInvRes26[dimnames(DistInvRes26)[[1]] %in% spp_intra$V1,]))


DistInvRes_0k26_lmer<-rbind(DistInvRes_0k26_Aust,DistInvRes_0k26_Neart, DistInvRes_0k26_Intrat)

#modelo misto
DistInvRes26_tempo<-lmer(distancia ~ tempo + (1|especie) ,data=DistInvRes_0k26_lmer,REML=F)
DistInvRes26_nulo<-lmer(distancia ~ 1 + (1|especie) ,data=DistInvRes_0k26_lmer,REML=F)
DistInvRes26_tempoSistema<-lmer(distancia ~ tempo + sistema + (1|especie) ,data=DistInvRes_0k26_lmer,REML=F)
DistInvRes26_sist<-lmer(distancia ~ sistema + (1|especie) ,data=DistInvRes_0k26_lmer,REML=F)
DistInvRes26_interacao<-lmer(distancia ~ tempo * sistema + (1|especie) ,data=DistInvRes_0k26_lmer,REML=F)

AICtab(DistInvRes26_tempo,DistInvRes26_nulo,DistInvRes26_tempoSistema,DistInvRes26_sist,DistInvRes26_interacao) 

dAIC df
DistInvRes26_tempo        0.0  4 
DistInvRes26_interacao    1.6  8 
DistInvRes26_tempoSistema 3.2  6 
DistInvRes26_nulo         6.4  3 
DistInvRes26_sist         9.6  5 

fixef(DistInvRes26_tempo)
(Intercept)     tempo26 
1319.904     226.975

#presente = 1319.904 
#futuro 1319.904 +  226.975 = 1546.879

#rcp85
DistInvRes_0k85_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                                  tempo= rep(c("0k", "85"), each=29), 
                                  especie = spp_nearticos$V1,
                                  distancia=c(DistInvRes0k[dimnames(DistInvRes0k)[[1]] %in% spp_nearticos$V1,],
                                              DistInvRes85[dimnames(DistInvRes85)[[1]] %in% spp_nearticos$V1,]))

DistInvRes_0k85_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                                 tempo= rep(c("0k", "85"), each=40), 
                                 especie = spp_austrais$V1,
                                 distancia=c(DistInvRes0k[dimnames(DistInvRes0k)[[1]] %in% spp_austrais$V1,],
                                             DistInvRes85[dimnames(DistInvRes85)[[1]] %in% spp_austrais$V1,]))

DistInvRes_0k85_Intrat<-data.frame(sistema=rep("intrat", each=22),
                                   tempo= rep(c("0k", "85"), each=11), 
                                   especie = spp_intra$V1,
                                   distancia=c(DistInvRes0k[dimnames(DistInvRes0k)[[1]] %in% spp_intra$V1,],
                                               DistInvRes85[dimnames(DistInvRes85)[[1]] %in% spp_intra$V1,]))


DistInvRes_0k85_lmer<-rbind(DistInvRes_0k85_Aust,DistInvRes_0k85_Neart, DistInvRes_0k85_Intrat)
teste<-DistInvRes_0k85_lmer[complete.cases(DistInvRes_0k85_lmer),]
table(teste$tempo, teste$sistema)
aust_neot neart_neot intrat
0k        25          8      8
85        12          0      3


#modelo misto
DistInvRes85_tempo<-lmer(distancia ~ tempo + (1|especie) ,data=DistInvRes_0k85_lmer,REML=F)
DistInvRes85_nulo<-lmer(distancia ~ 1 + (1|especie) ,data=DistInvRes_0k85_lmer,REML=F)
DistInvRes85_tempoSistema<-lmer(distancia ~ tempo + sistema + (1|especie) ,data=DistInvRes_0k85_lmer,REML=F)
DistInvRes85_sist<-lmer(distancia ~ sistema + (1|especie) ,data=DistInvRes_0k85_lmer,REML=F)
DistInvRes85_interacao<-lmer(distancia ~ tempo * sistema + (1|especie) ,data=DistInvRes_0k85_lmer,REML=F)


AICtab(DistInvRes85_tempo,DistInvRes85_nulo,DistInvRes85_tempoSistema,DistInvRes85_sist,DistInvRes85_interacao) 

dAIC df
DistInvRes85_tempo        0.0  4 
DistInvRes85_nulo         1.0  3 
DistInvRes85_interacao    2.4  8 
DistInvRes85_tempoSistema 3.2  6 
DistInvRes85_sist         4.6  5 

plot(DistInvRes85_tempo)
fixef(DistInvRes85_tempo)
(Intercept)     tempo85 
1322.9546    346.1966 


#quantas distancias foram calculadas 

#ENTRE INVERNADA E REPRODUÇÃO NO PRESENTE
sum(!is.na(DistInvRep_0k26_lmer$distancia[DistInvRep_0k26_lmer$tempo=="0k"])) # 48

#ENTRE INVERNADA E REPRODUÇÃO NO FUTURO RCP26
sum(!is.na(DistInvRep_0k26_lmer$distancia[DistInvRep_0k26_lmer$tempo=="26"])) # 40

#ENTRE INVERNADA E REPRODUÇÃO NO FUTURO RCP 85
sum(!is.na(DistInvRep_0k85_lmer$distancia[DistInvRep_0k85_lmer$tempo=="85"])) # 30



#ENTRE INVERNADA E RESIDÊNCIA NO PRESENTE
sum(!is.na(DistInvRes_0k85_lmer$distancia[DistInvRes_0k85_lmer$tempo=="0k"])) # 42

#ENTRE INVERNADA E RESIDÊNCIA NO FUTURO RCP26
sum(!is.na(DistInvRes_0k26_lmer$distancia[DistInvRes_0k26_lmer$tempo=="26"])) # 33

#ENTRE INVERNADA E RESIDÊNCIA NO FUTURO RCP 85
sum(!is.na(DistInvRes_0k85_lmer$distancia[DistInvRes_0k85_lmer$tempo=="85"])) # 16





#ENTRE RESIDÊNCIA E REPRODUÇÃO NO PRESENTE
sum(!is.na(DistRepRes_0k26_lmer$distancia[DistRepRes_0k26_lmer$tempo=="0k"])) # 28

#ENTRE RESIDÊNCIA E REPRODUÇÃO NO FUTURO RCP26
sum(!is.na(DistRepRes_0k26_lmer$distancia[DistRepRes_0k26_lmer$tempo=="26"])) # 18

#ENTRE RESIDÊNCIA E REPRODUÇÃO NO FUTURO RCP 85
sum(!is.na(DistRepRes_0k85_lmer$distancia[DistRepRes_0k85_lmer$tempo=="85"])) # 8






-------------------
  ##Maraísa Braga
  ##Leo Ré Jorge  
-------------------  


