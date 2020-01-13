----------------------------------
##ANÁLISE TAMANHO DAS ÁREAS##
----------------------------------  

nomes_index<-read.table("nomes_index.txt", h=T) #para a indexaçao 
resultados <- aves_migratorias_model_2$resultados_spp #RData originário do Modelo Bioclim

#listas das spp separadas por sistemas migratórios
spp_austrais<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_austrais.CSV", h=F)
spp_nearticos<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_nearticas.CSV", h=F)
spp_intra<-read.csv("D:/artigo aves migratorias/R/resultados analises/spp_intratropicais.csv", h=F)


#reproducao
Area0k.rep <- resultados[nomes_index$epoca=="0k" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",16,]
Area26.rep <- resultados[nomes_index$epoca=="26" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",16,]
Area85.rep <- resultados[nomes_index$epoca=="85" & nomes_index$fase=="REP" & nomes_index$modelo=="Ensemble",16,]


#invernada
Area0k.inv <- resultados[nomes_index$epoca=="0k" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",16,]
Area26.inv <- resultados[nomes_index$epoca=="26" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",16,]
Area85.inv <- resultados[nomes_index$epoca=="85" & nomes_index$fase=="INV" & nomes_index$modelo=="Ensemble",16,]



#residencia
Area0k.res <- resultados[nomes_index$epoca=="0k" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",16,]
Area26.res <- resultados[nomes_index$epoca=="26" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",16,]
Area85.res <- resultados[nomes_index$epoca=="85" & nomes_index$fase=="RES" & nomes_index$modelo=="Ensemble",16,]

#tyrannus melacholicus não tem área de invernada
resultados[nomes_index$fase=="INV",,74] <- NA

range(Area0k.res, na.rm=T)
range(Area0k.rep, na.rm=T)
range(Area0k.inv, na.rm=T)




## Diferenças entre as areas entre o presente e o futuro e entre os hemisferios
##ANÁLISE UTILIZANDO UMA SELEÇÃO DE MODELOS MISTOS


#CRIAR OS OBJETOS PARA A SELEÇÃO E SEPARAR POR SISTEMAS MIGRATÓRIOS
####REPRODUCAO

####rcp26
areaRep26_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                            epoca= rep(c("0k", "26"), each=29), 
                            especie = spp_nearticos$V1,
                            area=c(Area0k.rep[names(Area0k.rep) %in% spp_nearticos$V1],
                                   Area26.rep[names(Area26.rep) %in% spp_nearticos$V1]))


areaRep26_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                           epoca= rep(c("0k", "26"), each=40) , 
                           especie=spp_austrais$V1,
                           area=c(Area0k.rep[names(Area0k.rep) %in% spp_austrais$V1],
                                  Area26.rep[names(Area26.rep) %in% spp_austrais$V1]))

areaRep26_Intrat<-data.frame(sistema=rep("intratropical", each=22),
                             epoca= rep(c("0k", "26"), each=11) , 
                             especie=spp_intra$V1,
                             area=c(Area0k.rep[names(Area0k.rep) %in% spp_intra$V1],
                                    Area26.rep[names(Area26.rep) %in% spp_intra$V1]))

areaRep26_lmer<-rbind(areaRep26_Aust,areaRep26_Neart, areaRep26_Intrat)

teste<-areaRep26_lmer[complete.cases(areaRep26_lmer),]
table(teste$sistema , teste$epoca)




library(lme4)
library("bbmle")
Rep26_lmermodnulo<-lmer(area ~ 1 + (1|especie)  , data=areaRep26_lmer, REML=F)
Rep26_lmermodtempo<-lmer(area ~ epoca + (1|especie) , data=areaRep26_lmer,REML=F)
Rep26_lmermodtemposist<-lmer(area ~ epoca + sistema + (1|especie) , data=areaRep26_lmer,REML=F)
Rep26_lmermodsist<-lmer(area ~ sistema + (1|especie) , data=areaRep26_lmer,REML=F)
Rep26_lmermodinteracao<-lmer(area ~ epoca * sistema + (1|especie) , data=areaRep26_lmer,REML=F)

AICtab(Rep26_lmermodtempo,Rep26_lmermodnulo,Rep26_lmermodtemposist,Rep26_lmermodsist,Rep26_lmermodinteracao) 

dAIC df
Rep26_lmermodtemposist  0.0 6 
Rep26_lmermodinteracao  2.6 8 
Rep26_lmermodtempo      8.8 4 
Rep26_lmermodsist      30.6 5 
Rep26_lmermodnulo      39.4 3 


fixef(Rep26_lmermodtemposist)
(Intercept)              epoca26    sistemaneart_neot sistemaintratropical 
1154537.0            -295464.7              3022957.2            -624511.0 
#o intercepto é a area no presente para o sistema austral = 1.154.537,0
#sistema austral no futuro 1154537.0 - 295464.7 = 859.072,3
# 100-((859072.3*100)/1154537.0) perda de 25.6%

#sistema neártico no presente 1154537.0 + 3022957.2 =  4.177.494
#sistema neártico no futuro 1154537.0 + 3022957.2 -295464.7 =   3.882.030
100-(( 3882030*100)/4177494) #perda de 7%

#sistema intratropical no presente 1154537.0 - 624511.0  =  530026
#sistema intratropical no futuro 1154537.0 -624511.0  -295464.7 =   234561.3
 100-(( 234561.3*100)/423495.4) #perda de 44.6%


#para analisar o modelo "sistema" qual a diferença entre os sistemas migratórios
library(emmeans)
emmeans(Rep26_lmermodtemposist, pairwise ~ sistema, adjust = "Tukey")
$emmeans
sistema        emmean      SE  df asymp.LCL asymp.UCL
aust_neot     1006805  655050 Inf   -277070   2290679
neart_neot    4029762  583364 Inf   2886390   5173134
intratropical  382294 1404926 Inf  -2371311   3135898

Results are averaged over the levels of: epoca 
Degrees-of-freedom method: asymptotic 
Confidence level used: 0.95 

$contrasts
contrast                   estimate      SE  df z.ratio p.value
aust_neot - neart_neot     -3022957  877157 Inf -3.446  0.0016 
aust_neot - intratropical    624511 1550132 Inf  0.403  0.9144 
neart_neot - intratropical  3647468 1521227 Inf  2.398  0.0435 

Results are averaged over the levels of: epoca 
P value adjustment: tukey method for comparing a family of 3 estimates 


####rcp85
areaRep85_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                            epoca= rep(c("0k", "85"), each=29) , 
                            especie=spp_nearticos$V1,
                            area=c(Area0k.rep[names(Area0k.rep) %in% spp_nearticos$V1],
                                   Area85.rep[names(Area85.rep) %in% spp_nearticos$V1]))


areaRep85_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                           epoca= rep(c("0k", "85"), each=40) , 
                           especie=spp_austrais$V1,
                           area=c(Area0k.rep[names(Area0k.rep) %in% spp_austrais$V1],
                                  Area85.rep[names(Area85.rep) %in% spp_austrais$V1]))

areaRep85_Intrat<-data.frame(sistema=rep("intratropical", each=22),
                             epoca= rep(c("0k", "85"), each=11) , 
                             especie = spp_intra$V1,
                             area=c(Area0k.rep[names(Area0k.rep) %in% spp_intra$V1],
                                    Area85.rep[names(Area85.rep) %in% spp_intra$V1]))

areaRep85_lmer<-rbind(areaRep85_Aust,areaRep85_Neart, areaRep85_Intrat)



teste<-areaRep85_lmer[complete.cases(areaRep85_lmer),]

#quantas spp calculamos a área no presente
teste<-teste[teste$epoca=="0k",]
table(teste$sistema)
aust_neot    neart_neot intratropical 
      23            29             5  =57

      sum(complete.cases(areaRep85_Aust[areaRep85_Aust$epoca=="0k",])) # =  23
      sum(areaRep85_Aust[,4]==0, na.rm=T) # = 6
      sum(complete.cases(areaRep85_Neart[areaRep85_Neart$epoca=="0k",])) # = 29
      sum(areaRep85_Neart[,4]==0, na.rm=T) # = 9
      sum(complete.cases(areaRep85_Intrat[areaRep85_Intrat$epoca=="0k",])) # = 5
      sum(areaRep85_Intrat[,4]==0, na.rm=T) # = 3
      
      
##quantas spp perderam toda sua área
teste<-teste[teste$epoca=="85" & teste$area=="0",]
table(teste$sistema)
aust_neot    neart_neot intratropical 
6             9             3  = 18

1800/57 = 31%

library(lme4)
library("bbmle")
#MODELO MISTO
Rep85_lmermodtempo<-lmer(area ~ epoca + (1|especie) , data=areaRep85_lmer,REML=F)
Rep85_lmermodnulo<-lmer(area ~ 1 + (1|especie), data=areaRep85_lmer,REML=F)
Rep85_lmermodtemposist<-lmer(area ~ epoca + sistema + (1|especie) , data=areaRep85_lmer,REML=F)
Rep85_lmermodsist<-lmer(area ~ sistema + (1|especie) , data=areaRep85_lmer,REML=F)
Rep85_lmermodinteracao<-lmer(area ~ epoca * sistema + (1|especie) , data=areaRep85_lmer,REML=F)

AICtab(Rep85_lmermodtempo,Rep85_lmermodnulo,Rep85_lmermodtemposist,Rep85_lmermodsist,Rep85_lmermodinteracao) 

dAIC df
Rep85_lmermodinteracao  0.0 8 
Rep85_lmermodtemposist  1.5 6 
Rep85_lmermodtempo     10.2 4 
Rep85_lmermodsist      34.9 5 
Rep85_lmermodnulo      43.6 3 

plot(Rep85_lmermodinteracao)
fixef(Rep85_lmermodinteracao)
(Intercept)                      epoca85            sistemaneart_neot 
1176317.5                    -785348.0                    2996471.9 
sistemaintratropical    epoca85:sistemaneart_neot epoca85:sistemaintratropical 
-719194.5                    -548742.8                     443011.5 

  #o intercepto é a area no presente para o sistema austral = 1154537.0
  #sistema austral no futuro 1176317.5 - 785348.0   = 390969.5
  # 100-(( 390969.5*100)/1154537.0) perda de 66%
1154537.0 -  390969.5  =  763567.5 km²

  #sistema neártico no presente 1176317.5 + 2996471.9  =   4.172.789
  #sistema neártico no futuro 1176317.5 + 2996471.9  -785348.0 -548742.8 =  2.838.699
  # 100-((2838699*100)/  4172789) perda de 31%
4172789-2838699 =  1.334.090
  
  #sistema intratropical no presente 1176317.5 - 719194.5  =  457123
  #sistema intratropical no futuro 1176317.5 - 719194.5   -785348.0 + 443011.5 =  114786.5
  # 100-((114786.5*100)/ 457123) perda de 75%
  
  
#para analisar no modelo "interação" qual a diferença entre os sistemas migratórios
library(emmeans)
emmeans(Rep85_lmermodinteracao, pairwise ~ sistema + epoca , adjust = "Tukey")
$emmeans
sistema       epoca  emmean      SE  df asymp.LCL asymp.UCL
aust_neot     0k    1176318  600456 Inf      -554   2353190
neart_neot    0k    4172789  534744 Inf   3124710   5220869
intratropical 0k     457123 1287835 Inf  -2066986   2981232
aust_neot     85     390970  600456 Inf   -785902   1567842
neart_neot    85    2838699  534744 Inf   1790619   3886778
intratropical 85     114787 1287835 Inf  -2409323   2638896

Degrees-of-freedom method: asymptotic 
Confidence level used: 0.95 

$contrasts
contrast                            estimate      SE  df z.ratio p.value
aust_neot,0k - neart_neot,0k        -2996472  804051 Inf -3.727  0.0027 
aust_neot,0k - intratropical,0k       719194 1420938 Inf  0.506  0.9960 
aust_neot,0k - aust_neot,85           785348  219504 Inf  3.578  0.0047 
aust_neot,0k - neart_neot,85        -1662381  804051 Inf -2.068  0.3044 
aust_neot,0k - intratropical,85      1061531 1420938 Inf  0.747  0.9760 
neart_neot,0k - intratropical,0k     3715666 1394442 Inf  2.665  0.0825 
neart_neot,0k - aust_neot,85         3781820  804051 Inf  4.703  <.0001 
neart_neot,0k - neart_neot,85        1334091  195483 Inf  6.825  <.0001 
neart_neot,0k - intratropical,85     4058003 1394442 Inf  2.910  0.0421 
intratropical,0k - aust_neot,85        66154 1420938 Inf  0.047  1.0000 
intratropical,0k - neart_neot,85    -2381576 1394442 Inf -1.708  0.5264 
intratropical,0k - intratropical,85   342336  470784 Inf  0.727  0.9787 
aust_neot,85 - neart_neot,85        -2447729  804051 Inf -3.044  0.0282 
aust_neot,85 - intratropical,85       276183 1420938 Inf  0.194  1.0000 
neart_neot,85 - intratropical,85     2723912 1394442 Inf  1.953  0.3695 

P value adjustment: tukey method for comparing a family of 6 estimates 


####INVERNADA##
####RCP26
areaInv26_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                            epoca= rep(c("0k", "26"), each=29) , 
                            especie=spp_nearticos$V1,
                            area=c(Area0k.inv[names(Area0k.inv) %in% spp_nearticos$V1],
                                   Area26.inv[names(Area26.inv) %in% spp_nearticos$V1]))


areaInv26_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                            epoca= rep(c("0k", "26"), each=40) , 
                            especie=spp_austrais$V1,
                            area=c(Area0k.inv[names(Area0k.inv) %in% spp_austrais$V1],
                                   Area26.inv[names(Area26.inv) %in% spp_austrais$V1]))

areaInv26_Intrat<-data.frame(sistema=rep("intratropical", each=22),
                           epoca= rep(c("0k", "26"), each=11) , 
                           especie=spp_intra$V1,
                           area=c(Area0k.inv[names(Area0k.inv) %in% spp_intra$V1],
                                  Area26.inv[names(Area26.inv) %in% spp_intra$V1]))

areaInv26_lmer<-rbind(areaInv26_Aust,areaInv26_Neart, areaInv26_Intrat)

Inv26_lmermodnulo<-lmer(area ~ 1 + (1|especie), data=areaInv26_lmer,REML=F)
Inv26_lmermodtempo<-lmer(area ~ epoca + (1|especie) , data=areaInv26_lmer,REML=F)
Inv26_lmermodtemposist<-lmer(area ~ epoca + sistema + (1|especie) , data=areaInv26_lmer,REML=F)
Inv26_lmermodsist<-lmer(area ~ sistema + (1|especie) , data=areaInv26_lmer,REML=F)
Inv26_lmermodinteracao<-lmer(area ~ epoca * sistema + (1|especie) , data=areaInv26_lmer,REML=F)

AICtab(Inv26_lmermodtempo,Inv26_lmermodnulo,Inv26_lmermodtemposist,Inv26_lmermodsist,Inv26_lmermodinteracao) 

dAIC df
Inv26_lmermodtempo      0.0 4 
Inv26_lmermodtemposist  3.6 6 
Inv26_lmermodinteracao  6.3 8 
Inv26_lmermodnulo      37.2 3 
Inv26_lmermodsist      40.9 5 

fixef(Inv26_lmermodtempo)
(Intercept)     epoca26 
2159451.3   -965807.7 

#área no presente para todos os sistemas = 2.159.451,3   km²
#área no futuro para todos os sistemas 2159451.3   -965807.7 = 1.193.644
100-(( 1193644*100)/ 2159451.3) #perda de 44%
> 2159451-1193644
[1] 965.807

####rcp85
areaInv85_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                            epoca= rep(c("0k", "85"), each=29) , 
                            especie=spp_nearticos$V1,
                            area=c(Area0k.inv[names(Area0k.inv) %in% spp_nearticos$V1],
                                   Area85.inv[names(Area85.inv) %in% spp_nearticos$V1]))


areaInv85_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                           epoca= rep(c("0k", "85"), each=40) , 
                           especie=spp_austrais$V1,
                           area=c(Area0k.inv[names(Area0k.inv) %in% spp_austrais$V1],
                                  Area85.inv[names(Area85.inv) %in% spp_austrais$V1]))

areaInv85_Intrat<-data.frame(sistema=rep("intratropical", each=22),
                             epoca= rep(c("0k", "85"), each=11) , 
                             especie=spp_intra$V1,
                             area=c(Area0k.inv[names(Area0k.inv) %in% spp_intra$V1],
                                    Area85.inv[names(Area85.inv) %in% spp_intra$V1]))

areaInv85_lmer<-rbind(areaInv85_Aust,areaInv85_Neart, areaInv85_Intrat)

teste<-areaInv85_lmer[complete.cases(areaInv85_lmer),]
sum(complete.cases(areaInv85_Aust[areaInv85_Aust$epoca=="0k",])) # =  38
sum(areaInv85_Aust[,4]==0, na.rm=T) # = 13
sum(complete.cases(areaInv85_Neart[areaInv85_Neart$epoca=="0k",])) # = 24
sum(areaInv85_Neart[,4]==0, na.rm=T) # = 9
sum(complete.cases(areaInv85_Intrat[areaInv85_Intrat$epoca=="0k",])) # = 9
sum(areaInv85_Intrat[,4]==0, na.rm=T) # = 2


  #quantas spp calculamos a área no presente
  teste<-teste[teste$epoca=="0k",]
  table(teste$sistema)
  aust_neot    neart_neot intratropical 
  38            24             9 = 71
  
  #quantas spp zeraram as áreas
  teste<-teste[teste$epoca=="85" & teste$area=="0",]
  table(teste$sistema)
  aust_neot    neart_neot intratropical 
  13             9             2 = 24
   
  (24*100)/71  =  33%


Inv85_lmermodnulo<-lmer(area ~ 1 + (1|especie), data=areaInv85_lmer,REML=F)
Inv85_lmermodtempo<-lmer(area ~ epoca + (1|especie) , data=areaInv85_lmer,,REML=F)
Inv85_lmermodtemposist<-lmer(area ~ epoca + sistema + (1|especie) , data=areaInv85_lmer,REML=F)
Inv85_lmermodsist<-lmer(area ~ sistema + (1|especie) , data=areaInv85_lmer,REML=F)
Inv85_lmermodinteracao<-lmer(area ~ epoca * sistema + (1|especie) , data=areaInv85_lmer,REML=F)

AICtab(Inv85_lmermodtempo,Inv85_lmermodnulo,Inv85_lmermodtemposist,Inv85_lmermodsist,Inv85_lmermodinteracao)

dAIC df
Inv85_lmermodtempo      0.0 4 
Inv85_lmermodtemposist  3.6 6 
Inv85_lmermodinteracao  7.3 8 
Inv85_lmermodnulo      33.0 3 
Inv85_lmermodsist      36.7 5 

plot(Inv85_lmermodtempo)
fixef(Inv85_lmermodtempo)
(Intercept)     epoca85 
2159451    -1877197 

#área no presente para todos os sistemas = 2159451   km²
#área no futuro para todos os sistemas 2159451 - 1877197 = 282254 km²
100-(( 282254*100)/ 2159451) #perda de 86%
2159451-282254 = 1.877.197

#### RESIDENCIA##
####RCP26
areaRes26_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                            epoca= rep(c("0k", "26"), each=29) , 
                            especie=spp_nearticos$V1,
                            area=c(Area0k.res[names(Area0k.res) %in% spp_nearticos$V1],
                                   Area26.res[names(Area26.res) %in% spp_nearticos$V1]))


areaRes26_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                           epoca= rep(c("0k", "26"), each=40) , 
                           especie=spp_austrais$V1,
                           area=c(Area0k.res[names(Area0k.res) %in% spp_austrais$V1],
                                  Area26.res[names(Area26.res) %in% spp_austrais$V1]))

areaRes26_Intrat<-data.frame(sistema=rep("intratropical", each=22),
                             epoca= rep(c("0k", "26"), each=11) , 
                             especie=spp_intra$V1,
                             area=c(Area0k.res[names(Area0k.res) %in% spp_intra$V1],
                                    Area26.res[names(Area26.res) %in% spp_intra$V1]))

areaRes26_lmer<-rbind(areaRes26_Aust,areaRes26_Neart, areaRes26_Intrat)

#modelo misto
Res26_lmermodtempo<-lmer(area ~ epoca + (1|especie) , data=areaRes26_lmer,REML=F)
Res26_lmermodnulo<-lmer(area ~ 1 + (1|especie), data=areaRes26_lmer,REML=F)
Res26_lmermodtemposist<-lmer(area ~ epoca + sistema + (1|especie) , data=areaRes26_lmer,REML=F)
Res26_lmermodsist<-lmer(area ~ sistema + (1|especie) , data=areaRes26_lmer,REML=F)
Res26_lmermodinteracao<-lmer(area ~ epoca * sistema + (1|especie) , data=areaRes26_lmer,REML=F)

AICtab(Res26_lmermodtempo,Res26_lmermodnulo,Res26_lmermodtemposist,Res26_lmermodsist,Res26_lmermodinteracao) 

dAIC df
Res26_lmermodtempo      0.0 4 
Res26_lmermodinteracao  3.5 8 
Res26_lmermodtemposist  3.8 6 
Res26_lmermodnulo      22.7 3 
Res26_lmermodsist      26.5 5 

fixef(Res26_lmermodtempo)
(Intercept)     epoca26 
3529723    -1298340  

#área no presente para todos os sistemas = 3529723   km²
#área no futuro para todos os sistemas 3529723 - 1298340 =  2231383 km²
100-((  2231383*100)/ 3529723) #perda de 36.8%
3529723-2231383 = 1.298.340

####rcp85
areaRes85_Neart<-data.frame(sistema=rep("neart_neot", each=58),
                            epoca= rep(c("0k", "85"), each=29) , 
                            especie=spp_nearticos$V1,
                            area=c(Area0k.res[names(Area0k.res) %in% spp_nearticos$V1],
                                   Area85.res[names(Area85.res) %in% spp_nearticos$V1]))


areaRes85_Aust<-data.frame(sistema=rep("aust_neot", each=80),
                           epoca= rep(c("0k", "85"), each=40) , 
                           especie=spp_austrais$V1,
                           area=c(Area0k.res[names(Area0k.res) %in% spp_austrais$V1],
                                  Area85.res[names(Area85.res) %in% spp_austrais$V1]))

areaRes85_Intrat<-data.frame(sistema=rep("intratropical", each=22),
                             epoca= rep(c("0k", "85"), each=11) , 
                             especie=spp_intra$V1,
                             area=c(Area0k.res[names(Area0k.res) %in% spp_intra$V1],
                                    Area85.res[names(Area85.res) %in% spp_intra$V1]))

areaRes85_lmer<-rbind(areaRes85_Aust,areaRes85_Neart, areaRes85_Intrat)

teste<-areaRes85_lmer[complete.cases(areaRes85_lmer),]

#quantas spp calculamos a área no presente
teste<-teste[teste$epoca=="0k",]
table(teste$sistema)
aust_neot    neart_neot intratropical 
27            14            10  = total = 51


#quantas spp zeraram sua área?
teste<-teste[teste$epoca=="85" & teste$area=="0",]
table(teste$sistema)
aust_neot    neart_neot intratropical 
2             6             4  = total = 12

(12*100)/51 = 23 %

#modelo misto
Res85_lmermodtempo<-lmer(area ~ epoca + (1|especie) , data=areaRes85_lmer,REML=F)
Res85_lmermodnulo<-lmer(area ~ 1 + (1|especie), data=areaRes85_lmer,REML=F)
Res85_lmermodtemposist<-lmer(area ~ epoca + sistema + (1|especie) , data=areaRes85_lmer,REML=F)
Res85_lmermodsist<-lmer(area ~ sistema + (1|especie) , data=areaRes85_lmer,REML=F)
Res85_lmermodinteracao<-lmer(area ~ epoca * sistema + (1|especie) , data=areaRes85_lmer,REML=F)

AICtab(Res85_lmermodtempo,Res85_lmermodnulo,Res85_lmermodtemposist,Res85_lmermodsist,Res85_lmermodinteracao) 
dAIC df
Res85_lmermodtempo      0.0 4 
Res85_lmermodtemposist  3.7 6 
Res85_lmermodinteracao  6.7 8 
Res85_lmermodnulo      25.8 3 
Res85_lmermodsist      29.5 5 

plot(Res85_lmermodtempo)
fixef(Res85_lmermodtempo)
(Intercept)     epoca85 
3529723    -2797105 

#área no presente para todos os sistemas = 3529723   km²
#área no futuro para todos os sistemas 3529723 - 2797105  =  732618 km²
100-((732618*100)/ 3529723) #perda de 79%

3529723 - 732618 # =  2.797.105








####Graficos
####Cenario de baixa emissão

####reproducao
##baixa emissão
min(areaRep26_lmer$area, na.rm=T) #para determinar o ylim 0
max(areaRep26_lmer$area, na.rm=T) #para determinar o ylim 14 489 236

tiff("areas_0k26REP.tiff", width = 4200, height = 2200, compression = "lzw", res=300 )
options(scipen=10)
par(mar=c(6,6,5,5))
boxplot(area ~ epoca + sistema, data=areaRep26_lmer, border="white",names=rep(c("presente","futuro"),times=3), par(bty="l"), ylim=c(0,16000000), cex.lab=1.2, ylab="Tamanho das áreas (Km²)", xlab="Austral-Neotropical                                        Neártico-Neotropical                                        Intratropical")


setas_aust<-data.frame(especie=areaRep26_Aust[1:41,3], presente=areaRep26_Aust[1:41,4], futuro=areaRep26_Aust[42:82,4])


arrows(x0=1,y0=setas_aust$presente,
       x1=2,y1=setas_aust$futuro, 
       length=0, col="azure4",lwd=2)

setas_neart<-data.frame(especie=areaRep26_Neart[1:27,3], presente=areaRep26_Neart[1:27,4], futuro=areaRep26_Neart[28:54,4])

arrows(x0=3,y0=setas_neart$presente,
       x1=4,y1=setas_neart$futuro, 
       length=0, col="azure4",lwd=2)

setas_intrat<-data.frame(especie=areaRep26_Intrat[1:12,3], presente=areaRep26_Intrat[1:12,4], futuro=areaRep26_Intrat[13:24,4])
arrows(x0=5,y0=setas_intrat$presente,
       x1=6,y1=setas_intrat$futuro, 
       length=0, col="azure4",lwd=2)

library(lme4)
#média para areas de reprodução segundo o modelo mais explicativo
fixef(Rep26_lmermodtemposist)
(Intercept)              epoca26    sistemaneart_neot sistemaintratropical 
1154537.0            -295464.7            3320759.1            -731041.6 
#*LEMBRETE: Como o modelo nÃo inclui interação a diferença no presente e no fututo se mantem a mesma

#area média no presente sp austrais 1 154 537.0 
#area media no futuro 26 spp austrais 1154537.0 - 295464.7 = 859 072.3
#area média no presente dos sitema neartico 3320759.1  + 1154537.0 = 4475296
#area média no futuro dos sitema neartico 3320759.1  + 1154537.0 -295464.7 = 4179831
#area média no presente dos sitema intratropical -731041.6  + 1154537.0  =   423495.4
#area média no futuro dos sitema intratropical -731041.6  + 1154537.0 -295464.7 =  128030.7

arrows(x0=1,y0=1154537.0, x1=2,y1=859072.3, col="red", length=0, lwd=4) #média sistema austral
100-(859072.3*100/1154537.0 ) #25.6% de perda de área para as spp austrais 

arrows(x0=3,y0=4475296, x1=4,y1= 4179831, col="red", length=0, lwd=4) #média sistema neartico
100-(4179831*100/4475296) #6.6% de perda de area para spp nearticos

arrows(x0=5,y0=423495.4, x=6,y1= 128030.7, col="red", length=0, lwd=4) #média sistema intratroical
100-(128030.7*100/423495.4) # 69.8% de perda de area para spp nearticos


dev.off() 



#ALTA EMISSÃO

min(areaRep85_lmer$area, na.rm=T) #para determinar o ylim 0
max(areaRep85_lmer$area, na.rm=T) #para determinar o ylim  14 489 236

tiff("areas_0k85REP.tiff", width = 4200, height = 2200, compression = "lzw", res=300 )
options(scipen=10)
par(mar=c(6,6,5,5))
boxplot(area ~ epoca + sistema, data=areaRep85_lmer, border="white",names=rep(c("presente","futuro"),times=3), par(bty="l"), ylim=c(0,16000000), cex.lab=1.2, ylab="Tamanho das áreas (Km²)", xlab="Austral-Neotropical                                        Neártico-Neotropical                                        Intratropical")


setas_aust<-data.frame(especie=areaRep85_Aust[1:41,3], presente=areaRep85_Aust[1:41,4], futuro=areaRep85_Aust[42:82,4])


arrows(x0=1,y0=setas_aust$presente,
       x1=2,y1=setas_aust$futuro, 
       length=0, col="azure4",lwd=2)

setas_neart<-data.frame(especie=areaRep85_Neart[1:27,3], presente=areaRep85_Neart[1:27,4], futuro=areaRep85_Neart[28:54,4])

arrows(x0=3,y0=setas_neart$presente,
       x1=4,y1=setas_neart$futuro, 
       length=0, col="azure4",lwd=2)

setas_intrat<-data.frame(especie=areaRep85_Intrat[1:12,3], presente=areaRep85_Intrat[1:12,4], futuro=areaRep85_Intrat[13:24,4])
arrows(x0=5,y0=setas_intrat$presente,
       x1=6,y1=setas_intrat$futuro, 
       length=0, col="azure4",lwd=2)

library(lme4)
fixef(Rep85_lmermodinteracao)
(Intercept)                      epoca85            sistemaneart_neot         sistemaintratropical 
1176317.5                    -785348.0                    3304169.4                    -844408.7 
epoca85:sistemaneart_neot epoca85:sistemaintratropical 
-646166.3                     535429.5 

#LEMBRETE: Sempre somar a interação se olha para um sistema no futuro

#area média no presente sp austrais 1176317.54
#area media no futuro 85 spp austrais 1176317.54 -785348.0  = 390969.5

#area média no presente do sitema neartico 1176317.54 + 3304169.40  =  4480487
#area média no futuro dos sitema neartico    1176317.54 -785348.0 + 3304169.40 -646166.3 = 3048973

#area média no presente dos sitema intratropical 1176317.54 - 844408.66   =   331908.9
#area média no futuro dos sitema intratropical 1176317.54 - 844408.66 -785348.0 +  535429.5=  81990.38


emmeans(Rep85_lmermodinteracao, pairwise ~ sistema, adjust = "Tukey")

arrows(x0=1,y0=1176317.54, x1=2,y1=390969.5, col="red", length=0, lwd=4) #média sistema austral
100-(390969.5*100/1176317.5  ) #67% de perda de área para as spp austrais 

arrows(x0=3,y0=4480487, x1=4,y1= 3048973, col="red", length=0, lwd=4) #média sistema neartico
100-(3048973*100/4480487) #31.9% de perda de area para spp nearticos

arrows(x0=5,y0=331908.9, x=6,y1= 81990.38, col="red", length=0, lwd=4) #média sistema intratroical
100-(81990.38*100/331908.9) # 75.3% de perda de area para spp nearticos


dev.off() 

####Invernada 
##baixa emissão

min(areaInv26_lmer$area, na.rm=T) #para determinar o ylim 0
max(areaInv26_lmer$area, na.rm=T) #para determinar o ylim  15390230

tiff("areas_0k26INV.tiff", width = 4200, height = 2200, compression = "lzw", res=300 )
options(scipen=10)
par(mar=c(6,6,5,5))
boxplot(area ~ epoca + sistema, data=areaInv26_lmer, border="white",names=rep(c("presente","futuro"),times=3), par(bty="l"), ylim=c(0,16000000), cex.lab=1.2, ylab="Tamanho das áreas (Km²)", xlab="Austral-Neotropical                                        Neártico-Neotropical                                        Intratropical")


setas_aust<-data.frame(especie=areaInv26_Aust[1:41,3], presente=areaInv26_Aust[1:41,4], futuro=areaInv26_Aust[42:82,4])


arrows(x0=1,y0=setas_aust$presente,
       x1=2,y1=setas_aust$futuro, 
       length=0, col="azure4",lwd=2)

setas_neart<-data.frame(especie=areaInv26_Neart[1:27,3], presente=areaInv26_Neart[1:27,4], futuro=areaInv26_Neart[28:54,4])

arrows(x0=3,y0=setas_neart$presente,
       x1=4,y1=setas_neart$futuro, 
       length=0, col="azure4",lwd=2)

setas_intrat<-data.frame(especie=areaInv26_Intrat[1:12,3], presente=areaInv26_Intrat[1:12,4], futuro=areaInv26_Intrat[13:24,4])
arrows(x0=5,y0=setas_intrat$presente,
       x1=6,y1=setas_intrat$futuro, 
       length=0, col="azure4",lwd=2)

library(lme4)
fixef(Inv26_lmermodtempo)
(Intercept)     epoca26 
2159451.3   -965807.7 
  
#area média no presente para todas as spp 2159451.3 
#area media no futuro 26 para todas as spp 2159451.3  -965807.7 = 1193644


arrows(x0=1,y0=2159451.3, x1=2,y1=1193644, col="red", length=0, lwd=4) #média sistema austral

arrows(x0=3,y0=2159451.3, x1=4,y1= 1193644, col="red", length=0, lwd=4) #média sistema neartico

arrows(x0=5,y0=2159451.3, x=6,y1= 1193644, col="red", length=0, lwd=4) #média sistema intratroical
100-(1193644*100/2159451.3) # 44.7% de perda de area para todas as spp


dev.off() 


##alta emissão

min(areaInv85_lmer$area, na.rm=T) #para determinar o ylim 0
max(areaInv85_lmer$area, na.rm=T) #para determinar o ylim  15390230

tiff("areas_0k85INV.tiff", width = 4200, height = 2200, compression = "lzw", res=300 )
options(scipen=10)
par(mar=c(6,6,5,5))
boxplot(area ~ epoca + sistema, data=areaInv85_lmer, border="white",names=rep(c("presente","futuro"),times=3), par(bty="l"), ylim=c(0,16000000), cex.lab=1.2, ylab="Tamanho das áreas (Km²)", xlab="Austral-Neotropical                                        Neártico-Neotropical                                        Intratropical")


setas_aust<-data.frame(especie=areaInv85_Aust[1:41,3], presente=areaInv85_Aust[1:41,4], futuro=areaInv85_Aust[42:82,4])


arrows(x0=1,y0=setas_aust$presente,
       x1=2,y1=setas_aust$futuro, 
       length=0, col="azure4",lwd=2)

setas_neart<-data.frame(especie=areaInv85_Neart[1:27,3], presente=areaInv85_Neart[1:27,4], futuro=areaInv85_Neart[28:54,4])

arrows(x0=3,y0=setas_neart$presente,
       x1=4,y1=setas_neart$futuro, 
       length=0, col="azure4",lwd=2)

setas_intrat<-data.frame(especie=areaInv85_Intrat[1:12,3], presente=areaInv85_Intrat[1:12,4], futuro=areaInv85_Intrat[13:24,4])
arrows(x0=5,y0=setas_intrat$presente,
       x1=6,y1=setas_intrat$futuro, 
       length=0, col="azure4",lwd=2)

library(lme4)
fixef(Inv85_lmermodtempo)
(Intercept)     epoca85 
2159451    -1877197 
#area média no presente para todas as spp 2159451
#area media no futuro 26 para todas as spp 2159451 - 1877197  = 282254


arrows(x0=1,y0=2159451, x1=2,y1=282254, col="red", length=0, lwd=4) #média sistema austral

arrows(x0=3,y0=2159451, x1=4,y1= 282254, col="red", length=0, lwd=4) #média sistema neartico

arrows(x0=5,y0=2159451, x=6,y1= 282254, col="red", length=0, lwd=4) #média sistema intratroical
100-(282254*100/2159451) # 86.9% de perda de area para todas as spp


dev.off()


####residencia
##baixa emissão

min(areaRes26_lmer$area, na.rm=T) #para determinar o ylim 0
max(areaRes26_lmer$area, na.rm=T) #para determinar o ylim  15619260

tiff("areas_0k26RES.tiff", width = 4200, height = 2200, compression = "lzw", res=300 )
options(scipen=10)
par(mar=c(6,6,5,5))
boxplot(area ~ epoca + sistema, data=areaRes26_lmer, border="white",names=rep(c("presente","futuro"),times=3), par(bty="l"), ylim=c(0,16000000), cex.lab=1.2, ylab="Tamanho das áreas (Km²)", xlab="Austral-Neotropical                                        Neártico-Neotropical                                        Intratropical")


setas_aust<-data.frame(especie=areaRes26_Aust[1:41,3], presente=areaRes26_Aust[1:41,4], futuro=areaRes26_Aust[42:82,4])


arrows(x0=1,y0=setas_aust$presente,
       x1=2,y1=setas_aust$futuro, 
       length=0, col="azure4",lwd=2)

setas_neart<-data.frame(especie=areaRes26_Neart[1:27,3], presente=areaRes26_Neart[1:27,4], futuro=areaRes26_Neart[28:54,4])

arrows(x0=3,y0=setas_neart$presente,
       x1=4,y1=setas_neart$futuro, 
       length=0, col="azure4",lwd=2)

setas_intrat<-data.frame(especie=areaRes26_Intrat[1:12,3], presente=areaRes26_Intrat[1:12,4], futuro=areaRes26_Intrat[13:24,4])
arrows(x0=5,y0=setas_intrat$presente,
       x1=6,y1=setas_intrat$futuro, 
       length=0, col="azure4",lwd=2)

library(lme4)
fixef(Res85_lmermodtempo)
(Intercept)     epoca85 
3529723    -2797105  

#area média no presente para todas as spp 3529723 
#area media no futuro 26 para todas as spp 3529723    -2797105  = 732618


-------------------
  ##Maraísa Braga
  ##Leo Ré Jorge  
-------------------  