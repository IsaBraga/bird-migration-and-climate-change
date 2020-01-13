#testes
##sensibilidade:  (se o modelo detecta bem a presença da espécie-alvo), 
a=verdadeira presenca ou verdadeiro positivo
b=falsa presenca ou falso positivo
c=falsa ausencia ou falso negativo
d=verdadeira ausencia ou verdadeiro negativo

a e d são acertos
b e c são erros
b+d=todas as ausencias #(ausencia real)
a+c=todas as presencas #(presença real)

erro de sobreprevisão=b/(b+d) #falso positivo ou erro tipo 1
erro de omissão=c/(a+c) ##falso negativo ou erro-tipo II
sensibilidade=a/(a+c) ##se o modelo detecta bem a presença da espécie-alvo
especificidade=d/(b+d)  ##se o modelo detecta bem a ausência da espécie-alvo


$minimum
[1] 0.02342558

$objective
[1] 0.4345794

rep_model #presencas e ausencias modelada/projetada  da spp (varia de 0 a 0.6778846)
reproduccion_rst2 #presencas e ausencias reais da spp
rethreshold #projecao apenas no hemisferio adequado para reproducao

#como deveria ser
#erro de omissao original
sum((validate@data@values)==1, na.rm=T)/sum((shape_rst@data@values==1), na.rm=T) 
0.1308411
#erro de sobreprevisao original
sum((validate@data@values)==2, na.rm=T)/(sum((validate@data@values)==2, na.rm=T)+sum((validate@data@values)==0, na.rm=T))
0.006232104 #enviesado pelo nº de ausencias

0.1308411+0.006232104 # novo erro

sum((validate@data@values)==0, na.rm=T) #verdadeira ausencia = d
sum((validate@data@values)==1, na.rm=T) #falsa ausencia = c
sum((validate@data@values)==2, na.rm=T) #falsa presenca = b
sum((validate@data@values)==3, na.rm=T) #verdadeira presenca = a

Onde: sobreprevisão = número de células onde ocorreram falsas presenças / número de células de 
presenças do modelo; e omissão = número de células onde ocorreram falsas ausências/ número de células de presenças do modelo.

#alto valor de threshold aumenta falsas ausencias e diminui falsas presencas
#baixo valor de threshold diminui falsas ausencias e aumenta falsas presencas
adequab_rep <- optimise(reclass_rep, interval = c(0,max(rep_model@data@values, na.rm=T)), 
                        obj_model=rep_model, 
                        shape_rst=reproduccion_rst2, cortar_rep=rethreshold)$minimum  

threshold<-0.001
threshold<-0.01
threshold<-0.66
threshold<- 0.02342558
reclass_rep<-function(threshold, obj_model, shape_rst, cortar_rep){
  obj_model_reclass<-crop(reclassify(obj_model, c(0,threshold,0,threshold,+Inf,1))*2, cortar_rep)
  validate<-shape_rst+obj_model_reclass #valores de a=3, b=2, c=1 e d=0
  
  sensibilidade <-  sum((validate@data@values)==3, na.rm=T)/sum((shape_rst@data@values==1), na.rm=T)  ##a/(a+c)
  #sobrepredicao <- sum((validate@data@values)==2, na.rm=T)/sum((shape_rst@data@values==1), na.rm=T)  ##b/(a+c)
  sobreprevisao <- sum((validate@data@values)==2, na.rm=T)/(sum((validate@data@values)==3, na.rm=T)+sum((validate@data@values)==2, na.rm=T)) ##b/(a+b)
  erro <- min((1-sensibilidade) + sobreprevisao)
  #erro <- (1-sensibilidade)*2+sobrepredicao ## 1-sensibilidade = erro de omissao

  return(erro)
}

#sugestão do Leo: b/a+b

threshold= 0.02033323
erro=0.2967604

threshold=0.001
erro=0.2771349

threshold=0.01
erro= 0.2961021


####################################

threshold<-0.01931007
threshold<-0.2
threshold<-0.4
optimise(reclass,interval = c(0,max(inv_model@data@values, na.rm=T)),obj_model=inv_model,shape_rst=invernada_rst2, cortar=america_rst0)

reclass<-function(threshold, obj_model, shape_rst, cortar){
  obj_model_reclass<-crop(reclassify(obj_model, c(0,threshold,0,threshold,+Inf,1))*2, cortar)
  validate<-shape_rst+obj_model_reclass
  
  sensibilidade <-   sum((validate@data@values)==3, na.rm=T)/sum((shape_rst@data@values==1), na.rm=T)  ##a/(a+c)
  sobreprevisao <- sum((validate@data@values)==2, na.rm=T)/(sum((validate@data@values)==3, na.rm=T)+sum((validate@data@values)==2, na.rm=T)) ##b/(a+b)
  erro <- (1-sensibilidade)+sobreprevisao
  
  return(erro)
}




america_shp<-readShapePoly("D:/artigo aves migratorias/dados/mapas/America.shp")
