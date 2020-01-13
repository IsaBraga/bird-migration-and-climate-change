

  ##FUNCAO MODELO BIOCLIM PARA TYRANIDEOS MIGRATORIOS ##
 

library (SDMTools)
library (raster)
library (dismo)
library (rgdal)
library(maptools)
library(outliers)
library(rgeos)


###script para área de residencia, invernada e reproducao separadas

##corte da adequabilidade (threshold) calculado individualmente por spp e modelo

setwd("D:/artigo aves migratorias/R/script Bioclim")
source("eliminate_outlier.R") #19 variaveis, output=data_new
source("breeding_birds.R")  #separa os shapes entre as áreas de rep, inv e res para cada spp 
source("reclass.R") #reclassify considerando o threshold

#pasta com todos os shapes e com todos os modelos climaticos
setwd("D:/artigo aves migratorias/R/dados")

#objetos necessarios
clima_present<-dir(pattern="_Modern.grd")
clima_rcp26<-dir(pattern="_rcp26.grd")
clima_rcp85<-dir(pattern="_rcp85.grd")

shapes<-dir(pattern=".shp") #lista dos 82 shapefiles
sudamerica<-raster("sudamerica.grd") #raster da América do Sul preenchido com zeros
northamerica<-raster("northamerica.grd") #raster da América do Norte preenchido com zeros
americas_rst<-raster("americas_rst.grd") #raster da América preenchido com uns
america_rst0<-raster("america_rst0.grd") #raster da América preenchido com zeros
nomes_linhas<-read.table("nomes_linhas.txt", h=F)
spp_hem<-read.csv("spp_hem.csv", h=TRUE) 
 
#**ATENCAO!! As spp Tyrannus melancholicus e Pytocephalus rubinus devem ter os dados longitudes e latidudes máximas e mínimas, mid_point e número de células analisadas separdamente pois possuem áreas de reprodução nos hemisférios norte e sul.

  
  
#### FUNCAO ####
  
run_Bioclim_model<-function(shapes, clima_present, clima_rcp26, clima_rcp85){
  
  ## Criando objeto para os resultados
  resultados_spp<- array(dim=c(66,17,82))
  dimnames(resultados_spp)[[3]]<-gsub(".shp", "", shapes)
  dimnames(resultados_spp)[[2]]<-c("long_min","lat","long_max","lat","midpoint_x", "midpoint_y","n_cells", "sensibilidade", 
                                   "sobreprevisão", "range_ganho", "range_perdido", "range_mantido", 
                                   "n_cel_subpre", "n_cel_sobrep", "n_cel_acerto1", "area km2", "adequabil_corte")
  dimnames(resultados_spp)[[1]]<-nomes_linhas[,1]
  
  
  ## Rasters para a riqueza por modelo climatico
  Sraster_rep0k<-stack(america_rst0, america_rst0, america_rst0, america_rst0, america_rst0, america_rst0)
  Sraster_rep26<- Sraster_rep0k
  Sraster_rep85<- Sraster_rep0k
  
  Sraster_inv0k<-Sraster_rep0k
  Sraster_inv26<-Sraster_rep0k
  Sraster_inv85<-Sraster_rep0k
  
  Sraster_res0k<-Sraster_rep0k
  Sraster_res26<-Sraster_rep0k
  Sraster_res85<-Sraster_rep0k
  
  ##Rasters para a riqueza de especies por spp (Ensembles  dos 6 modelos climaticos) 
  S_ensemble_rep0k<- america_rst0
  S_ensemble_rep26<- america_rst0
  S_ensemble_rep85<- america_rst0
  
  S_ensemble_inv0k<-america_rst0
  S_ensemble_inv26<-america_rst0
  S_ensemble_inv85<-america_rst0
  
  S_ensemble_res0k<-america_rst0
  S_ensemble_res26<-america_rst0
  S_ensemble_res85<-america_rst0
  
  S_ensemble_REP<- brick(S_ensemble_rep0k, S_ensemble_rep26, S_ensemble_rep85)
  S_ensemble_INV<- brick(S_ensemble_inv0k, S_ensemble_inv26, S_ensemble_inv85)
  S_ensemble_RES<- brick(S_ensemble_res0k, S_ensemble_res26, S_ensemble_res85)
  S_ensemble_raster<-c(S_ensemble_REP, S_ensemble_INV, S_ensemble_RES)
    
    
  for(i in 1:length(shapes)){
  print(shapes[i]) #informar qual spp esta sendo modelada
    ensemble_REP_0k<- america_rst0 #raster para calcular o ensemble
    ensemble_REP_26<- america_rst0
    ensemble_REP_85<- america_rst0 
    ensemble_INV_0k<- america_rst0 #raster para calcular o ensemble
    ensemble_INV_26<- america_rst0
    ensemble_INV_85<- america_rst0 
    ensemble_RES_0k<- america_rst0 #raster para calcular o ensemble
    ensemble_RES_26<- america_rst0
    ensemble_RES_85<- america_rst0 
   
  ## Separar as areas de reprodução, invernada e residencia
    paxaro<- readShapePoly (shapes[i], delete_null_obj=TRUE)
    reproduccion<-breeding_birds(paxaro)
    if(!is.null(reproduccion)){
      reproduccion<-gUnaryUnion(reproduccion)}
    invernada<- nonbreeding_birds(paxaro) 
    if(!is.null(invernada)){
      invernada<-gUnaryUnion(invernada)}
    residencia<-resident_birds(paxaro)
    if(!is.null(residencia)){
      residencia<-gUnaryUnion(residencia)} 
    
    
    ## latitude e longitude min e max e tamanho das areas
    
    if(!is.null(reproduccion)){
      reproduccion_rst<-rasterize(reproduccion, americas_rst)
      reproduccion_point<-rasterToPoints(reproduccion_rst)
      if(nrow(reproduccion_point)>0){
        resultados_spp[1,c(1:4),i]<- c(apply(reproduccion_point[,c(1:2)], 2, min),apply(reproduccion_point[,c(1:2)], 2, max))
        resultados_spp[1,c(5,6),i]<- apply(reproduccion_point[,c(1:2)],2,median)
        resultados_spp[1,7,i]<-nrow(reproduccion_point)
      }    
    }
    
    if(!is.null(invernada)){
      invernada_rst<-rasterize(invernada, americas_rst )
      invernada_point<-rasterToPoints(invernada_rst)
      if(nrow(invernada_point)>0){
        resultados_spp[2,c(1:4),i]<- c(apply(invernada_point[,c(1:2)], 2, min),apply(invernada_point[,c(1:2)], 2, max))
        resultados_spp[2,c(5,6),i]<- apply(invernada_point[,c(1:2)],2,median)
        resultados_spp[2,7,i]<-nrow(invernada_point)
      }
    }
    
    
    if(!is.null(residencia)){
      residencia_rst<-rasterize(residencia, americas_rst )
      residencia_point<-rasterToPoints(residencia_rst)
      if(nrow(residencia_point)>0){
        resultados_spp[3,c(1:4),i]<- c(apply(residencia_point[,c(1:2)], 2, min),apply(residencia_point[,c(1:2)], 2, max))
        resultados_spp[3,c(5,6),i]<- apply(residencia_point[,c(1:2)],2,median)
        resultados_spp[3,7,i]<-nrow(residencia_point) #area do shape em numero de celulas
      }
    }  
    
  
    
    #recortes necessarios para as projecoes das areas de reproducao no continente adequado
    if(spp_hem[i,3]=="sul" & !is.na(spp_hem[i,3])){ 
      recorte<-sudamerica
    }
    
    if(spp_hem[i,3]=="sul_norte" & !is.na(spp_hem[i,3])){ 
      recorte<-america_rst0
    }
    
    if(spp_hem[i,3]=="norte" & !is.na(spp_hem[i,3])){ 
      recorte<-northamerica
    }
    
    
    ### Modelos climaticos ###
    
    for(j in 1:length(clima_present)){
      clima_Modern<-brick(clima_present[j])
      clima_LOW<-brick(clima_rcp26[j])
      clima_HIGH<-brick(clima_rcp85[j])
         
      
      #### MODELO AREA DE REPRODUCAO ####
      if(!is.null(reproduccion)){
        if(nrow(reproduccion_point)>1){
          
          #### PRESENTE ####
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_present[j]),"REP", sep=" "))
          
          # nicho climatico
          reproduccion_rst2<-merge(reproduccion_rst, recorte, ext=recorte@extent) #raster da area original sobre o raster da america do sul ou do norte
          data_new_rep<-eliminate_outliers(climate=clima_Modern, distribution=reproduccion)
          map_rep<- bioclim (data_new_rep)
          rep_model<-predict(map_rep, clima_Modern) #presenca da spp modelada/projetada
          adequab_rep <- optimise(reclass_rep, interval = c(0,max(rep_model@data@values, na.rm=T)), obj_model=rep_model, 
                                  shape_rst=reproduccion_rst2, cortar_rep=recorte)$minimum  #retorna $minimum que significa o valor do corte de adequabilidade (threshold) que minimiza o erro
          
          final_rep0k<-crop(reclassify(rep_model, c(0,adequab_rep,0,adequab_rep,+Inf,1))*2, recorte)+recorte
          
          
          ## validacao area de reproducao
          validate_repro<-reproduccion_rst2+final_rep0k #valores de a=3, b=2, c=1 e d=0
          #performance do modelo
          resultados_spp[(j-1)*9+4,8,i]<- sum(validate_repro@data@values==3, na.rm=T)/ sum(sum(validate_repro@data@values==3, na.rm=T, sum(validate_repro@data@values==1, na.rm=T))) #% acerto = sensibilidade = a/(a+c)
          resultados_spp[(j-1)*9+4,9,i]<- sum((validate_repro@data@values)==2, na.rm=T)/sum(sum((validate_repro@data@values)==3, na.rm=T, sum((validate_repro@data@values)==2, na.rm=T))) #% sobreprevisao = b/a+b
          resultados_spp[(j-1)*9+4,13,i]<-sum((validate_repro@data@values)==1, na.rm=T) #n_cel_subpre = c 
          resultados_spp[(j-1)*9+4,14,i]<-sum((validate_repro@data@values)==2, na.rm=T) #n_cel_sobrep = b
          resultados_spp[(j-1)*9+4,15,i]<-sum((validate_repro@data@values)==3, na.rm=T) #n_cel_acerto = a
          #resultados_spp[(j-1)*9+4,16,i]<-21014-(sum((validate_repro@data@values==1))+sum((validate_repro@data@values==2))+sum((validate_repro@data@values==3))) #n_cel_acerto0  
          
          resultados_spp[(j-1)*9+4,17,i]<- adequab_rep # valor do threshold
          
          # localizaçao, centroide e area
          rep_model_point<-rasterToPoints(final_rep0k)
          rep_model_point<-matrix(rep_model_point[rep_model_point[,3]==2,1:2],ncol=2)
          
          if(nrow(rep_model_point)>0){
            resultados_spp[(j-1)*9+4,c(1:4),i]<- c(apply(rep_model_point, 2, min),apply(rep_model_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+4,c(5,6),i]<- apply(rep_model_point,2,median) #centroide
            resultados_spp[(j-1)*9+4,7,i]<-sum((final_rep0k@data@values)==2, na.rm=T) #area do modelo
          }
          
         
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_rep0k@layers[[j]]<-merge(final_rep0k/2, america_rst0)+ Sraster_rep0k@layers[[j]]
          
          ## ENSEMBLE
          ensemble_REP_0k<-ensemble_REP_0k + merge(final_rep0k/2, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
         
          
          
         #### CENARIO DE BAIXA EMISSAO ####    
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_rcp26[j]),"REP", sep=" "))  
          
          rep_model26<-predict(map_rep,clima_LOW)
          final_rep26<-crop(reclassify(rep_model26, c(0,adequab_rep,0,adequab_rep,+Inf,1)), recorte)+recorte
          
          #localizaçao, centroide e area
          rep_model26_point<-rasterToPoints(final_rep26)
          rep_model26_point<-matrix(rep_model26_point[rep_model26_point[,3]==1,c(1:2)],ncol=2)
          
          if(nrow(rep_model26_point)>0){
            resultados_spp[(j-1)*9+5,c(1:4),i]<- c(apply(rep_model26_point, 2, min),apply(rep_model26_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+5,c(5,6),i]<- apply(rep_model26_point,2,median) #centroide 
          resultados_spp[(j-1)*9+5,7,i]<-sum((final_rep26@data@values)==1, na.rm=T) # area
          }
          
          #sobreposicao do range 26 x 0k
          sobreposicao_rep26<-final_rep0k+final_rep26
          resultados_spp[(j-1)*9+5,10,i]<-sum((sobreposicao_rep26@data@values)==1, na.rm=T)/sum((final_rep0k@data@values)==2, na.rm=T) # % range ganho
          resultados_spp[(j-1)*9+5,11,i]<-sum((sobreposicao_rep26@data@values)==2, na.rm=T)/sum((final_rep0k@data@values)==2, na.rm=T) # % range perdido
          resultados_spp[(j-1)*9+5,12,i]<-sum((sobreposicao_rep26@data@values)==3, na.rm=T)/sum((final_rep0k@data@values)==2, na.rm=T) # % range mantido
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_rep26@layers[[j]]<-merge(final_rep26, america_rst0)+Sraster_rep26@layers[[j]]
          ## ENSEMBLE
          ensemble_REP_26<-ensemble_REP_26 + merge(final_rep26, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
         
          
          
         #### CENARIO DE ALTA EMISSAO ####
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_rcp85[j]),"REP", sep=" "))
          
          rep_model85<-predict(map_rep,clima_HIGH)
          final_rep85<-crop(reclassify(rep_model85, c(0,adequab_rep,0,adequab_rep,+Inf,1)), recorte)+recorte
          
          #sobreposicao do range 85 x 0k
          sobreposicao_rep85<-final_rep0k+final_rep85
          resultados_spp[(j-1)*9+6,10,i]<-sum((sobreposicao_rep85@data@values)==1, na.rm=T)/sum((final_rep0k@data@values)==2, na.rm=T) # % range ganho
          resultados_spp[(j-1)*9+6,11,i]<-sum((sobreposicao_rep85@data@values)==2, na.rm=T)/sum((final_rep0k@data@values)==2, na.rm=T) # % range perdido
          resultados_spp[(j-1)*9+6,12,i]<-sum((sobreposicao_rep85@data@values)==3, na.rm=T)/sum((final_rep0k@data@values)==2, na.rm=T) # % range mantido
          
          # localizaçao, centroide e area
          rep_model85_point<-rasterToPoints(final_rep85)
          rep_model85_point<-matrix(rep_model85_point[rep_model85_point[,3]==1,c(1:2)], ncol=2)
          
          if(nrow(rep_model85_point)>0){
            resultados_spp[(j-1)*9+6,c(1:4),i]<- c(apply(rep_model85_point, 2, min),apply(rep_model85_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+6,c(5,6),i]<- apply(rep_model85_point,2,median) #centroide
          resultados_spp[(j-1)*9+6,7,i]<-sum((final_rep85@data@values)==1, na.rm=T) #area
          }
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_rep85@layers[[j]]<-merge(final_rep85, america_rst0)+Sraster_rep85@layers[[j]]
          
          ## ENSEMBLE
          ensemble_REP_85<-ensemble_REP_85 + merge(final_rep85, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
          
        
          #rm(reproduccion_rst2, data_new_rep, map_rep, rep_model, adequab_rep, final_rep0k, validate_repro, rep_model_point, clima_LOW, rep_model26,final_rep26, rep_model26_point, clima_HIGH,  rep_model85, final_rep85, rep_model85_point)
        } 
     }#fecha o if do modelo area de reproducao
      
      
      
      ##### MODELO AREA DE INVERNADA #####
      if(!is.null(invernada)){
        if( nrow(invernada_point)>1){
          
         #### PRESENTE ####
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_present[j]),"INV", sep=" "))
          
          # nicho climatico
          invernada_rst2<-merge(invernada_rst,america_rst0)
          data_new_inv<-eliminate_outliers(climate=clima_Modern, distribution=invernada)
          map_inv<- bioclim (data_new_inv)
          inv_model<-predict(map_inv,clima_Modern)
          adequab_inv <- optimise(reclass,interval = c(0,max(inv_model@data@values, na.rm=T)),obj_model=inv_model,shape_rst=invernada_rst2, cortar=america_rst0)$minimum
          final_inv0k<-crop(reclassify(inv_model, c(0,adequab_inv,0,adequab_inv,+Inf,1))*2, america_rst0)+america_rst0 
          
          # validacao area de invernada
          validate_inv<-invernada_rst2+final_inv0k
          # performance do modelo: acertos e sobrepredicao
          resultados_spp[(j-1)*9+7,8,i]<- sum((validate_inv@data@values)==3, na.rm=T)/sum((invernada_rst2@data@values)==1, na.rm=T) #% acerto = sensibilidade
          resultados_spp[(j-1)*9+7,9,i]<- sum((validate_inv@data@values)==2, na.rm=T)/ sum(sum((validate_inv@data@values)==3, na.rm=T),sum((validate_inv@data@values)==2, na.rm=T)) #% sobrepredicao 
          resultados_spp[(j-1)*9+7,13,i]<-sum((validate_inv@data@values)==1, na.rm=T) #n_cel_subpre
          resultados_spp[(j-1)*9+7,14,i]<-sum((validate_inv@data@values)==2, na.rm=T) #n_cel_sobrep
          resultados_spp[(j-1)*9+7,15,i]<-sum((validate_inv@data@values)==3, na.rm=T) #n_cel_acerto1
          #resultados_spp[(j-1)*9+7,16,i]<-21014-(sum((validate_inv@data@values==1))+sum((validate_inv@data@values==2))+sum((validate_inv@data@values==3))) #n_cel_acerto0 
          resultados_spp[(j-1)*9+7,17,i]<-adequab_inv
          
          # localizaçao, centroide e area
          inv_model_point<-rasterToPoints(final_inv0k)
          inv_model_point<-matrix(inv_model_point[inv_model_point[,3]==2,c(1:2)],ncol=2)
          
          if(nrow(inv_model_point)>0){
            resultados_spp[(j-1)*9+7,c(1:4),i]<- c(apply(inv_model_point, 2, min),apply(inv_model_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+7,c(5,6),i]<- apply(inv_model_point,2,median) #centroide
            resultados_spp[(j-1)*9+7,7,i]<-sum((final_inv0k@data@values)==2, na.rm=T) #area do modelo
          }
          
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_inv0k@layers[[j]]<-merge(final_inv0k/2, america_rst0)+Sraster_inv0k@layers[[j]]
          
          ## ENSEMBLE
          ensemble_INV_0k<-ensemble_INV_0k + merge(final_inv0k/2, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
          
          
         #### CENARIO DE BAIXA EMISSAO ####  
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_rcp26[j]),"INV", sep=" "))
          
         
          inv_model26<-predict(map_inv,clima_LOW)
          final_inv26<-crop(reclassify(inv_model26, c(0,adequab_inv,0,adequab_inv,+Inf,1)), america_rst0)+america_rst0
          
          inv_model26_point<-rasterToPoints(final_inv26)
          inv_model26_point<-matrix(inv_model26_point[inv_model26_point[,3]==1,c(1:2)],ncol=2)
          
          # localizaçao, centroide e area
          if(nrow(inv_model26_point)>0){
            resultados_spp[(j-1)*9+8,c(1:4),i]<- c(apply(inv_model26_point, 2, min),apply(inv_model26_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+8,c(5,6),i]<- apply(inv_model26_point,2,median) #centroide
          resultados_spp[(j-1)*9+8,7,i]<-sum((final_inv26@data@values)==1, na.rm=T) #area do modelo
          }
          
          #sobreposicao do range 26 x 0k
          sobreposicao_inv26<-final_inv0k+final_inv26
          resultados_spp[(j-1)*9+8,10,i]<-sum((sobreposicao_inv26@data@values)==1, na.rm=T)/sum((final_inv0k@data@values)==2, na.rm=T) # % range ganho
          resultados_spp[(j-1)*9+8,11,i]<-sum((sobreposicao_inv26@data@values)==2, na.rm=T)/sum((final_inv0k@data@values)==2, na.rm=T) # % range perdido
          resultados_spp[(j-1)*9+8,12,i]<-sum((sobreposicao_inv26@data@values)==3, na.rm=T)/sum((final_inv0k@data@values)==2, na.rm=T) # % range mantido
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_inv26@layers[[j]]<-merge(final_inv26, america_rst0)+Sraster_inv26@layers[[j]]
          
          ## ENSEMBLE
          ensemble_INV_26<-ensemble_INV_26 + merge(final_inv26, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
          
          
          #### CENARIO DE ALTA EMISSAO ####
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_rcp85[j]),"INV", sep=" "))
          
         
          inv_model85<-predict(map_inv,clima_HIGH)
          final_inv85<-crop(reclassify(inv_model85, c(0,adequab_inv,0,adequab_inv,+Inf,1)), america_rst0)+america_rst0
          
          #sobreposicao do range 85 x 0k
          sobreposicao_inv85<-final_inv0k+final_inv85
          resultados_spp[(j-1)*9+9,10,i]<-sum((sobreposicao_inv85@data@values)==1, na.rm=T)/sum((final_inv0k@data@values)==2, na.rm=T)
          resultados_spp[(j-1)*9+9,11,i]<-sum((sobreposicao_inv85@data@values)==2, na.rm=T)/sum((final_inv0k@data@values)==2, na.rm=T)
          resultados_spp[(j-1)*9+9,12,i]<-sum((sobreposicao_inv85@data@values)==3, na.rm=T)/sum((final_inv0k@data@values)==2, na.rm=T)
          
          # localizaçao, centroide e area
          inv_model85_point<-rasterToPoints(final_inv85)
          inv_model85_point<- matrix(inv_model85_point[inv_model85_point[,3]==1,c(1:2)], ncol=2)
          
          if(nrow(inv_model85_point)>0){
            resultados_spp[(j-1)*9+9,c(1:4),i]<- c(apply(inv_model85_point, 2, min),apply(inv_model85_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+9,c(5:6),i]<- apply(inv_model85_point,2,median) #centroide
          resultados_spp[(j-1)*9+9,7,i]<-sum((final_inv85@data@values)==1, na.rm=T) #area
          }
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_inv85@layers[[j]]<-merge(final_inv85, america_rst0)+Sraster_inv85@layers[[j]]
          
          ## ENSEMBLE
          ensemble_INV_85<-ensemble_INV_85 + merge(final_inv85, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
          
          
         # rm(clima_present, invernada_rst2, data_new_inv, map_inv, inv_model, adequab_inv, final_inv0k, validate_inv, inv_model_point,clima_LOW, inv_model26,final_inv26, inv_model26_point, clima_HIGH,  inv_model85, final_inv85, inv_model85_point)
        } 
      }#fecha o if do modelo area de invernada
      
      
      
      ##### MODELO AREA DE RESIDENCIA #####
      if(!is.null(residencia)){
        if( nrow(residencia_point)>1){
          
          #### PRESENTE ####
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_present[j]),"RES", sep=" "))
          
          # nicho climatico
          residencia_rst2<-merge(residencia_rst,america_rst0)
          data_new_res<-eliminate_outliers(climate=clima_Modern, distribution=residencia)
          map_res<- bioclim (data_new_res)
          res_model<-predict(map_res,clima_Modern)
          adequab_res <- optimise(reclass,interval = c(0,max(res_model@data@values, na.rm=T)),obj_model=res_model,shape_rst=residencia_rst2, cortar=america_rst0)$minimum
          final_res0k<-crop(reclassify(res_model, c(0,adequab_res,0,adequab_res,+Inf,1))*2, america_rst0)+america_rst0 
          
          # validacao area de residencia
          validate_res<-residencia_rst2+final_res0k
          # performance do modelo: sensibilidade e sobrepredicao
          resultados_spp[(j-1)*9+10,8,i]<- sum((validate_res@data@values)==3, na.rm=T)/sum((residencia_rst2@data@values)==1, na.rm=T) #% acerto = sensibilidade
          resultados_spp[(j-1)*9+10,9,i]<- sum((validate_res@data@values)==2, na.rm=T)/ sum(sum((validate_res@data@values)==3, na.rm=T),sum((validate_res@data@values)==2, na.rm=T)) #% sobrepredicao 
          resultados_spp[(j-1)*9+10,13,i]<-sum((validate_res@data@values)==1, na.rm=T) #n_cel_subpre
          resultados_spp[(j-1)*9+10,14,i]<-sum((validate_res@data@values)==2, na.rm=T) #n_cel_sobrep
          resultados_spp[(j-1)*9+10,15,i]<-sum((validate_res@data@values)==3, na.rm=T) #n_cel_acerto1
          #resultados_spp[(j-1)*9+10,16,i]<-21014-(sum((validate_res@data@values==1))+sum((validate_res@data@values==2))+sum((validate_res@data@values==3))) #n_cel_acerto0 
          resultados_spp[(j-1)*9+10,17,i]<-adequab_res
          
          # localizaçao, centroide e area
          res_model_point<-rasterToPoints(final_res0k)
          res_model_point<-matrix(res_model_point[res_model_point[,3]==2,c(1:2)],ncol=2)
          
          if(nrow(res_model_point)>0){
            resultados_spp[(j-1)*9+10,c(1:4),i]<- c(apply(res_model_point, 2, min),apply(res_model_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+10,c(5,6),i]<- apply(res_model_point,2,median) #centroide
            resultados_spp[(j-1)*9+10,7,i]<-sum((final_res0k@data@values)==2, na.rm=T) #area do modelo
          }
          
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_res0k@layers[[j]]<-merge(final_res0k/2, america_rst0)+Sraster_res0k@layers[[j]]
          
          ## ENSEMBLE
          ensemble_RES_0k<-ensemble_RES_0k + merge(final_res0k/2, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
          
          
         #### CENARIO DE BAIXA EMISSAO ####  
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_rcp26[j]),"RES", sep=" "))
          
         
          res_model26<-predict(map_res,clima_LOW)
          final_res26<-crop(reclassify(res_model26, c(0,adequab_res,0,adequab_res,+Inf,1)), america_rst0)+america_rst0
          
          res_model26_point<-rasterToPoints(final_res26)
          res_model26_point<-matrix(res_model26_point[res_model26_point[,3]==1,c(1:2)],ncol=2)
          
          # localizaçao, centroide e area
          if(nrow(res_model26_point)>0){
            resultados_spp[(j-1)*9+11,c(1:4),i]<- c(apply(res_model26_point, 2, min),apply(res_model26_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+11,c(5,6),i]<- apply(res_model26_point,2,median) #centroide
          resultados_spp[(j-1)*9+11,7,i]<-sum((final_res26@data@values)==1, na.rm=T) #area
          }
          
          #sobreposicao do range 26 x 0k
          sobreposicao_res26<-final_res0k+final_res26
          resultados_spp[(j-1)*9+11,10,i]<-sum((sobreposicao_res26@data@values)==1, na.rm=T)/sum((final_res0k@data@values)==2, na.rm=T) # % range ganho
          resultados_spp[(j-1)*9+11,11,i]<-sum((sobreposicao_res26@data@values)==2, na.rm=T)/sum((final_res0k@data@values)==2, na.rm=T) # % range perdido
          resultados_spp[(j-1)*9+11,12,i]<-sum((sobreposicao_res26@data@values)==3, na.rm=T)/sum((final_res0k@data@values)==2, na.rm=T) # % range mantido
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_res26@layers[[j]]<-merge(final_res26, america_rst0)+Sraster_res26@layers[[j]]
          
          ## ENSEMBLE
          ensemble_RES_26<-ensemble_RES_26 + merge(final_res26, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
          
          
          #### CENARIO DE ALTA EMISSAO ####
          print(paste(gsub(".shp", "", shapes[i]), gsub(".grd", "", clima_rcp85[j]),"RES", sep=" "))
          
         
          res_model85<-predict(map_res,clima_HIGH)
          final_res85<-crop(reclassify(res_model85, c(0,adequab_res,0,adequab_res,+Inf,1)), america_rst0)+america_rst0
          
          #sobreposicao do range 85 x 0k
          sobreposicao_res85<-final_res0k+final_res85
          resultados_spp[(j-1)*9+12,10,i]<-sum((sobreposicao_res85@data@values)==1, na.rm=T)/sum((final_res0k@data@values)==2, na.rm=T)
          resultados_spp[(j-1)*9+12,11,i]<-sum((sobreposicao_res85@data@values)==2, na.rm=T)/sum((final_res0k@data@values)==2, na.rm=T)
          resultados_spp[(j-1)*9+12,12,i]<-sum((sobreposicao_res85@data@values)==3, na.rm=T)/sum((final_res0k@data@values)==2, na.rm=T)
          
          # localizaçao, centroide e area
          res_model85_point<-rasterToPoints(final_res85)
          res_model85_point<- matrix(res_model85_point[res_model85_point[,3]==1,c(1:2)], ncol=2)
          
          if(nrow(res_model85_point)>0){
            resultados_spp[(j-1)*9+12,c(1:4),i]<- c(apply(res_model85_point, 2, min),apply(res_model85_point, 2, max)) #lat e long max e min
            resultados_spp[(j-1)*9+12,c(5:6),i]<- apply(res_model85_point,2,median) #centroide
          resultados_spp[(j-1)*9+12,7,i]<-sum((final_res85@data@values)==1, na.rm=T) #area
          }
          
          ## RIQUEZA POR MODELO CLIMATICO (6 rasters separados)
          Sraster_res85@layers[[j]]<-merge(final_res85, america_rst0)+Sraster_res85@layers[[j]]
          
          ## ENSEMBLE
          ensemble_RES_85<-ensemble_RES_85 + merge(final_res85, america_rst0) #para calcular o ensemble dos seis modelos climaticos por spp
          
         
          
         # rm(clima_residencia, residencia_rst2, data_new_res, map_res, res_model, adequab_res, final_res0k, validate_resid, res_model_point,clima_LOW, res_model26,final_res26, res_model26_point, clima_HIGH,  res_model85, final_res85, res_model85_point)
          
        } #fecha o if do modelo area de residencia
      }
    }#fecha o for dos modelos climaticos
    rm(recorte)
    
    
    
     #### ENSEMBLE ####
     
  #agrupar os Ensembles da spp para validação e obteção dos resultados no for a seguir
  ensemble_REP<- brick(ensemble_REP_0k, ensemble_REP_26, ensemble_REP_85)
  #ensemble_REP@data@names<- c("ensemble_0k", "ensemble_26", "ensemble_85")
  ensemble_INV<- brick(ensemble_INV_0k, ensemble_INV_26, ensemble_INV_85)
  #ensemble_INV@data@names<- c("ensemble_0k", "ensemble_26", "ensemble_85")
  ensemble_RES<- brick(ensemble_RES_0k,ensemble_RES_26, ensemble_RES_85 )
  #ensemble_RES@data@names<- c("ensemble_0k", "ensemble_26", "ensemble_85")
  ensemble_raster<-c(ensemble_REP, ensemble_INV, ensemble_RES)
  
  
  #agrupar os rasters originais de cada spp para validar o Ensemble no for
  if(is.null(reproduccion)){
    reproduccion_rst<-NULL ;
   if (!is.null(reproduccion))
      reproduccion_rst<-merge(reproduccion_rst,america_rst0)
  }
  if(is.null(invernada)){
    invernada_rst2<-NULL 
  }
  if(is.null(residencia)){
    residencia_rst2<-NULL
  }
  
  raster_original<-list(reproduccion_rst , invernada_rst2, residencia_rst2)
  ensemble_raster_FINAL<-ensemble_raster #modelo para receber os ensembles já reclassificados para salvar o stack de cada spp
    
    
  #for para calcular e obter dados dos Ensembles (Ensemble com no mínimo 3 modelos climáticos sobrepostos)
    for (k in 1: length(ensemble_raster)) {
      if(!is.null(raster_original[[k]])){
    
        #### PRESENTE ####
    ensemble_raster0k<-reclassify(ensemble_raster[[k]]$layer.1, c(-Inf,2,0,2.9,6,1,6,+Inf, NA))
    S_ensemble_raster[[k]]$layer.1<-S_ensemble_raster[[k]]$layer.1 + ensemble_raster0k  #riqueza de spp
    
    # localizaçao, centroide e area 0k ensemble
    ensemble_point0k<-rasterToPoints(ensemble_raster0k)
    ensemble_point0k<-matrix(ensemble_point0k[ensemble_point0k[,3]==1,c(1:2)],ncol=2)
    ensemble_cellsize0k <- area(cut(ensemble_raster0k, breaks=c(0,9:1)), na.rm=TRUE,weights=FALSE)
    if(!is.null(ensemble_raster0k)){
      resultados_spp[(k-1)*3+58,16,i]<-sum(ensemble_cellsize0k@data@values, na.rm=T) #area em Km2
    }                  
    
    if(nrow(ensemble_point0k)>0){
      resultados_spp[(k-1)*3+58,c(1:4),i]<- c(apply(ensemble_point0k, 2, min),apply(ensemble_point0k, 2, max,na.rm=T)) #lat e long max e min
      resultados_spp[(k-1)*3+58,c(5,6),i]<- apply(ensemble_point0k,2,median) #centroide
      resultados_spp[(k-1)*3+58,7,i]<-sum((ensemble_raster0k@data@values)==1,na.rm=T) #area do modelo em celulas
     
    }
    
    
    ## validacao 
    validate_ensemble<-raster_original[[k]]+(ensemble_raster0k*2) #a=3, b=2, c=1
    
    #performance do Ensemble
    resultados_spp[(k-1)*3+58,8,i]<- sum((validate_ensemble@data@values)==3, na.rm=T)/sum(sum((validate_ensemble@data@values)==3, na.rm=T), sum((validate_ensemble@data@values)==1, na.rm=T))  #% acerto = sensibilidade = a/(a+c)
    resultados_spp[(k-1)*3+58,9,i]<- sum((validate_ensemble@data@values)==2,na.rm=T)/sum(sum((validate_ensemble@data@values)==3,na.rm=T), sum((validate_ensemble@data@values)==2,na.rm=T)) #% sobrepredicao = b/ (a+b) 
  
 
    resultados_spp[(k-1)*3+58,13,i]<-sum((validate_ensemble@data@values)==1,na.rm=T) #n_cel_subpre
    resultados_spp[(k-1)*3+58,14,i]<-sum((validate_ensemble@data@values)==2,na.rm=T) #n_cel_sobrep
    resultados_spp[(k-1)*3+58,15,i]<-sum((validate_ensemble@data@values)==3,na.rm=T) #n_cel_acerto1
    #resultados_spp[20,16,i]<-21014-(sum((validate_ensemble@data@values)==1,na.rm=T)+sum((validate_ensemble@data@values)==2,na.rm=T)+sum((validate_ensemble@data@values)==3,na.rm=T)) #n_cel_acerto0  
    
    #para salvar o ensemble de casa spp
    ensemble_raster_FINAL[[k]]$layer.1<- ensemble_raster0k
    
    
     #### CENARIO DE BAIXA EMISSAO ####
    ensemble_raster26<-reclassify(ensemble_raster[[k]]$layer.2, c(-Inf,2,0,2.9,6,1,6,+Inf, NA),na.rm=T)
    S_ensemble_raster[[k]]$layer.2<-S_ensemble_raster[[k]]$layer.2 + ensemble_raster26  #riqueza de spp
    
    # localizaçao, centroide RCP 26 ensemble
    ensemble_point26<-rasterToPoints(ensemble_raster26)
    ensemble_point26<-matrix(ensemble_point26[ensemble_point26[,3]==1,1:2],ncol=2)
    ensemble_cellsize26 <- area(cut(ensemble_raster26, breaks=c(0,9:1)), na.rm=TRUE,weights=FALSE)
    
    if(!is.null(ensemble_raster26)){
    resultados_spp[(k-1)*3+59,16,i]<-sum(ensemble_cellsize26@data@values, na.rm=T) #area em Km2
    }
    
     if(nrow(ensemble_point26)>0){
      resultados_spp[(k-1)*3+59,c(1:4),i]<- c(apply(ensemble_point26, 2, min),apply(ensemble_point26, 2, max)) #lat e long max e min
      resultados_spp[(k-1)*3+59,c(5,6),i]<- apply(ensemble_point26,2,median) #centroide
    resultados_spp[(k-1)*3+59,7,i]<-sum((ensemble_raster26@data@values)==1,na.rm=T) #area do modelo em celulas
    }
    
    
    #sobreposicao do range ensemble 26 x 0k
    sobreposicao_Ens26<-(ensemble_raster0k*2)+ensemble_raster26
    resultados_spp[(k-1)*3+59,10,i]<-sum((sobreposicao_Ens26@data@values)==1,na.rm=T)/sum((ensemble_raster0k@data@values)==1,na.rm=T) #range ganho
    resultados_spp[(k-1)*3+59,11,i]<-sum((sobreposicao_Ens26@data@values)==2,na.rm=T)/sum((ensemble_raster0k@data@values)==1,na.rm=T) #range perdido
    resultados_spp[(k-1)*3+59,12,i]<-sum((sobreposicao_Ens26@data@values)==3,na.rm=T)/sum((ensemble_raster0k@data@values)==1,na.rm=T) #% range mantido
    
    #para salvar o ensemble de casa spp
    ensemble_raster_FINAL[[k]]$layer.2<- ensemble_raster26
  
    
    #### CENARIO DE ALTA EMISSAO ####
    ensemble_raster85<-reclassify(ensemble_raster[[k]]$layer.3, c(-Inf,2,0,2.9,6,1,6,+Inf, NA),na.rm=T)
    S_ensemble_raster[[k]]$layer.3 <- S_ensemble_raster[[k]]$layer.3 + ensemble_raster85  #riqueza de spp
    
    # localizaçao, centroide e area 85 ensemble
    ensemble_point85<-rasterToPoints(ensemble_raster85)
    ensemble_point85<-matrix(ensemble_point85[ensemble_point85[,3]==1,1:2],ncol=2)
    ensemble_cellsize85 <- area(cut(ensemble_raster85, breaks=c(0,9:1)), na.rm=TRUE,weights=F)
    if(!is.null(ensemble_raster85)){
      resultados_spp[(k-1)*3+60,16,i]<-sum(ensemble_cellsize85@data@values, na.rm=T) #area em Km2
    } 
    
    if(nrow(ensemble_point85)>0){
      resultados_spp[(k-1)*3+60,c(1:4),i]<- c(apply(ensemble_point85, 2, min),apply(ensemble_point85, 2, max)) #lat e long max e min
      resultados_spp[(k-1)*3+60,c(5,6),i]<- apply(ensemble_point85,2,median) #centroide
    resultados_spp[(k-1)*3+60,7,i]<-sum((ensemble_raster85@data@values)==1,na.rm=T) #area do modelo em celulas
    
    }
    
    #sobreposicao do range ensemble 85 x 0k
    sobreposicao_Ens85<-(ensemble_raster0k*2)+ensemble_raster85
    resultados_spp[(k-1)*3+60,10,i]<-sum((sobreposicao_Ens85@data@values)==1,na.rm=T)/sum((ensemble_raster0k@data@values)==1,na.rm=T) #range ganho
    resultados_spp[(k-1)*3+60,11,i]<-sum((sobreposicao_Ens85@data@values)==2,na.rm=T)/sum((ensemble_raster0k@data@values)==1,na.rm=T) #range perdido
    resultados_spp[(k-1)*3+60,12,i]<-sum((sobreposicao_Ens85@data@values)==3,na.rm=T)/sum((ensemble_raster0k@data@values)==1,na.rm=T) #range mantido
    
    #para salvar o ensemble de casa spp
    ensemble_raster_FINAL[[k]]$layer.3<- ensemble_raster85
    
    
      } #fecha o if das areas (se a area existe)
    } #fecha o for do Ensemble
    
  
    ##para salvar o ensemble de casa spp
    ensemble_raster_FINAL<-stack(ensemble_raster_FINAL)
    names(ensemble_raster_FINAL)<- c("REP_0k", "REP_26", "REP_85","INV_0k", "INV_26", "INV_85", "RES_0k", "RES_26", "RES_85" )
    writeRaster(ensemble_raster_FINAL, paste("./ensemble_rasters/", gsub(".shp", "", shapes[i])), overwrite=TRUE)
    
    
  } #fecha o for dos 82 shapefiles
  
  #nomes para os layers do raster brick com todos os modelos climaticos separados 
  names(Sraster_rep0k) <-c("CCSM", "FGOALS", "GISS", "IPSL", "MIROC", "MRI")
  names(Sraster_rep26) <- names(Sraster_rep0k)
  names(Sraster_rep85) <- names(Sraster_rep0k)
  names(Sraster_inv0k) <- names(Sraster_rep0k)
  names(Sraster_inv26) <- names(Sraster_rep0k)
  names(Sraster_inv85) <- names(Sraster_rep0k)
  names(Sraster_res0k) <- names(Sraster_rep0k)
  names(Sraster_res26) <- names(Sraster_rep0k)
  names(Sraster_res85) <- names(Sraster_rep0k)
  
  S_raster_rep <- stack(Sraster_rep0k, Sraster_rep26, Sraster_rep85)
  S_raster_inv <- stack(Sraster_inv0k, Sraster_inv26, Sraster_inv85)
  S_raster_res <- stack(Sraster_res0k, Sraster_res26, Sraster_res85)
  S_raster <- c(S_raster_rep, S_raster_inv, S_raster_res)
  names(S_raster) <- c("S_raster_rep","S_raster_inv","S_raster_res")
  
  #nomes para os layers do raster stack Ensembles da riqueza
  names(S_ensemble_raster)<- c("S_ensemble_REP", "S_ensemble_INV", "S_ensemble_RES")
  names(S_ensemble_raster[[1]])<- c("S_ensemble_rep0k", "S_ensemble_rep26", "S_ensemble_rep85")
  names(S_ensemble_raster[[2]])<- c("S_ensemble_inv0k", "S_ensemble_inv26", "S_ensemble_inv85")
  names(S_ensemble_raster[[3]])<- c("S_ensemble_res0k", "S_ensemble_res26", "S_ensemble_res85")
   
   
  return(list(resultados_spp=resultados_spp, S_ensemble_raster=S_ensemble_raster, S_raster=S_raster))
  
} #fecha a funcao run_Bioclim_model
 



################
#Rodar modelo: 19 variáveis de interesse (bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19).
#função  eliminate_outlier para 19 variáveis

aves_migratorias_model_2<-run_Bioclim_model(shapes=shapes, clima_present=clima_present, clima_rcp26=clima_rcp26, clima_rcp85=clima_rcp85) 


save(aves_migratorias_model_2, file="model2.Rdata")

#################  

## Maraísa Braga 
## Sara Varela  



