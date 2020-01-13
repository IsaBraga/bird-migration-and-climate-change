---------------------------------------------------------------------------
    ##separa shapes em breeding, non-breeding e permanent resident areas ##
---------------------------------------------------------------------------  


breeding_birds<- function (shape_bird)
{
  if (( 2%in%shape_bird$SEASONA) & #####2BREEDING 
      (1%in%shape_bird$ORIGIN| 2%in%shape_bird$ORIGIN) & #####1NATIVE AND 2REINTRODUCED POLYGONS
      any(shape_bird$PRESENC==1)){ 
    
    bird<-shape_bird[shape_bird$SEASONA==2 & shape_bird$PRESENC==1 & shape_bird$ORIGIN==1 | shape_bird$ORIGIN==2,]
   
     return(bird)
  }else{
    return(NULL)
  }
}



nonbreeding_birds<- function (shape_bird)
{
  if (( 3%in%shape_bird$SEASONA) & #####2BREEDING 
      (1%in%shape_bird$ORIGIN| 2%in%shape_bird$ORIGIN) & #####1NATIVE AND 2REINTRODUCED POLYGONS
      any(shape_bird$PRESENC==1)){  #####1PRESENTE
    
    bird<-shape_bird[shape_bird$SEASONA==3 & shape_bird$PRESENC==1 & shape_bird$ORIGIN==1 | shape_bird$ORIGIN==2,]
    
    return(bird)
  }else{
    return(NULL)
  }
}


resident_birds<- function (shape_bird)
{
  if (( 1%in%shape_bird$SEASONA) & #####2BREEDING 
      (1%in%shape_bird$ORIGIN| 2%in%shape_bird$ORIGIN) & #####1NATIVE AND 2REINTRODUCED POLYGONS
      any(shape_bird$PRESENC==1)){  #####1PRESENTE
    
    bird<-shape_bird[shape_bird$SEASONA==1 & shape_bird$PRESENC==1 & shape_bird$ORIGIN==1 | shape_bird$ORIGIN==2,]
    
    return(bird)
  }else{
    return(NULL)
  }
}




#Sara Varela
#Maraisa Braga
