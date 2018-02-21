DTLZ1 <- function(individual,nObj,nVar){
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  gSigma <- 0
  for (subIndex in nObj:nVar) {
    gSigma <- gSigma + ((individual[subIndex] - 0.5)^2 - cos(20*pi*((individual[subIndex] - 0.5)))) 
  }
  g <- 100*(nVar- nObj + 1 + gSigma)
  
  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- 0.5* (1+g)
    
    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- obj[objectiveIndex] * individual[cosIndex]
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  (1 - individual[nObj - objectiveIndex + 1])
  }
  
  return(obj)
}

DTLZ2 <- function(individual,nObj,nVar){
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  g <- 0
  for (subIndex in nObj:nVar) {
    g <- g + (individual[subIndex] - 0.5)^2
  }
  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- (1+g)
    
    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- obj[objectiveIndex] * cos(individual[cosIndex] * pi / 2)
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  sin(individual[nObj - objectiveIndex + 1] * pi / 2)
  }
  return(obj)
}

DTLZ3 <- function(individual,nObj,nVar){
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  gSigma <- 0
  for (subIndex in nObj:nVar) {
    gSigma <- gSigma + ((individual[subIndex] - 0.5)^2 - cos(20*pi*((individual[subIndex] - 0.5)))) 
  }
  g <- 100*(nVar- nObj + 1 + gSigma)
  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- (1+g)
    
    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- obj[objectiveIndex] * cos(individual[cosIndex] * pi / 2)
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  sin(individual[nObj - objectiveIndex + 1] * pi / 2)
  }
  return(obj)
}

DTLZ4 <- function(individual,nObj,nVar,alpha){
  obj <- matrix(rep(1,nObj),nrow = nObj,ncol = 1)
  g <- 0
  for (subIndex in nObj:nVar) {
    g <- g + (individual[subIndex] - 0.5)^2
  }
  for(objectiveIndex in 1:nObj){
    obj[objectiveIndex] <- (1+g)
    
    if( (nObj-objectiveIndex) > 0){
      for(cosIndex in 1:(nObj-objectiveIndex)){
        obj[objectiveIndex] <- (obj[objectiveIndex]^alpha) * cos(individual[cosIndex] * pi / 2)
      }
    }
    if(objectiveIndex > 1)
      obj[objectiveIndex] <- obj[objectiveIndex] *  sin(individual[nObj - objectiveIndex + 1] * pi / 2)
  }
  return(obj)
}
