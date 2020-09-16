###
# Arguments
# fn: fn function.. returns a row vector.. 

#source("utils.R")
#source("F_misc.R")

P_evaluate <- function(fn, Population){
  results<-apply(Population, FUN=function(x){fn(x)}, MARGIN=1)
  results <- t(results)
  results
}

rvea6 <- function(fn, varcnt, fncnt, lowerbound, upperbound, opt, popsize=21, maxgen=110, p1=2, p2=0,  
                 alpha=2, fr=0.1, FE=0){
  
  if (opt == 0){
    optimize_func <- Min
  }
  else if (optimize_func == 1){
    optimize_func <- Max
  }
  
  Evaluations <- maxgen*popsize  # max number of fitness evaluations
  
  #reference vector initialization
  list[popsize,Vs] <- F_weight(p1,p2,fncnt)
  
  
  for (i in 1:popsize){
    Vs[i,] <- Vs[i,]/norm(Vs[i,])
  }
  V <- Vs
  
  maxgen <- floor(Evaluations/popsize)
  
  #calculat neighboring angle for angle normalization
  cosineVV <- V %*% t(V)
  list[scosineVV, neighbor] <- Sort_descend(cosineVV)
  acosVV <- acos(scosineVV[,2])
  refV <- (acosVV)
  
  # initialize population
  lowerbound <- R(lowerbound); upperbound <- R(upperbound)
  Boundary <- C(upperbound, lowerbound)
  Coding <- "Real"
  Population <- rand(popsize,varcnt)
  Population <- Population*repmat(upperbound,popsize,1)+(1-Population)*repmat(lowerbound,popsize,1)
  FunctionValue <- P_evaluate(fn, Population)
  
  all_population <- Population
  all_functionvalues <- FunctionValue
  
  for (Gene in 0 : (maxgen - 1) ){
    #random mating and reproduction
    MatingPool <- F_mating(Population)
    
    Offspring <- P_generator(MatingPool,Boundary,Coding,popsize);  
    
    FE <- FE + size(Offspring, 1)
    
    Population <- C(Population, Offspring)
    val <- P_evaluate(fn, Offspring)
    FunctionValue <- C(FunctionValue, val)
    
    #APD based selection
    theta0 <-  (Gene/(maxgen))^alpha*(fncnt)
    Selection <- F_select(FunctionValue,V, theta0, refV, optimize_func)
    Population <- Population[Selection,]
    FunctionValue <- FunctionValue[Selection,]
    
    p <- Population
    f <- FunctionValue
    if (size(p,1) != popsize){
      delta <- popsize - size(p, 1)
      if (delta > size(p,1)){
        s <- sample(size(p,1), delta, replace=T)
      }
      else{
        s <- sample(size(p, 1), delta, replace = F)
      }
      p <- C(Population, Population[s,])
      f <- C(FunctionValue, FunctionValue[s,])
    }
    all_population <- C(all_population, p)
    all_functionvalues <- C(all_functionvalues, f)
    
    #reference vector adaption
    if (Gene %% ceiling(maxgen*fr) == 0){
      if (is.null(FunctionValue)){
        print("null")
      }
      #update the reference vectors
      Zmin <- Min(FunctionValue,1)
      Zmax <- Max(FunctionValue,1)
      V <- Vs
      V <- V*repmat((Zmax - Zmin)*1.0,popsize,1)
      for (i in 1:popsize){
        V[i,] <- V[i,] / norm(V[i,])
      }
      #update the neighborning angle value for angle normalization
      cosineVV <- V%*%t(V)
      list[scosineVV, neighbor] <- Sort_descend(cosineVV)
      acosVV <- acos(scosineVV[,2])
      refV <- (acosVV)
      printf('Progress %4s%%',as.character(Gene/maxgen*100));
    }
    
    
    
  }
  
  FunctionValue <- P_evaluate(fn, Population)
  list[FrontValue,] <- P_sort(FunctionValue, 'first')
  NonDominated  <- (FrontValue == 1)
  
  Population    <- Population[NonDominated,]
  FunctionValue <- FunctionValue[NonDominated,]
  return (list(paramvalues=Population, objfnvalues=FunctionValue, 
               numsols=size(Population, 1), 
               stored_params=all_population, popFit=all_functionvalues) )
  
}