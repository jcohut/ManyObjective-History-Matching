# The selection function in RVEA
#function [Selection] = F_select(FunctionValue, V, theta0, refV)
# Completed!
F_select <- function(FunctionValue, V, theta0, refV, optimize_func){
  
  # disable this entire segment
  #   values <- readMat("function_value.mat")
  #   FunctionValue <- values$FunctionValue
  #   V <- values$V
  #   theta0 <- values$theta0
  #   refV <- values$refV
  
  NM <- size(FunctionValue)
  N <- NM[1]
  M <- NM[2]
  VN <- size(V, 1)
  
  # only name is Zmin, but it depends upon the optimize_func.. can be min or max
  Zmin <- optimize_func(FunctionValue ,1)
  
  #Translation
  FunctionValue <- (FunctionValue - repmat(Zmin, R(size(FunctionValue,1), 1)) )
  
  #Solutions associattion to reference vectors
  
  div <- repmat(sqrt(Sum(FunctionValue^2,2)), R(1, M))
  uFunctionValue <- FunctionValue / div
  
  # Matrix multiplication
  cosine <- uFunctionValue %*% t(V) #calculate the cosine values between each solution and each vector
  acosine <- acos(cosine)
  
  # call max with argument to give index too..
  list[maxc, maxcidx] <- Max(cosine, 2, T)
  
  # class <- data.frame(c = rep(NA, VN)) #classification
  class <- as.list(rep(NA, VN))
  
  for (i in 1:N){
    # empty at first
    #     if (is.na(class[maxcidx[[i]], 'c'])){
    #       class[maxcidx[[i]], 'c'] <- R(i)
    #     }
    #     else { # append
    #       class[maxcidx[[i]], 'c'] <- R(class[maxcidx[[i]], 'c'], i)
    #     }
    if (is.na(class[maxcidx[i]])){
      class[maxcidx[i]] <- i
    }
    else {
      class[[maxcidx[i]]] <- c(class[[maxcidx[i]]], i)
    }
  }
  
  Selection <- NULL
  for (k in 1:VN){
    if (!is.na(class[k])){
      sub <- class[[k]]
      subFunctionValue <- FunctionValue[sub,]
      
      #APD calculation
      subacosine <- acosine[sub, k]
      subacosine <- subacosine/refV[k]# angle normalization
      D1 <- sqrt(Sum(subFunctionValue^2,2))# Euclidean distance from solution to the ideal point
      D <- D1*(1 + theta0[1]*(subacosine))# APD
      
      list[mind, mindidx] <- Min(D, 1, T)
      Selection <- C(Selection, sub[mindidx])
    }
  }
  return(Selection)
}


#return randomly matched mating pool
# function [MatingPool] <- F_mating(Population)
# Completed !!
F_mating <- function(Population){
  ND <- size(Population)
  N <- ND[1]
  D <- ND[2]
  MatingPool <- zeros(N,D)
  RandList <- R(sample(N))
  MatingPool <- Population[RandList, ]
  if (N %% 2 == 1){
    MatingPool <- C(MatingPool, MatingPool[1,])
  }
  return (MatingPool)
}

#Function to generate uniformly distributed weight vectors
# function [N,W] <- F_weight(p1,p2,M)
# Completed !!
F_weight <- function(p1, p2, M){
  list[N,W] <- T_weight(p1,M)
  if (p2 > 0){
    list[N2,W2] <- T_weight(p2,M)
    N <- N+N2
    W <- C(W, W2*0.5+(1 - 0.5)/(M))
  }
  return (list(N, W))
}

# function [N,W] <- T_weight(H,M)
T_weight <- function(H, M){
  N <- nchoosek(H+M-1,M-1)
  Temp <- nchoosek(1:(H+M-1),M-1)-repmat(0:(M-2),nchoosek(H+M-1,M-1),1)-1
  W <- zeros(N,M)
  W[,1] <- Temp[,1]-0
  if ((M-1)>=2)
    for (i in 2: (M-1)) {
      W[,i] <- Temp[,i]-Temp[,i-1]
    }
  W[,size(W, 2)] <- H-Temp[,size(Temp, 2)]
  W <- W/H
  return (list(N, W))
}


# checked!
P_sort <- function(FunctionValue, operation = ""){
  
  # Efficient non-dominated sort on sequential search strategy, TEVC, 2014,
  # Xingyi Zhang, Ye Tian, Ran Cheng and Yaochu Jin
  # Copyright 2014 BCMA Group, Written by Mr Ye Tian and Prof Xingyi Zhang
  # Contact: xyzhanghust@gmail.com
  
  if (operation == 'half'){
    kinds <- 2
  } else if (operation == 'first'){
    kinds <- 3
  } else {
    kinds <- 1
  }
  
  NM <- size(FunctionValue)
  N <- NM[1]; M <- NM[2]
  MaxFront <- 0
  Sorted <- zeros(1,N)
  
  list[FunctionValue,rank] <- Sortrows(FunctionValue)
  FrontValue <- zeros(1,N) + Inf
  
  while ((kinds == 1 && sum(Sorted)<N) || (kinds == 2 && sum(Sorted)<N/2) || (kinds == 3 && MaxFront<1)){
    MaxFront <- MaxFront + 1
    ThisFront <- as.logical(zeros(1, N))
    for (i in 1:N){
      if (!Sorted[i]){
        x <- 0
        for (j in 1:N){
          if (ThisFront[j]){
            x <- 2
            for (j2 in 2 : M){
              if (FunctionValue[i,j2] < FunctionValue[j,j2]){
                x <- 0
                break
              }
            }
            if (x == 2){
              break
            }
          }
        }
        if (x != 2){
          ThisFront[i] <- T
          Sorted[i] <- T
        }
      }
    }
    # Potentially problematic?
    # index <- 1 * (1 %in% ThisFront)
    FrontValue[rank[ThisFront]] <- MaxFront
  }
  return(list(FrontValue,MaxFront))
}

# checked!
P_generator <- function(MatingPool,Boundary,Coding,MaxOffspring = 0){
  # This function includes the SBX crossover operator and the polynomial
  # mutation operator.
  
  ND <- size(MatingPool)
  N <- ND[1]
  D <- ND[2]
  if (MaxOffspring < 1 || MaxOffspring > N){
    MaxOffspring <- N
  }
  
  if (Coding == 'Real'){
    ProC <- 1 # SBC probability
    ProM <- 1/D # Polynomial mutation probability
    DisC <- 30 # SBX distribution index
    DisM <- 20 # Polynomial mutation distribution index 
    Offspring <- zeros( N, D)
    
    for (i in seq(1, N, 2) ) {
      beta <- zeros(1, D)
      miu  <- rand(1, D)
      beta[miu<=0.5] <- (2*miu[miu<=0.5])^(1/(DisC+1))
      beta[miu>0.5]  <- (2-2*miu[miu>0.5])^(-1/(DisC+1))
      beta <- beta * (-1)^round(rand(1, D))
      beta[rand(1, D)>ProC] <- 1
      Offspring[i,] <- (MatingPool[i,] +MatingPool[i+1,])/2 + beta * (MatingPool[i,]-MatingPool[i+1,])/2
      Offspring[i+1,] <- (MatingPool[i,]+MatingPool[i+1,])/2 - beta * (MatingPool[i,]-MatingPool[i+1,])/2
    }
    
    Offspring <- Offspring[1:MaxOffspring,]
    
    if (MaxOffspring == 1){
      MaxValue <- Boundary[1,]
      MinValue <- Boundary[2,]
    } else {
      # repmat defined at utils
      MaxValue <- repmat(Boundary[1,], MaxOffspring,1);
      MinValue <- repmat(Boundary[2,], MaxOffspring,1);
    }
    k    <- rand(MaxOffspring, D)
    miu  <- rand(MaxOffspring, D)
    Temp <- (k<=ProM & miu<0.5)
    
    Offspring[Temp] <- Offspring[Temp]+(MaxValue[Temp]-MinValue[Temp]) * ((2*miu[Temp]+(1-2*miu[Temp])*(1-(Offspring[Temp]-MinValue[Temp])/(MaxValue[Temp]-MinValue[Temp]))^(DisM+1))^(1/(DisM+1))-1)
    Temp <- (k<=ProM & miu>=0.5)
    Offspring[Temp] <- Offspring[Temp]+(MaxValue[Temp]-MinValue[Temp]) * (1-(2*(1-miu[Temp])+2*(miu[Temp]-0.5)*(1-(MaxValue[Temp]-Offspring[Temp])/(MaxValue[Temp]-MinValue[Temp]))^(DisM+1))^(1/(DisM+1)))
    
    Offspring[Offspring>MaxValue] <- MaxValue[Offspring>MaxValue]
    Offspring[Offspring<MinValue] <- MinValue[Offspring<MinValue]
    
  }
  return (Offspring)
}