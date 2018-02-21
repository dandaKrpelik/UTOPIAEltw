#taget vector - linear combinaion
# score  - use some (VADS, weighted minmax)

#source('./DTLZ.R')

getInitPop <-function(dim, npop)
{
  pop = matrix( runif(dim*npop)*10 ,nrow=npop, ncol=dim)
}

cross <- function(pop, F_par = 0.7, C_par = 0.7)
{
  npop = nrow(pop)
  dim = ncol(pop)
  
  new_pop = matrix(,nrow=2*npop, ncol = dim)
  new_pop[1:npop,] = pop
  
  for(i in 1:npop)
  {
    permut = sample(seq(1,npop))
    o_ind = permut[1:3]
    for(j in 1:3)
    {
      if (o_ind[j] == i)
      {
        o_ind[j] = permut[4]
      }
    }
    
    trial_chrom = F_par*(pop[o_ind[1],] - pop[o_ind[2],]) + pop[o_ind[3],]
    posi = sample(1:dim)[1]
    
    offspring = pop[i,]
    
    goOn = TRUE
    v = 1
    r = runif(1)
    while(goOn)
    {
      offspring[posi] = min(abs(trial_chrom[posi]),10)     ### UGLY BOUNDS
      
      posi = posi +1
      v = v+1
      if(posi > dim) posi = 1
      if(r > C_par ** (v-1)) goOn = FALSE
    }
    
    new_pop[npop+i,] = offspring
    
  }
  return(new_pop)
}

morder <- function(a)
{
  n = length(a)
  o = order(a)
  out = rep(0,n)
  for (i in 1:n)
  {
    out[o[i]] = i
  }
  return(out)
}


solveProblem <- function(fobj, dim, npop, maxIter = 50, maxMe=1, no_objs=1, no_T=1)
{
  
  target = runif(no_objs*no_T)
  target = matrix(target, nrow=no_T)
  
  pop = getInitPop(dim, npop)
  
  for (iteration in 1:maxIter)
  {
      
      new_pop = cross(pop)
      
      fit = matrix(,nrow=2*npop, ncol=no_objs)
      for (i in 1:(2*npop))
      {
        y = maxMe*fobj(new_pop[i,])
        fit[i,] = y
      }
      
      S = matrix(,nrow = 2*npop, ncol = no_T)
      for (i in 1:(2*npop))
      {
        for (j in 1:no_T)
        {
          S[i,j] = sum(fit[i,]*target[j,])
        }
      }
      
      R = matrix(, nrow = 2*npop, ncol = no_T)
      for(k in 1:no_T)
      {
        R[,k] = order(S[,k])
      }
    
      hold_fit = fit
      for (i in 1:npop)
      {
        pop[i,] = new_pop[R[i,1], ]
        fit[i] = hold_fit[R[i,1]]
      }
  }  
  
  fit = fit[1:npop]
  return(list(pop,maxMe*fit))
}


f9 <- function(x)
{
  xc = x[1]
  xe = x[2]
  
  val1 = 3-0.2*xc+0.3*xe
  val2 = 3+0.2*xc-0.1*xe
  
  return(min(val1,val2))
}

f10 <- function(x)
{
  xc = x[1]
  xe = x[2]
  
  #print(c(xc,xe))
  val = sin(xc-xe)/sqrt(xc**2+xe**2)
  return(val)
}


f9minmax <- function(x)
{
  xc = x[1]
  f9max <- function(y)
  {
    xe = y[1]
    return (f9(c(xc,xe)))
  }
  
  L = solveProblem(f9max,1,5,15,-1)
  fit = L[[2]]
  return(max(fit))
}



f10minmax <- function(x)
{
  xc = x[1]
  f10max <- function(y)
  {
    xe = y[1]
    return (f10(c(xc,xe)))
  }
  
  L = solveProblem(f10max,1,10,25,-1)
  fit = L[[2]]
  return(max(fit))
}

getArgOpt <- function(f, z, ind, npop, maxit, maxMe)
{
  x0 = rep(z,2)
  objf <- function(x)
  {
    xx = x0
    xx[ind] = x[1]
    return (f(xx))
  }
  
  L = solveProblem(objf, 1, npop, maxit, maxMe)
  fit = L[[2]]
  
  #return(maxMe*min(fit))
  return(L[[1]][ order(fit)[1],  ][1])
}


getArgMin <- function(f, npop, maxit)
{
  L = solveProblem(f,2,npop,maxit)
  fit = L[[2]]
  return(L[[1]][ order(fit)[1],  ][1]) 
}


testMinMax <- function(f,npop, maxit)
{
  x_c = runif(1)*10
  x_emax = getArgOpt(f, x_c, 2 , npop, maxit, -1)
  A = c(x_emax)
  A_C = c()
  
  for (iteration in 1:maxit)
  {
    nfnc <- function(x)
    {
      sols = rep(0,iteration)
      for (i in 1:iteration)
      {
        sols[i] = f(c(x[1],A[i]))
      }
      return (max(sols,na.rm=TRUE))
    }
    
    L = solveProblem(nfnc, 1, npop, maxit)
    FI = L[[2]]
    best = order(FI)[1]
    x_cmin = L[[1]][best][1]
    
    
    #maximum in A
    xes = rep(0,iteration)
    for (i in 1:iteration)
    {
      xes[i] = f(c(x_cmin,A[i]))
    }
    #argxemax = A[ order(xes)[1]  ]
    xemax = max(xes, na.rm = TRUE)
    
    
    
    #opti Maximum
    x_emax = getArgOpt(f, x_cmin, 2 , npop, maxit, -1)
    fx_emax = f(c(x_cmin, x_emax))
    
    #print("HERE")
    #print(x_cmin)
    #print(x_emax)
    #print(xemax)
    #print(xes)
    #print(fx_emax)
    
    if(xemax < fx_emax)
    {
      A = append(A, x_emax)
      A_C = rbind(A_C, c(x_cmin, x_emax, fx_emax))
    }
    #print(A)
  }
 
  KK = nrow(A_C)
  for (i in KK)
  {
    for(j in KK)
    {
      tval = f(c(A_C[i,1],A_C[j,2]))
      if (tval > A_C[i,3])
      {
        A_C[i,2] = A_C[j,2]
        A_C[i,3] = tval
      }
    }
  }
  
   
  o = order(A_C[,3])
  
  return(A_C[o[1],]) 
}
