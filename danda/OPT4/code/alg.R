#taget vector - linear combinaion
# score  - use some (VADS, weighted minmax)

#source('./DTLZ.R')

getInitPop <-function(dim, npop)
{
  pop = matrix( runif(dim*npop) ,nrow=npop, ncol=dim)
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
      offspring[posi] = abs(trial_chrom[posi])
      
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

mtruncate <- function(R_in)
{
  R = R_in
  n = nrow(R)
  nt = ncol(R)
  
  desired_n = n/2
  
  passedIndeces = c()
  
  n_passed = 0
  for(i in 1:nt)
  {
    dn = desired_n - n_passed
    ord = order(R[,i])
    
    minsame = dn
    maxsame = dn
    vallim = R[ord[dn],i]
    
    if (dn > 1)
    {
      for (j in (dn-1):1)
      {
        if (R[ord[j],i] < vallim){ break }
        minsame = j
      }
    }
      
    if (dn < n)
    {
      for (j in (dn+1):n)
      {
        if (R[ord[j],i] > vallim){ break }
        maxsame = j
      }
    }
    
    if (minsame > 1)
    {
      passedIndeces = append(passedIndeces, ord[1:(minsame-1)])
     M = matrix(rep(2*n,(minsame-1)*(nt-i+1) ) ,nrow=minsame-1, ncol = nt-i+1)
      R[ord[1:(minsame-1)],i:nt] = M
    }
    
    n_passed = length(passedIndeces)
    if (minsame == maxsame){ break }
  }
  
  #pi2 = passedIndeces
  
  if(n_passed < desired_n){
    for (i in minsame:(desired_n - n_passed))
    {
      passedIndeces = append(passedIndeces, ord[i])
      #pi2 = append(pi2, ord[i])
    }
  }
  #return(list(passedIndeces,pi2))
  return(passedIndeces)
}


getWrapDTLZ <- function(x, dim, no_objs)
{
  return (function(x){DTLZ1(x,no_objs, dim)})
}


solveProblem <- function(fobj, dim, npop, no_objs, no_T, maxit=25)
{
  
  target = runif(no_objs*no_T)
  target = matrix(target, nrow=no_T)
  
  pop = getInitPop(dim, npop)
  
  
  for(iterate in 1:maxit)
  {
    
    new_pop = cross(pop)
    
    fit = matrix(,nrow=2*npop, ncol=no_objs)
    for (i in 1:(2*npop))
    {
      y = fobj(new_pop[i,])
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
      R[,k] = morder(S[,k])
    }
    #RR = matrix(, nrow = 2*npop, ncol = no_T)
    for(k in 1:(2*npop))
    {
      R[k,] = sort(R[k,])
    }
    
    npop_indeces = mtruncate(R)
    
    hold_fit = fit
    
    for (i in 1:npop)
    {
      pop[i,] = new_pop[npop_indeces[i],]
    }
  }
  return(pop)
}





