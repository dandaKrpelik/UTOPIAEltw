M = matrix(c(1,4,1,1,1,4,4,1,1)/6, nrow = 3, byrow = TRUE)
x0 = c(1,1,1)/3


T = function(state)
{
  return (M %*% state)
}

Tmin = function(state, eps)
{
  out = (1-eps)*T(state) + eps*min(state)
  return (out)
}


Tmax = function(state, eps)
{
  out = (1-eps)*T(state) + eps*max(state)
  return (out)
}

iterate <- function(n,x0, eps = 0.2, k =1)
{
  states = rbind(x0)
  statesMin = rbind(x0)
  statesMax = rbind(x0)
  
  for (i in 1:n)
  {
    ind = nrow(states)
    x = t(T(states[ind,]))
    x1 = t(Tmin(statesMin[ind,],eps))
    x2 = t(Tmax(statesMax[ind,], eps))
    states = rbind(states,x)
    statesMin = rbind(statesMin,x1)
    statesMax = rbind(statesMax,x2)
  }
  
  plot(0:n, states[,k], ylim = c(0,1), main=paste(k))
  lines(0:n, statesMin[,k])
  lines(0:n, statesMax[,k])
  
}

mnorm = function(a)
{
  return(sqrt(sum(a**2)))
}

getLimits = function(eps)
{
  conveps = 1e-9
  
  x0 = c(1,0,0)
  x01 = x0
  
  cnt = 0
  while(TRUE)
  {
    cnt = cnt + 1
    x = Tmin(x0, eps)
    x1 = Tmax(x01, eps)
      
    print(mnorm(x-x0))
    if((mnorm(x-x0) < conveps) & (cnt > 5))
    {
      break
    }
    x0 = x
    x01 = x1
  }
  
  return(list(x,x1))
}