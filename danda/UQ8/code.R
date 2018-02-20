rdirichlet <- function(n, p)
{
	out = rgamma(n, p, 1)
	sm = sum(out)
	out = out / sm
	return(out)
}


sampleQ <- function(n,k)
{
	p = rep(2/k,k)
	out = c()
	for (i in 1:n)
	{
		out = rbind(out, rdirichlet(k,p))
	}	
	return(out)
}
