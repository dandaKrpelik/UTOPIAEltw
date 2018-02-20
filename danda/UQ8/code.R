rdirichlet <- function(n, p)
{
	out = rgamma(n, p, 1)
	sm = sum(out)
	out = out / sm
	return(out)
}


sampleQ <- function(n)
{
	p = rep(2/n,n)
	out = rdirichlet(n,p)
	return(out)
}

