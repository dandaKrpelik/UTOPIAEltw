rdir0 <- function(p)
{
	n = length(p)
	out = rep(0,n)
	for (i in 1:n)
	{
		out[i] = rgamma(1, p[i], 1)
	}
	sm = sum(out)
	out = out / sm
	return(out)
}

rdirichlet <- function(n, p)
{
	k = length(p)
	p = rep(2/k,k)
	out = c()
	for (i in 1:n)
	{
		out = rbind(out, rdir0(p))
	}	
	return(out)
}


sampleQ <- function(n,k)
{
	return(rdirichlet(n, rep(2/k,k)))
}


weight <- function(x, p_target)
{
	n = length(x)
	p_q = rep(2/n,n)
	out = 1
	for (i in 1:n)
	{
		out = out * x[i] ** (p_target[i] - p_q[i])
	}
	return(out)
}

weightT <- function(x, t_target)
{
	p_target = 2*t_target
	return(weight(x, p_target))
}


effSS <- function(weights)
{
	n = length(weights)
	nom = 0
	den = 0

	for (i in 1:n)
	{
		nom = nom + weights[i]
		den = den + weights[i] ** 2
	}
	nom = nom ** 2
	return (nom / den)
}


objFun <- function(x)
{
	return(x[1] + 2*x[2] + 5*x[3] + 4*x[4] - 3*x[5])
}

ISest <- function(x, ts)
{
	n = nrow(x)
	w = apply(x,1,function(z){weightT(z,ts)})
	
	nom = 0
	for(i in n)
	{
		nom = nom + w[i]*objFun(x[i,])
	}
	
	return (nom/sum(w))
}


plotESS <- function(ns)
{
	tot = length(ns)
	y = rep(0,tot)
	for (i in 1:tot)
	{
		n = ns[i]
		samples = sampleQ(n,5)
		w = apply(samples, 1 , function(z){weightT(z,c(0.6,0.1,0.1,0.1,0.1))})
		y[i] = effSS(w)
		
	}
	plot(ns,y)
}


findMin <- function(xs, eps = 1e-4)
{
	n = nrow(xs)
	k = ncol(xs)

	t0 = rep(1/k,k)

	#equality approx
	ui = rep(1,k)
	ci = c(1-eps, -(1+eps))
	ui = rbind(ui, rep(-1,k))

	#lower bound
	ui = rbind(ui, diag(k))
	ci = append(ci, rep(0.1,k))

	#upper bound
	ui = rbind(ui, -diag(k))
	ci = append(ci, -rep(0.6,k))


	print(ui)
	print(ci)

	fnc <- function(t){ISest(xs, t)}
	res = constrOptim(t0, fnc, NULL, ui,ci, control = list(maxit = 1000))
	return (res)
}

samples = sampleQ(50,5)
