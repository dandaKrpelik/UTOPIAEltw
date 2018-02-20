library(ReliabilityTheory)

gr = graph.formula(s -- 1 -- 2 -- 3 -- t, 1 -- 4:5 -- 6 -- t,
		s -- 7 -- 8 -- t, s -- 9 -- 10 -- 11 -- t, 7 -- 10 -- 8)

V(gr)$compType[match( c("1","6","11") , V(gr)$name  )] <- "G1" 
V(gr)$compType[match( c("2","3","9") , V(gr)$name  )] <- "G2"
V(gr)$compType[match( c("4","5","10") , V(gr)$name  )] <- "G3"
V(gr)$compType[match( c("7","8") , V(gr)$name  )] <- "G4" 
#V(gr)$compType 

M_vect = c(3,3,3,2)
ctypes= 4

sig <- computeSystemSurvivalSignature(gr)


N1 = c(10,10,15,20)

sampleT1 = rexp(N1[1],0.55)
sampleT2 = rweibull(N1[2],2.2,1.8)
sampleT3 = exp(rnorm(N1[3], 0.4,0.9))
sampleT4 = rgamma(N1[4], scale=0.9, shape=3.2)

prior <- function(t, typ)
{
	out = list(c(0.001,0.999),c(1,100))
	if (typ == 1)
	{
		if (t < 1.5)
		{
			out[[1]][1] = 0.7 #y
			out[[1]][2] = 0.99
			
			out[[2]][1] = 20
			out[[2]][2] = 35
		}else{ 
			out[[1]][1] = 0.3 #y
			out[[1]][2] = 0.5
			
			out[[2]][1] = 5
			out[[2]][2] = 15			
		}
	}
	if (typ == 2)
	{
		if (t < 2)
		{
			out[[1]][1] = 0.6 #y
			out[[1]][2] = 0.8
			
			out[[2]][1] = 5
			out[[2]][2] = 15
		}else{
			out[[1]][1] = 0.1 #y
			out[[1]][2] = 0.3
			
			out[[2]][1] = 1
			out[[2]][2] = 5			
		}
	}
	if (typ == 3)
	{
		if (t < 2.5)
		{
			out[[1]][1] = 0.8 #y
			out[[1]][2] = 0.9
			
			out[[2]][1] = 10
			out[[2]][2] = 20
		}else{
			out[[1]][1] = 0.4 #y
			out[[1]][2] = 0.5
			
			out[[2]][1] = 5
			out[[2]][2] = 10			
		}
	}
	if (typ == 4)
	{
		out[[1]][1] = 0.01 #y
		out[[1]][2] = 0.99
		
		out[[2]][1] = 1
		out[[2]][2] = 30
	}

	return(out)
}





Pli <- function(l, s, n, y0, n0, M)
{
	out = choose(M, l)
	nom = beta(l+n0*y0+s, M-l+n0*(1-y0)+n - s)
	den = beta(n0*y0+s, n0*(1-y0)+n-s)
	return (out*nom/den)
}

Pl <- function(l, s, n , y0, n0 , M)
{
	out = 1
	for (i in 1:ctypes)
	{
		out = out * Pli(l[[i]], s[i], n[i], y0[i], n0[i], M[i])
	}
	return(out)
}



preciseSurvival <- function(t)
{
	s = c(sum(sampleT1 > t),sum(sampleT2 > t),sum(sampleT3 > t),sum(sampleT4 > t))
	
	y0 = rep(0,ctypes)
	n0 = rep(0,ctypes)
	
	for (i in 1:ctypes)
	{
		pr = prior(t,i)
		y0[i] = (pr[[1]][1]+pr[[1]][2])/2
		n0[i] = (pr[[2]][1]+pr[[2]][2])/2
	}
	
	out = 0
	pls = c()
	for (i in 1:nrow(sig))
	{
		#pls = append(pls,Pl(sig[i,1:ctypes], s, N1, y0, n0, M_vect) )
		out = out + sig[i,ctypes+1]*Pl(sig[i,1:ctypes], s, N1, y0, n0, M_vect)
	}
	return (out)
}


