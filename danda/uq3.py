import numpy as np
import matplotlib.pyplot as plt


def sample(n):
	#return np.random.rand(n)*2-1
	return np.linspace(-1,1,n)

def evalPoly(x, phi0, phi_1, i):
	c_ = (i+0.0)/(i+1)
	c0 = (2*i+1.0)/(i+1)

	return c0*x*phi0 - c_*phi_1

def getZ(m, q):
	phi_1 = np.ones(m)
	x = sample(m)
	phi0 = x

	Z = np.zeros((m,q+1))
	Z[:,0] = phi_1
	Z[:,1] = phi0

	for i in range(2,q+1):
		Z[:,i] = evalPoly(x, phi0, phi_1, i-1)
		phi_1 = phi0
		phi0 = Z[:,i]

	return Z, x
	

def solve(x, fnc, Z):
	A = Z.transpose().dot(Z)
	b = Z.transpose().dot(fnc(x))
	
	d = np.diag(A)
	dd = 1/d
	B = np.diag(dd)

	A = B.dot(A)
	b = B.dot(b)

	return np.linalg.solve(A,b)

def f0(x):
	return np.tanh(3*x)

def buildSur(m,q):
	Z,X = getZ(m,q)
	c = solve(X, f0, Z)


	def eval(x):
		out = c[0] + c[1]*x

		p_1 = 1
		p_0 = x
		for i in range(2, q+1):
			p = evalPoly(x, p_0, p_1, i-1)
			p_1 = p_0
			p_0 = p
			out += c[i] * p

		return out

	m = c[0]
	v = 0
	for i in range(1,q+1):
		v += (np.abs(c[i]) **2) * 1./(2*i+1)

	return eval, m, v


def compare(m,q, show = True):
	x = np.linspace(-1,1,50)
	y = f0(x)
	

	sf,m,v = buildSur(m,q)
	yy = sf(x)	

	print('mean est: %f, var est: %f' % (m,v))

	if (show):
		plt.figure()
		plt.plot(x,y)
		plt.plot(x,yy)
		plt.show(0)

	return x,y,yy


def plotConvergenceINM(q, maxm, step = 2):
	minm = max(q+1,2)
	ms = np.arange(minm,maxm,step)
	diffs = []
	for i in range(len(ms)):
		m = ms[i]

		_diffs = []		
		for k in range(50):
			x,y,yy = compare(m,q, False)
			_diffs.append(sum(np.abs(y-yy))/len(x))
		diffs.append(np.array(_diffs).mean())

	plt.figure()
	plt.plot(ms,diffs)
	plt.show(0)
