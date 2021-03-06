import numpy as np
import matplotlib.pyplot as plt


def ftest(x,y):
	return np.sin(x)*y**2+(np.exp(y)-1)

def showMe():
	x = np.linspace(-5,5,100)
	y = np.linspace(-5,5,100)
	y = ftest(x,y)


def sample(n):
	s1 = np.random.randn(n,2)
	s2 = np.random.randn(n,2)
	
	return s1,s2

def estimateF(f, x):
	n = len(x)
	y = f(x[:,0],x[:,1])
	Ef = y.mean()
	Vf = y.var()

	return Ef,Vf

def calcSS(f, s1, s2, Ef, Vf):
	sums = [0,0]

	Ef2 = Ef ** 2

	n = len(s1)
	assert n == len(s2), 't'

	for i in range(n):
		x1 = s1[i][0]
		x2 = s1[i][1]
		y1 = s2[i][0]
		y2 = s2[i][1]

		sums[0] += f(x1,x2)*f(x1,y2) - Ef2
		sums[1] += f(x1,x2)*f(y1,x2) - Ef2

	sums[0] /= Vf*n
	sums[1] /= Vf*n

	return sums[0],sums[1], 1- sums[0]-sums[1]

def test(n, f = ftest):
	sam1, sam2 = sample(n)
	Ef, Vf = estimateF(ftest, sam1)
	vals1 = calcSS(ftest, sam1, sam2, Ef, Vf)	

	Ef, Vf = estimateF(ftest, sam2)
	vals2 = calcSS(ftest, sam1, sam2, Ef, Vf)	
	
	return vals1,vals2

def converge(nmax = 150,step = 10, f = ftest):
	ns = np.arange(5,nmax,step)
	y = []

	for i in range(len(ns)):
		n = ns[i]
		v1,v2 = test(n, f)
		y.append([v1,v2])

	y = np.array(y)

	plt.figure()
	plt.plot(ns, y[:,0,0], label = 's1')
	plt.plot(ns, y[:,0,1], label = 's2')
	plt.plot(ns, y[:,0,2], label = 's12')
	plt.title('By guide')
	plt.legend()
	plt.show(0)

	plt.figure()
	plt.plot(ns, y[:,1,0], label = 's1')
	plt.plot(ns, y[:,1,1], label = 's2')
	plt.plot(ns, y[:,1,2], label = 's12')
	plt.title('Mix est')
	plt.legend()
	plt.show(0)

	return ns,y
