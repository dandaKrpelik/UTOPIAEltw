import numpy as np
import matplotlib.pyplot as plt

import scipy.interpolate as sin


def fnc(x):
	return 5 + x + np.cos(x)

def scfnc(x):
	M = x[0]
	Y = x[1]
	P = x[2]

	b = 5.
	h = 15.

	return 1-4.*M/(b*h*h*Y) - (P ** 2) / (b*b*h*h*Y*Y)

def sample(n):
	return np.random.randn(n)*4

def sampleuni(n):
	return np.random.rand(n)*30 - 15


def scsample(n):
	y = np.random.randn(n)*0.5+5
	M = np.random.randn(n)*400 + 2000
	P = np.random.randn(n)*100 + 500
	Y = np.exp(y)

	out = []
	for i in range(n):
		out.append([M[i], Y[i], P[i] ])

	return out


def scsampleuni(n):
	Y = np.random.rand(n)*np.exp(7.5)
	M = np.random.rand(n)*6*400 +2000-1200
	P = np.random.rand(n)*6*100 + 500 - 300

	out = []
	for i in range(n):
		out.append([M[i], Y[i], P[i]])
	return out

def RBF(n, suni = sampleuni):
	#samps = sample(n)
	samps = suni(n)

	ys = []
	for i in range(n):
		ys.append(fnc(samps[i]))

	return sin.Rbf(samps, ys), samps


def MCstats(n, f, sam):
	vals = []
	samps = sam(n)

	for i in range(n):
		vals.append(f(samps[i]))

	vals = np.array(vals)

	mean = vals.mean()
	var = vals.var()

	return mean, var

def MCstatsTR(n, f, sam):
	vals = []
	osamps = sam(n)

	samps = []
	for i in range(n):
		if ((osamps[i] > -16) & (osamps[i] < 16)):
			samps.append(osamps[i])

	n = len(samps)

	for i in range(n):
		vals.append(f(samps[i]))

	vals = np.array(vals)

	mean = vals.mean()
	var = vals.var()

	return mean, var

def compare(nsam):

	n = 100
	x = np.linspace(-15,15,n)

	rbf, points = RBF(nsam)

	orig = []
	surr = []

	for i in range(n):
		orig.append(fnc(x[i]))
		surr.append(rbf(x[i]))

	
	svals = []
	for j in range(len(points)):
		svals.append(fnc(points[j]))

	plt.figure()
	plt.plot(x,orig)
	plt.plot(x,surr)
	plt.scatter(points, svals)
	plt.show(0)


def compareStats(n, nsam, sam = sample):
	rbf, points = RBF(nsam)

	m1,v1 = MCstats(n, fnc, sam)
	m2,v2 = MCstats(n, rbf, sam)
	m3,v3 = MCstatsTR(n, rbf, sam)

	print('orinal function: mean = %f, var = %f' % (m1,v1))
	print('surrogate function: mean = %f, var = %f' % (m2,v2))
	print('-||- trunc.: mean = %f, var = %f' % (m3,v3))


	x = np.linspace(-15,15,n)
	orig = []
	surr = []

	for i in range(n):
		orig.append(fnc(x[i]))
		surr.append(rbf(x[i]))

	
	svals = []
	for j in range(len(points)):
		svals.append(fnc(points[j]))

	plt.figure()
	plt.plot(x,orig)
	plt.plot(x,surr)
	plt.scatter(points, svals)
	plt.show(0)


def plotConvergency(maxSam, step = 5, sam = sample):
	MCN = 10000	

	ns = np.round(np.linspace(step, maxSam, maxSam/step))


	meanval1 = []
	minval1 = []
	maxval1 = []

	meanval2 = []
	minval2 = []
	maxval2 = []

	for i in range(len(ns)):
		print('calculating %d of %d' % (i, len(ns)))
		nsam = int(ns[i])


		ms1 = []
		vs1 = []

		ms2 = []
		vs2 = []

		for k in range(25):
			rbf, points = RBF(nsam)

			m1,v1 = MCstats(MCN, fnc, sam)
			m2,v2 = MCstats(MCN, rbf, sam)
			m3,v3 = MCstatsTR(MCN, rbf, sam)		

			
			ms1.append(np.abs(m1 - m2))
			ms2.append(np.abs(m1 - m3))

			vs1.append(np.abs(v1 - v2))
			vs2.append(np.abs(v1 - v3))


		ms1 = np.array(ms1)
		ms2 = np.array(ms2)
		vs1 = np.array(vs1)
		vs2 = np.array(vs2)

		meanval1.append(ms1.mean())
		meanval2.append(ms2.mean())

		m1 = ms1.mean()
		m2 = ms2.mean()

		sd1 = np.mean(np.sqrt(vs1))
		sd2 = np.mean(np.sqrt(vs2))

		minval1.append(m1 - sd1)
		maxval1.append(m1 + sd1)

		minval2.append(m2 - sd2)
		maxval2.append(m2 + sd2)


	
	plt.figure()

	plt.plot(ns, meanval1, label='mean1')
	plt.plot(ns, meanval2, label='meanTR')

	plt.plot(ns,minval1, label='min1')
	plt.plot(ns,maxval1, label='max1')


	plt.plot(ns,minval2, label='minTR')
	plt.plot(ns,maxval2, label='maxTR')

	plt.show(0)

		

