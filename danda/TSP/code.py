import numpy as np

def readData():
	out = []
	with open("att532.tsp",'r') as file:
		for line in file:
			spl = line.split()
			out.append([int(spl[1]),int(spl[2])])

	return out

def L1norm(x,y):
	out = 0
	dim = len(x)
	for i in range(dim):
		out += np.abs(x[i]-y[i])
	return out

def NNheuri(data, distance = L1norm):
	n = len(data)
	path = [np.random.randint(n)]
	to_choose = [i for i in range(n)]
	for i in range(1,n):
		#print((i,path[-1]))
		to_choose.remove(path[-1])
		no_left = len(to_choose)
		#dist = np.zeros(no_left)
		closest = 0
		mindist = -1			
		for j in range(no_left):
			d = distance(data[path[-1]], data[to_choose[j]])
			if (mindist < 0 or d < mindist):		
				#print("changing %d for %d" % (closest,j))
				closest = j
				mindist = d

		#print(j)
		path.append(to_choose[j])

	return path
	
def switch


def iterate(data, path):
	
	
def calcTotalDistance(data,path):
	n = len(data)
	tot_dist = L1norm(data[0],data[-1])
	for i in range(1,n):
		tot_dist += L1norm(data[i-1],data[i])
	return tot_dist	

if __name__ == '__main__':
	data = readData()

	#data_toshuffle = [[x[0],x[1]]]

	path = NNheuri(data)

	n = len(data)
	print(calcTotalDistance(data, path))
	print(calcTotalDistance(data, [i for i in range(n)]))
