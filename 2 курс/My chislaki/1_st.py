import math

pogU = 0.00000027
pogV = 0.00000037
pogF = 0.00000033


def geron(w, u):
	return (w + u/w)/2

def f_mysi(x):
	return math.sin(4.5*x + 0.6)/(math.sqrt(1 + x - 12*x**2))

def znamen(x):
	return 1 + x - 12*x**2

def makulaka(x):
	
	sum_sin = 0
	k = 0

	while math.pow(4.5*x + 0.6, 2*k + 1)/math.factorial(2*k + 1) > pogU:
		sum_sin += math.pow(-1, k)*math.pow(4.5*x + 0.6, 2*k + 1)/math.factorial(2*k + 1)
		k+=1

	w = 1
	zn = znamen(x)
	
	while (w - geron(w,zn))/sum_sin > pogF:
		w = geron(w, zn)

	return (sum_sin/(geron(w, zn)))


i = 0.1
print('x          python           makloran')
while i <= 0.2:
	print(round(i, 2), "\t", round(f_mysi(round(i, 2)), 10), "\t", round(makulaka(round(i, 2)), 10))
	i += 0.01





