import math

def func(x): 
	a = math.tan(0.5*x + 0.2) - x**2
	return a
	
def der(x):
	a = 0.5*(1/((math.cos(0.5*x + 0.2))**2)) - 2*x
	return a

e = 0.0001
t = 100
x0 = 1

while t > e:
	xk = x0 - func(x0)/der(x0)
	t = abs(x0 - xk)
	x0 = xk	

print(x0)
