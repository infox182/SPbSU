import math
import numpy as np

e = 0.0000001
A = [[8., 1., 1., 1.], [1., 10., 1., 1.], [1., 1., 12., 1.], [1., 1., 1., 14.]]
b = [11., 13., 15., 17.]

def sum(ar, val):
    a = [0, 0, 0, 0]
    for i in range(len(ar)):
        a[i] = ar[i] + val
    return a

def diff(ar1, ar2):
	a = [0, 0, 0, 0]
	for i in range(len(ar1)):
		a[i] = ar1[i] - ar2[i]
	return a

def mult(ar, val):
	a = [0, 0, 0, 0]
	for i in range(len(ar)):
		a[i] = ar[i]*val
	return a

def ln(a):
    rez = 0
    for i in range(len(a)):
        rez += a[i]**2
    return (math.sqrt(rez))

def div(ar, val):
	a = [0, 0, 0, 0]
	for i in range(len(ar)):
		a[i] = ar[i]/val
	return a

def vmult(a, b):
    rez = 0
    for i in range (len(a)):
        rez += a[i]*b[i]
    return rez

for i in range(len(A)):
    b[i] = b[i]/A[i][i]
    A[i] = div(A[i], A[i][i])
    A[i][i] -= 1
    A[i] = mult(A[i], -1)

x0 = sum(b, 0)
t = 100
xk = [0, 0, 0, 0]

while t > e:
    xk[0] = vmult(x0, A[0]) + b[0]
    xk[1] = xk[0]*A[1][0] + x0[2]*A[1][2] + x0[3]*A[1][3] + b[1]
    xk[2] = xk[0]*A[2][0] + xk[1]*A[2][1] + x0[3]*A[2][3] + b[2]
    xk[3] = xk[0]*A[3][0] + xk[1]*A[3][1] + xk[2]*A[3][2] + b[3]
    t =  ln(diff(x0, xk))
    x0 = sum(xk, 0)

print(xk)
