import numpy as np
import math
A = np.array([[8, 1, 1, 1], [1, 10, 1, 1], [1, 1, 12, 1], [1, 1, 1, 14]])
b = np.array([10, 12, 14, 16])

def mvmult(M, b):
    b1 = np.zeros(len(b))
    for i in range(len(b)):
        b1 += M[:, i]*b[i]
    return b1

def mmult(M, N):
    K1 = np.zeros((len(M), len(M)))
    for i in range(len(M)):
        K1[:, i] = mvmult(M, N[:, i])
    return K1


def ln(a):
    rez = 0
    for i in range(len(a)):
        rez += a[i]**2
    return (math.sqrt(rez))

def vmult(a, b):
    ab = np.zeros((len(a), len(a)))
    for i in range(len(a)):
        for j in range(len(a)):
            ab[i, j] = a[j]*b[i]
    return ab


for i in range(len(A)):
    e = np.zeros(len(A) - i)
    e[0] = 1
    y = A[i:, i]
    alpha = ln(y)
    w = y + alpha*e
    w = w/ln(w)
    U = np.eye(len(w)) - 2*vmult(w, w)
    Z = np.eye(len(A))
    Z[i:, i:] = U
    A = mmult(Z, A)
    b = mvmult(Z, b)

B = A

x4 = b[3]/B[3, 3]
x3 = (b[2] - x4*B[2, 3])/B[2, 2]
x2 = (b[1] - x4*B[1, 3] - x3*B[1, 2])/B[1, 1]
x1 = (b[0] - x4*B[0, 3] - x3*B[0, 2] - x2*B[0, 1])/B[0, 0]
x = [x1, x2, x3, x4]

print (x)
