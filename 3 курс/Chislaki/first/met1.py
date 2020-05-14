import numpy as np
import random
from operator import itemgetter
import matplotlib.pyplot as plt

def mngs(A, b, x0, eps):
    """
    this method should numerically find min(y),
    where y = 1/2*x.T*A*x + b.T*x
    :param A: matrix NxN
    :param b: matrix Nx1
    :param x0: matrix Nx1
    :param eps: accuracy (see test_met1())
    :return: list of x, list of y
    """
    def func(x):
        return (1 / 2 * x.T @ A @ x + b.T @ x)

    def grad(x):
        return(A @ x + b)

    x_k = x0
    X = [x0]
    Y = [func(x0)]
    
    while np.linalg.norm(A @ x_k + b) > eps :
        q = grad(x_k)
        mu = -(np.dot(q.T, q)/np.dot(q.T, A @ q))
        x_k = x_k + mu*q
        X.append(x_k)
        Y.append(func(x_k))

    return X, Y


def mps(A, b, x0, eps):
    """
    this method should numerically find min(y),
    where y = 1/2*x.T*A*x + b.T*x
    :param A: matrix NxN
    :param b: matrix Nx1
    :param x0: matrix Nx1
    :param eps: accuracy (see test_met1())
    :return: list of x, list of y
    """
    def func(x):
        return (1 / 2 * x.T @ A @ x + b.T @ x)
    def grad(x):
        return(A @ x + b)

    n_dim = b.size
    e_ort = np.eye(n_dim)

    x_k = x0
    X = [x0]
    Y = [func(x0)]

    k = True
    while k:
        for ort in e_ort:
            q = grad(x_k)
            mu = -(np.dot(ort.T, q)/np.dot(ort.T, A @ ort))
            x_k = x_k + mu*ort
            X.append(x_k)
            Y.append(func(x_k))
            if np.linalg.norm(A @ x_k + b) < eps :
                k = False
                break
    return X, Y

# n_dim = 3
# N = 13
# A = np.array([[4,  1,      1],
#                 [1,  6+.2*N, -1],
#                 [1,  -1,      8+.2*N]],
#                 dtype='float')[:n_dim, :n_dim]
# b = np.array([1, -2, 3], dtype='float')[:n_dim]
# x0 = np.array([0, 0, 0], dtype='float')[:n_dim]
# X, Y = mps(A, b, x0 , 1e-6)
# X2,Y2 = mngs(A,b,x0,1e-6)
# print('MPS',X[-1],Y[-1])
# print('MNGS',X2[-1],Y2[-1])
# x1 = np.linalg.solve(A, -b)
# y1 = (1/2 * x1.T @ A @ x1 + b.T @ x1)
# print('REAL',x1,y1)
# eps_y = 1e-6
# print(np.linalg.norm(y1 - Y[-1]) < eps_y)
# print(np.linalg.norm(y1 - Y2[-1]) < eps_y)
# plt.figure()
# plt.xlabel('номер итерации')
# plt.ylabel('точность')
# plt.plot(-np.log10([y - y1 for y in Y]))
# plt.show()