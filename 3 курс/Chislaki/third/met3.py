import numpy as np
from enum import Enum
from scipy import linalg

class ApproxType(Enum):
    algebraic = 0
    legendre = 1
    harmonic = 2


def func(x):
    """
    this method should implement VECTORIZED target function
    """
    return x*np.tan(x)

def pol (a,x):
    return sum([a[i]*x**i for i in range(len(a))])

def pol_leg(n,x): #pol_leg(n+1) 
    if n == -1:
        return 1
    elif n == 0:
        return x
    else:
        return ((2*n + 1)/(n + 1))*x*pol_leg(n-1,x) - n/(n+1)*pol_leg(n-2,x)

def pos_pol_leg(a,x):
    return sum([a[i]*pol_leg(i-1,x) for i in range(len(a))])

def approx(X0, Y0, X1, approx_type: ApproxType, dim):
    """
    this method should perform approximation on [-1; 1] interval
    :param X0: X-values (1 x N0)
    :param Y0: Y-values (1 x N0)
    :param X1: approximation points (1 x N1)
    :param approx_type:
        0 - algebraic polynomes (1, x, x^2, ...)
        1 - legendre polynomes
        2 - harmonic
    :param dim: dimension
    :return Y1: approximated Y-values (1 x N1)
    :return a: vector (1 x dim) of approximation coefficients
    :return P: (for approx_type 0 and 1) coefficients of approximation polynome P (1 x dim)
    """
    if approx_type is ApproxType.algebraic:
        Q = np.array([[i**k for k in range(dim + 1)] for i in X0 ])
        H = Q.T @ Q
        b = Q.T @ Y0
        ans = linalg.solve(H,b)
        return pol(ans,X1), \
               np.eye(1, dim)[0], \
               ans
    
    if approx_type is ApproxType.legendre:
        Q = np.array([[pol_leg(k,i) for k in range(-1,dim)] for i in X0 ])
        H = Q.T @ Q
        b = Q.T @ Y0
        ans = linalg.solve(H,b)

        G = np.eye(dim+1)
        for i in range(1,dim):
            G[i+1][0] = -i/(i + 1)*G[i - 1][0]
            for j in range(dim):
                G[i+1][j+1] = (2*i + 1)/(i + 1)*G[i][j] -i/(i + 1)*G[i - 1][j + 1]

        return pos_pol_leg(ans,X1), \
               np.eye(1, dim)[0], \
               ans @ G


