import numpy as np
from scipy import linalg

# 13 function
def func(x):
    """
    this method should implement VECTORIZED target function
    """
    return x*np.log1p(x)


def interpol(X, Y):
    """
    this method should find polynomial interpolation
    :param X: X-values (1xN)
    :param Y: Y-values (1xN)
    :return: coefficients of N-1-degree polynome P (1xN)
    """
    A = np.array([[i**k for k in range(len(X))] for i in X ])
    b = np.array(Y)
    ans = linalg.solve(A,b)
    
    return np.flip(ans)
