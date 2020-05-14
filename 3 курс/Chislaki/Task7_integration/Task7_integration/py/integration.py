import numpy as np
from scipy import integrate
from scipy.special import comb

MAXITER = 12

def moments(max_s: int, xl: float, xr: float, a: float = 0.0, b: float = 1.0, alpha: float = 0.0, beta: float = 0.0):

    assert alpha * beta == 0, \
        f'alpha ({alpha}) and/or beta ({beta}) must be 0'

    if alpha == 0 and beta == 0:
        return [(xr ** s - xl ** s) / s for s in range(1, max_s + 2)]

    mu = np.zeros(max_s + 1)
    gamma = alpha
    k = 1
    c = a
    ot = 1

    if alpha == 0:
        k = -1
        gamma = beta
        c = b
        ot = -1
        xr,xl = xl,xr

    for j in range(max_s + 1):

        coefs = np.array([k** (j - i) * comb(j, i) * c ** i / (j - i - gamma + 1) for i in range(j + 1)])
        r = np.array([(ot*(xr - c)) ** (j - i - gamma + 1) for i in range(j + 1)])
        l = np.array([(ot*(xl - c)) ** (j - i - gamma + 1) for i in range(j + 1)])

        mu[j] = sum(coefs * (r - l))
    return mu

def quad(f, xl: float, xr: float, nodes, *params):
    """
    small Newton—Cotes formula
    f: function to integrate
    xl : left limit
    xr : right limit
    nodes: nodes within [xl, xr]
    *params: parameters of the variant — a, b, alpha, beta)
    """
    n = len(nodes)
    mu = moments(len(nodes) - 1, xl, xr, *params)
    X = [[x**s for x in nodes] for s in range(n)]
    A = np.linalg.solve(X, mu)
    # raise NotImplementedError
    # return # small formula result over [xl, xr]
    return [f(x) for x in nodes] @ A



def quad_gauss(f, xl: float, xr: float, n: int, *params):
    """
    small Gauss formula
    f: function to integrate
    xl : left limit
    xr : right limit
    n : number of nodes
    *params: parameters of the variant — a, b, alpha, beta)
    """
    mu = (moments(2 * n - 1, xl, xr, *params))
    mu_js = np.zeros((n,n))
    mu_ns = np.zeros(n)
    
    for s in range(n):
        for j in range(n):
            mu_js[j][s] = mu[j+s]
        mu_ns[s] = -mu[n+s]
    
    a_j = np.linalg.solve(mu_js, mu_ns)
    a_j = np.append(a_j, 1)[::-1]

    nodes = np.roots(a_j)
    X = [[x**s for x in nodes] for s in range(n)]
    A = np.linalg.solve(X, mu[:n])
    # raise NotImplementedError
    # return # small formula result over [xl, xr]
    return [f(x) for x in nodes] @ A

def composite_quad(f, xl: float, xr: float, N: int, n: int, *params):
    """
    composite Newton—Cotes formula
    f: function to integrate
    xl : left limit
    xr : right limit
    N : number of steps
    n : number of nodes od small formulae
    *params: parameters of the variant — a, b, alpha, beta)
    """
    mesh = np.linspace(xl, xr, N + 1)
    return sum(quad(f, mesh[i], mesh[i + 1], equidist(n, mesh[i], mesh[i + 1]), *params) for i in range(N))


def composite_gauss(f, a: float, b: float, N: int, n: int, *params):
    """
    composite Gauss formula
    f: function to integrate
    xl : left limit
    xr : right limit
    N : number of steps
    n : number of nodes od small formulae
    *params: parameters of the variant — a, b, alpha, beta)
    """
    mesh = np.linspace(a, b, N + 1)
    return sum(quad_gauss(f, mesh[i], mesh[i + 1], n, *params) for i in range(N))


def equidist(n: int, xl: float, xr: float):
    if n == 1:
        return [0.5 * (xl + xr)]
    else:
        return np.linspace(xl, xr, n)


def runge(s1: float, s2: float, L: float, m: float):
    """ estimate m-degree error for s2 """
    # raise NotImplementedError
    return abs(s2-s1)/(L**m - 1)


def aitken(s1: float, s2: float, s3: float, L: float):
    """
    estimate convergence order
    s1, s2, s3: consecutive composite quads
    return: convergence order estimation
    """
    temp1 = s3-s2 
    temp2 = s2-s1
    if temp1*temp2 < 0:
       return -1
    else:
        return -np.log(temp1/temp2)/np.log(L)
    # raise NotImplementedError


def doubling_nc(f, xl: float, xr: float, n: int, tol: float, *params):
    """
    compute integral by doubling the steps number with theoretical convergence rate
    f : function to integrate
    xl : left limit
    xr : right limit
    n : nodes number in the small formula
    tol : error tolerance
    *params : arguments to pass to composite_quad function
    """

    iter = 0
    N = 4

    s1 = composite_quad(f, xl, xr, N, n, *params)
    s2 = composite_quad(f, xl, xr, N*2, n, *params)
    err = runge(s1, s2, 2, 3)
    
    while iter < MAXITER:
        if err <= tol:
            return N*2, s2, err

        N *= 2
        s1 = s2
        s2 = composite_quad(f, xl, xr, 2 * N, n, *params)
        err = runge(s1, s2, 2, 3)
        iter += 1

    if iter == MAXITER:
        print("Convergence not reached!")
        return 0, 0, 10*tol, -100


def doubling_nc_aitken(f, xl: float, xr: float, n: int, tol: float, *params):
    """
    compute integral by doubling the steps number with Aitken estimation of the convergence rate
    f : function to integrate
    xl : left limit
    xr : right limit
    n : nodes number in the small formula
    tol : error tolerance
    *params : arguments to pass to composite_quad function
    """
    # required local variables to return
    # S : computed value of the integral with required tolerance
    # N : number of steps for S
    # err : estimated error of S
    # m : estimated convergence rate by Aitken for S
    # iter : number of iterations (steps doubling)
    iter = 0
    N = 4
    
    s1 = composite_quad(f,xl,xr,N,n,*params)
    s2 = composite_quad(f,xl,xr,N*2,n,*params)
    s3 = composite_quad(f,xl,xr,N*4,n,*params)
    m = aitken(s1,s2,s3,2)
    err = runge(s2,s3,2,m)
    
    while iter < MAXITER:
        if err <= tol:
            return N*4, s3, err, m
        N *= 2
        s1 = s2
        s2 = s3
        s3 = composite_quad(f,xl,xr,N*4,n,*params)
        m = aitken(s1,s2,s3,2)
        err = runge(s2,s3,2,m)
        iter += 1
    if iter == MAXITER:
        print("Convergence not reached!")
        return 0, 0, 10*tol, -100


def doubling_gauss(f, xl: float, xr: float, n: int, tol: float, *params):
    """
    compute integral by doubling the steps number with theoretical convergence rate
    f : function to integrate
    xl : left limit
    xr : right limit
    n : nodes number in the small formula
    tol : error tolerance
    *params : arguments to pass to composite_quad function
    """
    # required local variables to return
    # S : computed value of the integral with required tolerance
    # N : number of steps for S
    # err : estimated error of S
    # iter : number of iterations (steps doubling)
    iter = 0
    N = 4
    
    s1 = composite_gauss(f,xl,xr,N,n,*params)
    s2 = composite_gauss(f,xl,xr,N*2,n,*params)
    err = runge(s1,s2,2,3)

    while iter < MAXITER:
        if err <= tol:
            return N*2, s2, err 
        N *= 2
        s1 = s2
        s2 = composite_gauss(f, xl, xr,N*2, n, *params)
        err = runge(s1,s2,2,3)
        iter += 1
    if iter == MAXITER:
        print("Convergence not reached!")
        return 0, 0, 10*tol


def doubling_gauss_aitken(f, xl: float, xr: float, n: int, tol: float, *params):
    """
    compute integral by doubling the steps number with Aitken estimation of the convergence rate
    f : function to integrate
    xl : left limit
    xr : right limit
    n : nodes number in the small formula
    tol : error tolerance
    *params : arguments to pass to composite_quad function
    """
    # required local variables to return
    # S : computed value of the integral with required tolerance
    # N : number of steps for S
    # err : estimated error of S
    # m : estimated convergence rate by Aitken for S
    # iter : number of iterations (steps doubling)
    iter = 0
    N = 4
    
    s1 = composite_gauss(f,xl,xr,N,n,*params)
    s2 = composite_gauss(f,xl,xr,N*2,n,*params)
    s3 = composite_gauss(f,xl,xr,N*4,n,*params)
    m = aitken(s1,s2,s3,2)
    err = runge(s2,s3,2,m)
    
    while iter < MAXITER:
        if err <= tol:
            return N*4, s3, err, m
        N *= 2
        s1 = s2
        s2 = s3
        s3 = composite_gauss(f,xl,xr,N*4,n,*params)
    
        m = aitken(s1, s2, s3, 2)
        
        err = runge(s2,s3,2,m)
        iter += 1
    
    if iter == MAXITER:
        print("Convergence not reached!")
        return 0, 0, 10*tol, -100


def optimal_nc(f, xl: float, xr: float, n: int, tol: float, *params):
    """ estimate the optimal step with Aitken and Runge procedures
    f : function to integrate
    xl : left limit
    xr : right limit
    n : nodes number in the small formula
    tol : error tolerance
    *params : arguments to pass to composite_quad function
    """
    # required local variables to return
    # S : computed value of the integral with required tolerance
    # N : number of steps for S
    # err : estimated error of S
    # iter : number of iterations (steps doubling)
    
    N = 4
    lens = abs(xr - xl)
    h = lens/4
    s1 = composite_quad(f,xl,xr,N,n,*params)
    s2 = composite_quad(f,xl,xr,N*2,n,*params)
    s3 = composite_quad(f,xl,xr,N*4,n,*params)
    m = aitken(s1,s2,s3,2)

    Rh = runge(s1,s2,2,m)
    h_optim = 0.8*h*(tol/abs(Rh))**(1/m)
    N = int(np.ceil(lens/h_optim))
    
    s1 = composite_quad(f, xl, xr, N, n, *params)
    s2 = composite_quad(f, xl, xr, N * 2, n, *params)
    s3 = composite_quad(f, xl, xr, N * 4, n, *params)
    m = aitken(s1, s2, s3, 2)
    err = runge(s1, s2, 2, m)

    return N,s2,err
