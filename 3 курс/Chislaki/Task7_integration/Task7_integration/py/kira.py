import math
import numpy as np
import matplotlib.pyplot as plt


def f(x):
    A = []
    for i in range(x.__len__()):
        A.append(
            4.5 * math.cos(7 * x[i]) * math.exp(-2 / 3 * x[i]) + 1.4 * math.sin(1.5 * x[i]) * math.exp(-x[i] / 3) + 3)
    return A


def mu(n, a, b, alpha):
    B = np.array([])
    for j in range(n):
        B = np.append(B, ((b - 2.1) ** (j - alpha + 1) - (a - 2.1) ** (j - alpha + 1)) / (j - alpha + 1))
    return B


def X(x):
    size = x.__len__()
    ans = np.zeros((size, size))
    ans[0][:] = np.ones(size)
    for i in range(1, size):
        for j in range(size):
            ans[i][j] = ans[i - 1][j] * x[j]
    return ans


count = 8


def IKF(a, b, n):
    alpha = 2 / 5

    x = np.linspace(a - 2.1, b - 2.1, n)

    A = np.linalg.solve(X(x), mu(n, a, b, alpha))
    F = f(np.linspace(a, b, n))
    ans = 0
    for i in range(F.__len__()):
        ans += F[i] * A[i]

    return ans


def GaussKF(a, b, n):
    alpha = 2 / 5

    MU_j = mu(2 * n, a, b, alpha)
    MU_js = np.zeros((n, n))
    MU_jn = np.zeros(n)
    for s in range(n):
        for j in range(n):
            MU_js[s][j] = MU_j[j + s]
        MU_jn[s] = -MU_j[n + s]

    a_j = np.linalg.solve(MU_js, MU_jn)
    a_j = np.append(a_j, 1)[::-1]

    x = np.roots(a_j)

    A = np.linalg.solve(X(x), MU_j[0: n])

    x = x + ([2.1] * x.__len__())
    F = f(x)

    ans = 0
    for i in range(F.__len__()):
        ans += F[i] * A[i]

    return ans


def SKF(a, b, k):
    space = np.linspace(a, b, k)
    ans = 0
    for i in range(space.__len__() - 1):
        if i == 0:
            ans = IKF(space[i], space[i + 1], k)
        elif i % 2 == 0:
            ans += IKF(space[i], space[i + 1], k)
        else:
            ans += GaussKF(space[i], space[i + 1], k)
    return ans


def parametersOfMethod(h, L, a, b, func, Epsilon):
    h1 = h
    h2 = h / L
    h3 = h2 / L

    n1 = math.ceil((b - a) // h1)
    n2 = math.ceil((b - a) // h2)
    n3 = math.ceil((b - a) // h3)

    S1 = func(a, b, n1)
    S2 = func(a, b, n2)
    S3 = func(a, b, n3)

    # Ричардсон
    # print((S3 - S2) / (S2 - S1))
    m = -math.log((S3 - S2) / (S2 - S1)) / math.log(L)
    # Рунге
    Rh1 = (S2 - S1) / (1 - L ** (-m))
    Rh2 = (S2 - S1) / (L ** m - 1)
    hOpt = h1 * ((Epsilon) / abs(Rh1)) ** (1 / m)
    nOpt = math.ceil((b - a) / hOpt)

    return np.array([Rh1, Rh2, m, hOpt, nOpt])


a = 2.1
b = 3.3
c = 2.8
h = 0.2
L = 2
Epsilon = 0.00001

k = int(parametersOfMethod(h, L, a, b, IKF, Epsilon)[4])
print(IKF(a, b, k))