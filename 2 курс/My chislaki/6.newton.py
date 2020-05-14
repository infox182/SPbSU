import math


e = 0.0001

def f1(x, y):
    return math.sin(y + 0.5) - x - 1

def f2(x, y):
    return y + math.cos(x - 2)

def gn(x,y):
    return (x + 1 -math.sin(y + 0.5) + math.cos(y + 0.5)*(y + math.cos(x - 2)))/(math.cos(y + 0.5)*math.sin(x - 2) -1)

def hn(x, y):
    return - y - math.cos(x - 2) + math.sin(x - 2)*gn(x,y) 

x = 0
y = 0
xk = 0
yk = 0

while math.sqrt(math.pow(f1(x, y), 2) + math.pow(f2(x,y), 2)) > e:
    xk = x + gn(x, y)
    yk = y + hn(x, y)
    x = xk
    y = yk

print (x, y)
