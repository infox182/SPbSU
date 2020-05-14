#deform mnogogr
class vector(object):
    def __init__(self, x1, x2):
        
        self.x1 = x1
        self.x2 = x2

    def __repr__(self):
        return "({0}, {1})".format(self.x1, self.x2)

    def __add__(self, other):
        x1 = self.x1 + other.x1
        x2 = self.x2 + other.x2
        return vector(x1, x2)

    def __sub__(self, other):
        x1 = self.x1 - other.x1
        x2 = self.x2 - other.x2
        return vector(x1, x2)

    def __rmul__(self, other):
        x1 = self.x1 * other
        x2 = self.x2 * other
        return vector(x1, x2)

    def __truediv__(self, other):
        x1 = self.x1 / other
        x2 = self.x2 / other
        return vector(x1, x2)

    def c(self):
        return (self.x1, self.x2)
        
# функция
def f1(point):
    x1, x2 = point
    return 100*(x2 - x1**2)**2 + 5*(1 - x1)**2

def f2(point):
    x1, x2 = point
    return (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2

def mnogogr(f,alpha=1, beta=0.5, gamma=2, maxiter=10):
    
    # инизиализация
    v1 = vector(0, 0)
    v2 = vector(1.0, 0)
    v3 = vector(0, 1)

    for i in range(maxiter):
        adict = {v1:f(v1.c()), v2:f(v2.c()), v3:f(v3.c())}
        points = sorted(adict.items(), key=lambda x: x[1])
        
        b = points[0][0]
        g = points[1][0]
        w = points[2][0]
        
        
        mid = (g + b)/2

        #отражение
        xr = mid + alpha * (mid - w)
        if f(xr.c()) < f(g.c()):
            w = xr
        else:
            if f(xr.c()) < f(w.c()):
                w = xr
            c = (w + mid)/2
            if f(c.c()) < f(w.c()):
                w = c
        if f(xr.c()) < f(b.c()):

            #растяжение
            xe = mid + gamma * (xr - mid)
            if f(xe.c()) < f(xr.c()):
                w = xe
            else:
                w = xr
        if f(xr.c()) > f(g.c()):
            
            #сжатие
            xc = mid + beta * (w - mid)
            if f(xc.c()) < f(w.c()):
                w = xc

        # новые точки
        v1 = w
        v2 = g
        v3 = b
    return b

print("Результаты (метод деформируемого многогранника)")
rf1 = mnogogr(f1)
rf2 = mnogogr(f2)

print("Min f1 %s"%(rf1))
print("Min f2 %s"%(rf2))