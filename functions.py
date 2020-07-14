import math

cBoundVDP = [[-3., 3.], [3., -3.]]
cBoundBogdanov = [[-4., 4.], [4., -4.]]
cBoundZaslavskii = [[-1., 1.], [1., -1.]]
cBoundDuffing = [[-2., 2.], [2., -2.]]
cBoundJulia = [[-3., 3.], [3., -3.]]
cBoundHenon = [[-5., 10.], [5., -10.]]
cBoundIkeda = [[-4.0, 4.0], [4.0, -4.0]]
cBoundDLM = [[0., 1.], [1., 0.]]
cBoundSympM = [[-5., 5.], [5., -5.]]
cBoundRM = [[-20., 20.], [20., -20.]]
cBoundTNK = [[-2., 2.], [2., -2.]]
cBoundPJA = [[-3., 3.], [3., -3.]]
cBoundGBM = [[-10., 10.], [10., -10.]]
cBoundRK = [[0.0, 2.0], [2.0, 0.0]]
cBoundHM = [[-2., 2.], [2., -2.]]
cBoundKR = [[0.0, 1.0], [1.0, 0.0]]
cBoundZS = [[2.0, 5.0], [5.0, 2.0]]
cBoundPH = [[0.0, 5.0], [5.0, 0.0]]

def henonMap(coord):    
    x = coord[0]
    y = coord[1]

    a = 1.4
    b = 0.3
    
    updated_coord = []
    updated_coord.append(1 - a * coord[0] ** 2 + b * coord[1])
    updated_coord.append(coord[0])

    return updated_coord

def rKN(x, fx, n, hs):
    k1 = []
    k2 = []
    k3 = []
    k4 = []
    xk = []
    for i in range(n):
        k1.append(fx[i](x)*hs)
    for i in range(n):
        xk.append(x[i] + k1[i]*0.5)
    for i in range(n):
        k2.append(fx[i](xk)*hs)
    for i in range(n):
        xk[i] = x[i] + k2[i]*0.5
    for i in range(n):
        k3.append(fx[i](xk)*hs)
    for i in range(n):
        xk[i] = x[i] + k3[i]
    for i in range(n):
        k4.append(fx[i](xk)*hs)
    for i in range(n):
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6
    return x

def fa1(x):
    return 0.9*(1 - x[1]*x[1])*x[0] - x[1]

def fb1(x):
    return x[0]

def fc1(x):
    return 0.5

def sympMapF1(x):
    return (x[0] + x[1]) % (2 * math.pi)

def sympMapF2(x):
    nu = 0.5
    return ((x[1] - nu * math.sin(x[0]+x[1]))) % (2 * math.pi)

def VDP1(coord):
    f = [fa1, fb1, fc1]
    x = [coord[0], coord[1], 0]
    hs = 0.01
    for i in range(50):
        x = rKN(x, f, 2, hs)
    updated_coord = []
    updated_coord.append(x[0])
    updated_coord.append(x[1])
    return updated_coord

def sympMap(coord):
    f = [sympMapF1, sympMapF2]
    x = [coord[0], coord[1]]
    hs = 0.01
    for i in range(100):
        x = rKN(x, f, 2, hs)
    updated_coord = []
    updated_coord.append(x[0])
    updated_coord.append(x[1])
    return updated_coord

def rsf1(x):
    return -x[1]
def rsf2(x):
    a = 0.2
    return x[0]+a*x[1]

def rosslerMap(coord):
    f = [rsf1, rsf2]
    x = [coord[0], coord[1]]
    hs = 0.01
    for i in range(20):
        x = rKN(x, f, 2, hs)
    updated_coord = []
    updated_coord.append(x[0])
    updated_coord.append(x[1])
    return updated_coord

def duff1(x):
    return x[1]
def duff2(x):
    delta = 0.2
    beta = -1.
    alpha = 1.
    return -delta*x[1] - beta*x[0] - alpha*x[0]**3

def duffingMap(coord):
    f = [duff1, duff2]
    x = [coord[0], coord[1]]
    hs = 0.02
    for i in range(50):
        x = rKN(x, f, 2, hs)
    updated_coord = []
    updated_coord.append(x[0])
    updated_coord.append(x[1])
    return updated_coord

####################################################################

def bogdanovMap(coord):
    eps = 0.0025
    k = 1.44
    mu = -0.1
    updated_coord = [0, 0]
    x_n = coord[0]
    y_n = coord[1]

    updated_coord[1] = (y_n+eps*y_n+k*x_n*(x_n-1)+mu*x_n*y_n)
    updated_coord[0] = (x_n+updated_coord[1])
    return updated_coord

def zaslavskiiMap(coord):
    x_n = coord[0]
    y_n = coord[1]
    eps = 5
    nu = 0.2
    r = 2
    mu = 1 - (math.e**(-r))/r

    updated_coord = [0, 0]
    updated_coord[0] = (x_n + nu*(1+mu*y_n)+eps*nu*mu*math.cos(2*math.pi*x_n))%1
    updated_coord[1] = (math.e**(-r)*(y_n+eps*math.cos(2*math.pi*x_n)))
    return updated_coord

def tau(X, Y):
    C1 = 0.4
    C3 = 6.
    return C1 - (C3 / (1 + X ** 2 + Y ** 2))

def Ricker(coord):
    #PARAMETRES
    a = 1.8
    r = 2.541
    K = 5
    m = 0.5


    newCoord = []
    newCoord.append(coord[0] * math.exp(r * (1.0 - coord[0] / K) - a * math.pow(coord[1], (1.0 - m))))
    newCoord.append(coord[0] * (1.0 - math.exp( (-a) * math.pow(coord[1], (1.0 - m)))))
    return newCoord


def ikedaMap(coord):
    updated_coord = []
    d = 2.
    C1 = 0.4
    C2 = 0.9
    C3 = 6.
    x = coord[0]
    y = coord[1]
    updated_coord.append(d - C2 * (x * math.cos(tau(x, y)) - y * math.sin(tau(x, y))))
    updated_coord.append(C2 * (x * math.sin(tau(x, y)) + y * math.cos(tau(x, y))))
    return updated_coord


def julia(coord):
    updated_coord = []
    re = -0.12
    im = 0.74
    x = coord[0]
    y = coord[1]
    updated_coord.append(x ** 2 - y ** 2 + re)
    updated_coord.append(2 * x * y + im)

    return updated_coord

def tinkerbellMap(coord):
    updated_coord = []
    a = 0.9
    b = -0.6013
    c = 2.0
    d = 0.3
    x = coord[0]
    y = coord[1]

    updated_coord.append(x**2-y**2+a*x+b*y)
    updated_coord.append(2*x*y+c*x+d*y)
    return updated_coord

def pjattractor(coord):
    updated_coord = []
    a = -2.7
    b = -0.09
    c = -0.86
    d = -2.2
    x = coord[0]
    y = coord[1]

    updated_coord.append(math.sin(a*y)-math.cos(b*x))
    updated_coord.append(math.sin(c*x)-math.cos(d*y))
    return updated_coord

def gbmMap(coord):
    updated_coord = []
    x = coord[0]
    y = coord[1]

    updated_coord.append(1-y+math.fabs(x))
    updated_coord.append(x)
    return updated_coord

def doubleLogMap(coord):
    updated_coord = []
    x = coord[0]
    y = coord[1]
    a = 0.38
    b = 4.11
    updated_coord.append((1 - a) * x + a * b * y * (1 - y))
    updated_coord.append((1 - a) * y + a * b * x * (1 - x))

    return updated_coord

def parasite_host(coord):
    x = coord[0]
    y = coord[1]

    lam = 2.0
    beta = 3.0
    c = 0.2
    p = 2.0 
    b = 1.0
    k = 2.0

    updated_coord = []
    updated_coord.append(lam * x * (k * x)/(b + x**p) * math.exp(-c * y))
    #updated_coord.append(beta * x * (1 - math.exp(-c * y)))
    updated_coord.append(beta * x * (k * x)/(b + x**p) * (1 - math.exp(-c * y)))

    return updated_coord


def template(coord):
    x = coord[0]
    y = coord[1]

    updated_coord = []
    updated_coord.append()
    updated_coord.append()

    return updated_coord

def zaslavskiiMap(coord):
    x_n = coord[0]
    y_n = coord[1]
    eps = 5
    nu = 0.2
    r = 2.0
    mu = 1 - (math.e**(-r))/r

    updated_coord = [0, 0]
    updated_coord[0] = (x_n + nu*(1+mu*y_n)+eps*nu*mu*math.cos(2*math.pi*x_n))%1
    updated_coord[1] = (math.e**(-r)*(y_n+eps*math.cos(2*math.pi*x_n)))
    return updated_coord

def zeraoulia_sprott(coord):
    x = coord[0]
    y = coord[1]
    a = 1
    b = 1

    updated_coord = []
    updated_coord.append(-a * x * (1 - y**2)**(-1))
    updated_coord.append(x + b * y)
    return updated_coord


def f(coord):
    return henonMap(coord)