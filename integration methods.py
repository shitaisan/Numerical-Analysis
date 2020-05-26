import pandas as pd
import numpy as np

a = 0.4
b = 4.0
eps = 0.00001


def Ci(x):
    a = -x*x/4
    sum = 0
    n = 1
    while abs(a) > eps:
        sum += a
        a *= (-n*x*x/((4*n+2)*(n+1)*(n+1)))
        n += 1
    return sum 

def Ln(x, nodes):
    res = 0
    for xi in nodes:
        li = 1
        for xj in nodes:
            if xj != xi:
                li *= (x-xj)/(xi-xj)
        res += li*Ci(xi)
    return res


def zadanie1():
    x = [round(a + (b-a) / 20 * i, 2) for i in range(21)]
    f = [Ci(each) for each in x]
    return pd.DataFrame({'f': f}, index=x)

def zadanie2_1():
    nodes = [round(a + (b-a) / 10 * i, 2) for i in range(11)]
    f = zadanie1()
    errors = [abs(f[x]-Ln(x, nodes)) for x in f.index]
    pd.Series(errors, f.index).plot()
    
def zadanie2_2(isChebs):
    ns = [5, 8, 14, 18]
    if (isChebs): 
        n_nodes = [[round((a+b+(b-a)*np.cos((2*k+1)*np.pi/(2*n)))/2, 2) for k in range(n)] for n in ns]
    else: 
        n_nodes = [[round(a + (b-a) / (n-1) * i, 2) for i in range(n)] for n in ns]
    f = zadanie1()
    return pd.DataFrame([pd.Series([abs(f.loc[x]['f']-Ln(x, nodes)) for x in f.index], index = f.index) for nodes in n_nodes], index=ns)

def zadanie2_3(isChebs):
    ns = [i for i in range(100, 500, 10)]
    if (isChebs): 
        n_nodes = [[(a+b)/2. + (b-a)/2.*np.cos(np.pi*(2*k+1)/(2.*n))  for k in range(n)] for n in ns]
    else: 
        n_nodes = [[a + (b-a) / (n-1) * i for i in range(n)] for n in ns]
    f = zadanie1()
    max_errors = pd.Series([max([abs(f.loc[x]['f']-Ln(x, nodes)) for x in f.index]) for nodes in n_nodes], [len(nodes) for nodes in n_nodes])
    max_errors.plot()


def rightrects(f, c, d, n):
    h = (d-c)/n
    z = [c+i*h for i in range(1, n+1)]    
    return h*sum(f(z[i+1]) for i in range(n-1))

def centralrects(f, c, d, n):
    h = (d-c)/n
    z = [c+i*h for i in range(1, n+1)]    
    return h*sum(f((z[i+1]+z[i])/2) for i in range(n-1))
    
def simpson(f, c, d, n):
    h = (d-c)/n
    z = [c+i*h for i in range(1, n+1)]    
    return h/6*sum(f(z[i])+f(z[i+1])+4*f((z[i]+z[i+1])/2) for i in range(n-1))
    
   
def gauss(f, c, d, n):
    h = (d-c)/n
    z = [c+i*h for i in range(1, n+1)]    
    return h/2*sum(f(z[i-1]+h/2*(1-1/np.sqrt(3)))+f(z[i-1]+h/2*(1+1/np.sqrt(3))) for i in range(1, n))
    


def methodtable(method, eps):
    X = [round(a + (b-a) / 10 * i, 2) for i in range(11)]
    df = pd.DataFrame(columns=['N', 'S'], index = X)
    f = lambda x: (np.cos(x)-1)/x
    for x in X:
        n = 1
        cur, prev = eps, 0
        while (abs(cur-prev)>=eps):
            n *= 2
            prev = cur
            if method == 'rightrects':
                cur = rightrects(f, 0, x, n)
            elif method == 'centralrects':
                cur = centralrects(f, 0, x, n)
            elif method == 'simpson':
                cur = simpson(f, 0, x, n)
            elif method == 'gauss':
                cur = gauss(f, 0, x, n)
        df.loc[x] = {'N': n, 'S': cur}
    return df

def Newton():
    n = 10
    X = [round(a + (b-a) / n * i, 2) for i in range(n+1)]
    f = lambda x: (np.cos(x)-1)/x
    F = lambda k: Ci(X[0])+k*(Ci(X[n])-Ci(X[0]))/n
    df = pd.DataFrame(columns = ['z', 'F'])
    for i in range(n+1):
        z = X[i]
        while np.abs(Ci(z)-F(i))>eps:
            z = z - (Ci(z)-F(i))/f(z)
        df.loc[i] = {'z': Ci(z), 'F': z}
    return df

def chord():
    n = 10
    X = [round(a + (b-a) / n * i, 2) for i in range(n+1)]
    f = lambda x: (np.cos(x)-1)/x
    F = lambda k: Ci(X[0])+k*(Ci(X[n])-Ci(X[0]))/n
    df = pd.DataFrame(columns = ['z','F'])
    for i in range(n+1):
        z0= X[i]
        zk = z0 - (Ci(z0)-F(i))/f(z0)
        while np.abs(Ci(zk)-F(i))>eps:
            zk = (z0*(Ci(zk)-F(i))-zk*(Ci(z0)-F(i)))/(Ci(zk)-Ci(z0))
        df.loc[i] = {'z': Ci(zk), 'F': zk}
    return df


def secant():
    n = 10
    X = [round(a + (b-a) / n * i, 2) for i in range(n+1)]
    F = lambda k: Ci(X[0])+k*(Ci(X[n])-Ci(X[0]))/n
    df = pd.DataFrame(columns = ['z','F'])
    for i in range(n+1):
        z1 = X[i]
        z2 = X[(i+1)%(n+1)]
        z3 = z2 - (z2-z1)*(Ci(z2)-F(i))/(Ci(z2)-Ci(z1))
        while np.abs(Ci(z3)-F(i))>eps:
            z1 = z2
            z2 = z3
            z3 = z2 - (z2-z1)*(Ci(z2)-F(i))/(Ci(z2)-Ci(z1))
        df.loc[i] = {'z': Ci(z3), 'F': z3}
    return df



def zadanie3(method):
    return {eps: methodtable(method, eps) for eps in [0.001, 0.0001, 0.00001, 0.000001]}

def zadanie4(method):
    if method =='Newton':
        s = Newton()
    elif method == 'chord':
        s = chord()
    elif method == 'secant':
        s = secant()
    pd.Series([Ci(x) for x in s['F']], index = s['F'] ).plot()
    pd.Series([s.loc[i]['F'] for i in s.index], index = [s.loc[i]['z'] for i in s.index]).plot().set_xlabel(None)
    return s

