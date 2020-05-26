import pandas as pd
import os
os.chdir('D:/_h3110fri3nd/git/kursoCHM/res')

alpha = 3
beta = 1
gamma = 1
y0 = lambda x: 0.8


q = lambda x: 1 + x
p = lambda x: 1 + x**gamma
u = lambda x: x**alpha*(1-x)**beta
F = lambda x: -x**5+17*x**3+3*x**2-6*x


ns = [10, 20, 30, 40]
ms = [2, 3, 4, 5, 6]


# Upper relaxation
dfUR = {n: pd.DataFrame(columns=['eps', 'minW', 'minK', 'max|y-u|' ], index = ms) for n in ns}
for n in ns:
    for m in ms:
        h = 1/n
        eps = h**m
        mink = 10**6
        a = [p(h*i) for i in range(0, n+1)]
        g = [q(h*i) for i in range(0, n+1)]
        f = [F(i*h) for i in range(0, n+1)]
        for w in [round(1+i/100, 2) for i in range(1, 100)]:
            y = [y0(i*h) for i in range(0, n+1)]
            k = 0
            r = [eps for i in range(1, n)]
            while (max(r[1:n])>=eps):
               k+=1
               for i in range(1, n):
                   y[i] = (1-w)*y[i]+w*(f[i]*h*h+a[i]*y[i-1]+a[i+1]*y[i+1])/(a[i]+a[i+1]+g[i]*h*h)
               r = [abs(-a[i]*y[i-1]+(a[i]+a[i+1]+h*h*g[i])*y[i]-a[i+1]*y[i+1] - f[i]*h*h) for i in range(1, n)]    
            print ("n = ", n, "m = ", m, "w = ", w, "k = ", k)
            if k < mink:
                mink = k
                minw = w
                yminw = [y[i] for i in range(0, n+1)]
        print('mink = ', mink)
        dfUR[n].loc[m] = {'eps':eps, 'minW': minw, 'minK': mink, 'max|y-u|': max(abs(yminw[i]-u(i*h)) for i in range(1, n))}
        print('\n')
        
for key in dfUR.keys():
    dfUR[key].to_excel('UR_y0=08''%s''.xlsx' % (key))
        
# Seidel
dfSeidel = {n: pd.DataFrame(columns=['eps', 'k', 'max|y-u|' ], index = ms) for n in ns}
for j in range(len(ns)):
    n = ns[j]
    for m in ms:
        h = 1/n
        eps = h**m
        a = [p(h*i) for i in range(0, n+1)]
        g = [q(h*i) for i in range(0, n+1)]
        f = [F(i*h) for i in range(0, n+1)]
        y = [y0(i*h) for i in range(0, n+1)]
        k = 0
        r = [eps for i in range(1, n)]
        while (max(r[1:n])>=eps):
            k+=1
            for i in range(1, n):
                y[i] = (f[i]*h*h+a[i]*y[i-1]+a[i+1]*y[i+1])/(a[i]+a[i+1]+g[i]*h*h)
            r = [abs(-a[i]*y[i-1]+(a[i]+a[i+1]+h*h*g[i])*y[i]-a[i+1]*y[i+1] - f[i]*h*h) for i in range(1, n)]      
        print ("n = ", n, "m = ", m, "k = ", k)
        dfSeidel[n].loc[m] = {'eps':eps, 'k': k, 'max|y-u|': max(abs(y[i]-u(i*h)) for i in range(1, n))}
for key in dfSeidel.keys():
    dfSeidel[key].to_excel('Seidel_y0=08''%s''.xlsx' % (key))
    
# Jakobi
dfJakobi = {n: pd.DataFrame(columns=['eps', 'k', 'max|y-u|' ], index = ms) for n in ns}
for j in range(len(ns)):
    n = ns[j]
    for m in ms:
        h = 1/n
        eps = h**m
        a = [p(h*i) for i in range(0, n+1)]
        g = [q(h*i) for i in range(0, n+1)]
        f = [F(i*h) for i in range(0, n+1)]
        y = [y0(i*h) for i in range(0, n+1)]
        k = 0
        r = [eps for i in range(1, n)]
        while (max(r[1:n])>=eps):
            k+=1
            yprev = [y[i] for i in range(0, n+1)]
            for i in range(1, n):
                y[i] = (f[i]*h*h+a[i]*yprev[i-1]+a[i+1]*yprev[i+1])/(a[i]+a[i+1]+g[i]*h*h)
            r = [abs(-a[i]*y[i-1]+(a[i]+a[i+1]+h*h*g[i])*y[i]-a[i+1]*y[i+1] - f[i]*h*h) for i in range(1, n)]      
        print ("n = ", n, "m = ", m, "k = ", k)
        dfJakobi[n].loc[m] = {'eps':eps, 'k': k, 'max|y-u|': max(abs(y[i]-u(i*h)) for i in range(1,n))}
for key in dfJakobi.keys():
    dfJakobi[key].to_excel('Jakobis_y0=08''%s''.xlsx' % (key))
    
pd.Series([abs(dfUR[key].loc[6]['max|y-u|']) for  key in dfUR.keys()], index = dfUR.keys()).plot()
pd.Series([abs(dfSeidel[key].loc[6]['max|y-u|']) for  key in dfSeidel.keys()], index = dfSeidel.keys()).plot()
pd.Series([abs(dfJakobi[key].loc[6]['max|y-u|']) for  key in dfJakobi.keys()], index = dfJakobi.keys()).plot()