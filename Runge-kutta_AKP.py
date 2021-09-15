import numpy as np
import matplotlib.pyplot as plt
#%%
G = 6.7*10**(-11)
M = 2*10**(30)
f =  lambda m,x,y,vx,vy,t: (-m*G*M*x/(m*x**2+y**2)**(3/2))
g =  lambda m,x,y,vx,vy,t: (-y*G*M/(m*x**2+y**2)**(3/2))

def rk4sol(f,g,mu,a,b,x0,y0,vx0,vy0,t0,h):
    x = np.arange(a, b+h, h)
    y = np.arange(a, b+h, h)
    vx = np.arange(a, b+h, h)
    vy = np.arange(a, b+h, h)
    t = np.arange(a, b+h, h)
    n = len(x)
    x = np.zeros(n)
    y = np.zeros(n)
    vx = np.zeros(n)
    vy = np.zeros(n)
    t = np.zeros(n)
    x[0] = x0
    y[0] = y0
    vx[0] = vx0
    vy[0] = vy0
    t[0] = t0
    m = 'Órbita de Marte con mu = ' + str(mu)
    for i in range(0,n-1):
        k1 = h*vx[i]
        l1 = h*f(mu,x[i],y[i],vx[i],vy[i],t[i])
        q1 = h*vy[i]
        m1 = h*g(mu,x[i],y[i],vx[i],vy[i],t[i])
        
        k2 = h*(vx[i] + l1/2)
        l2 = h*f(mu,x[i] + k1/2,y[i] + q1/2,vx[i] + l1/2,vy[i] + m1/2,t[i] + h/2)
        q2 = h*(vy[i] + m1/2)
        m2 = h*g(mu,x[i] + k1/2,y[i] + q1/2,vx[i] + l1/2,vy[i] + m1/2,t[i] + h/2)
        
        k3 = h*(vx[i] + l2/2)
        l3 = h*f(mu,x[i] + k2/2,y[i] + q2/2,vx[i] + l2/2,vy[i] + m2/2,t[i] + h/2)
        q3 = h*(vy[i] + m2/2)
        m3 = h*g(mu,x[i] + k2/2,y[i] + q2/2,vx[i] + l2/2,vy[i] + m2/2,t[i] + h/2)
        
        k4 = h*(vx[i] + l3)
        l4 = h*f(mu,x[i] + k3,y[i] + q3,vx[i] + l3,vy[i] + m3,t[i] + h)
        q4 = h*(vy[i] + m3)
        m4 = h*g(mu,x[i] + k3,y[i] + q3,vx[i] + l3,vy[i] + m3,t[i] + h)
    
        x[i+1] = x[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        y[i+1] = y[i] + (1/6)*(q1 + 2*q2 + 2*q3 + q4)
        vx[i+1] = vx[i] + (1/6)*(l1 + 2*l2 + 2*l3 + l4)
        vy[i+1] = vy[i] + (1/6)*(m1 + 2*m2 + 2*m3 + m4)
        
        
    plt.scatter(x, y, s=0.05)
    plt.scatter(x[0],y[0],color='r',label='punto inicial')
    plt.scatter(x[-1],y[-1],color='g',label='punto final')
    plt.scatter(0,0,color='k')
    plt.legend(loc='lower right')
    plt.title(m)
    plt.show()
