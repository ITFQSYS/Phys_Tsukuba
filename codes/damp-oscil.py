import numpy as np
import matplotlib.pyplot as plt

k=25
c=1
m=1
h=.05
t_max=h+5


def f(X):
    A=np.array([[0,1],[-k/m,-2*c/m]])
    return A.dot(X)


t=np.arange(0,t_max,h)

x=np.arange(0,t_max,h)
x=0*x
v=np.arange(0,t_max,h)
v=0*v

x_e=np.arange(0,t_max,h)#厳密解

v[0]=2
for i in range(t.size-1):
    X=np.array([x[i],v[i]])
    X.transpose()
    
    #Runge-Kutta
    k_1 = h*f(X)
    k_2 = h*f(X+k_1/2)
    k_3 = h*f(X+k_2/2)
    k_4 = h*f(X+k_3)
    X+=(k_1+2*k_2+2*k_3+k_4)/6

    x[i+1]=X[0]
    v[i+1]=X[1]

    x_e[i]=1/np.sqrt(6)*np.exp(-t[i])*(np.sin(2*np.sqrt(6)*t[i]))

i=t.size - 1
x_e[i]=1/np.sqrt(6)*np.exp(-t[i])*(np.sin(2*np.sqrt(6)*t[i]))

plt.title("h="+str(h))
plt.xlabel("time [s]")
plt.ylabel("x [m]")
plt.plot(t,x,label="Runge-Kutta")
plt.plot(t,x_e,label="Exact")
plt.legend()
plt.show()