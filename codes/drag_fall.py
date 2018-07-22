import numpy as np
import matplotlib.pyplot as plt

m = .5
g = 9.8
k = 0.5
h = .05  # 刻み幅
t_max = 10


def f(x, t):
    return g - k/m * x


t = np.arange(0, t_max, h)

v_1 = m * g / k * (1 - np.exp(-k / m * t))  # 解析解

v_2 = np.arange(0, t_max, h) #Euler法
v_3 = np.arange(0, t_max, h) #Runge-Kutta法

v_2[0] = 0.0


for i in range(v_2.size - 1):
    
    #Euler
    v_2[i + 1] = v_2[i] + f(v_2[i], t[i])*h
    
    #Runge-Kutta
    k_1 = h*f(v_3[i],t)
    k_2 = h*f(v_3[i]+k_1/2,t[i]+h/2)
    k_3 = h*f(v_3[i]+k_2/2,t[i]+h/2)
    k_4 = h*f(v_3[i]+k_3,t[i]+h)
    v_3[i+1]=v_3[i]+(k_1+2*k_2+2*k_3+k_4)/6

plt.title("numerical calculation\t($h=$"+str(h)+")")
plt.xlabel("time [s]")
plt.ylabel("velocity [$m\cdot s^{-1}$]")
plt.plot(t, v_1,label="Exact solution")
plt.plot(t, v_2,label="Euler")
plt.plot(t, v_3,label="Lunge-Kutta")
plt.legend()
plt.show()
