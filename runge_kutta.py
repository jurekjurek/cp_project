import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants
# import oct2py

# rk4 only works for 1st order DE

G = constants.G
c = constants.speed_of_light

m = 1 #1.989 * 10**30              # mass of the sun in kg
M = 2 #10**6 * m                   # mass of BH

L = 100

A = 1 # G*M
B = 0.005 #L**2 / m**2
C = 0.8 #3*G*M*L**2 / (m*c**2)


r0 = 1
v0 = 1
t = np.linspace(0,1000,1000)

t0 = t[0]

h = len(t)

# function to be solved: dy/dx = x + y
def f1(t,r,v):
    return v

def f2(t,r,v, k=1, m=1):
    return (-k/m)*r

# du/dx = f(x)
def v_(t,r,v):
    # r = np.sqrt(x**2+y**2)
    return -A/r**2 - B/(r**3)-C/(r**4)

def r_(t,r,v):
    return v

'''
Zweiter Versuch
'''

r = np.zeros(len(t))
v = np.zeros(len(t))

r[0] = r0
v[0] = v0


# def rk4(h = 0.05):
#     for i in range(0,len(t)-1):
#         k1r = (f1(t[i], r[i], v[i]))
#         k1v = (f2(t[i], r[i], v[i]))
#         k2r = (f1(t[i]+ h/2, r[i] + h*k1r/2, v[i] + h*k1v/2))
#         k2v = (f2(t[i]+ h/2, r[i] + h*k1r/2, v[i] + h*k1v/2))
#         k3r = (f1(t[i]+ h/2, r[i] + h*k2r/2, v[i] + h*k2v/2))
#         k3v = (f2(t[i]+ h/2, r[i] + h*k2r/2, v[i] + h*k2v/2))
#         k4r = (f1(t[i]+ h, r[i] + h*k3r/2, v[i] + h*k3v/2))
#         k4v = (f2(t[i]+ h, r[i] + h*k3r/2, v[i] + h*k3v/2))
#         print(k1r, k2r, k3r, k4r)
#         r[i+1]  = r[i] + (k1r+2*k2r+2*k3r+k4r)*h/6
#         v[i+1]  = v[i] + (k1v+2*k2v+2*k3v+k4v)*h/6
#         # rn = r0 + kr
#         # vn = v0 + kv
#     return v, r

def my_rk4(h = 0.2):
    r = np.zeros(len(t))
    v = np.zeros(len(t))
    r[0] = r0
    v[0] = v0
    for i in range(0,len(t)-1):
        k1r = (r_(t[i], r[i], v[i]))
        k1v = (v_(t[i], r[i], v[i]))
        k2r = (r_(t[i]+ h/2, r[i] + h*k1r/2, v[i] + h*k1v/2))
        k2v = (v_(t[i]+ h/2, r[i] + h*k1r/2, v[i] + h*k1v/2))
        k3r = (r_(t[i]+ h/2, r[i] + h*k2r/2, v[i] + h*k2v/2))
        k3v = (v_(t[i]+ h/2, r[i] + h*k2r/2, v[i] + h*k2v/2))
        k4r = (r_(t[i]+ h, r[i] + h*k3r/2, v[i] + h*k3v/2))
        k4v = (v_(t[i]+ h, r[i] + h*k3r/2, v[i] + h*k3v/2))
        print(k1r, k2r, k3r, k4r)
        r[i+1]  = r[i] + (k1r+2*k2r+2*k3r+k4r)*h/6
        v[i+1]  = v[i] + (k1v+2*k2v+2*k3v+k4v)*h/6
        if r[i] <= 0:
            r[i] = 0
            r = r[:i]
            v = v[:i]
            break
    return v, r


# v_res, r_res = rk4()
# print(np.shape(v_res), np.shape(r_res))

my_v, my_r = my_rk4()
print(len(my_v), len(my_r))

plt.figure()
plt.plot(t[:len(my_r)], my_r)
plt.xlabel('Time')
plt.ylabel('Solution of r(t)')
plt.show()

plt.figure()
plt.plot(t[:len(my_r)], my_v)
plt.xlabel('Time')
plt.ylabel('Solution of v(t)')
plt.show()






