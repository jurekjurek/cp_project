import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants

# rk4 only works for 1st order DE

G = constants.G
c = constants.speed_of_light

m = 1.989 * 10**30              # mass of the sun in kg
M = 10**6 * m                   # mass of BH

L = 100

A = G*M
B = L**2 / m**2
C = 3*G*M*L**2 / (m*c**2)


r0 = 3
v0 = 0
t = np.linspace(0,100)

t0 = t[0]

h = len(t)

# function to be solved: dy/dx = x + y
def f1(t,r,v):
    return v

def f2(t,r,v, k=1):
    return -k/m*r

# du/dx = f(x)
def u(x,y):
    r = np.sqrt(x**2+y**2)
    return -A/r**2 + B/(r**3)+C/(r**4)

# or
# f = lambda x: x+y

# RK-4 method
def rk4(t0,r0,v0,n, f1, f2):
    
    # Calculating step size
    # h = (xn-x0)/n
    
    print('\n--------SOLUTION--------')
    print('-------------------------')    
    # print('x0\ty0\tyn')
    print('-------------------------')
    for i in range(n):
        k1r = h * (f1(t0, r0, v0))
        k1v = h * (f2(t0, r0, v0))
        k2r = h * (f1(t0+ h/2, r0 + h*k1r/2, v0 + h*k1v/2))
        k2v = h * (f2(t0+ h/2, r0 + h*k1r/2, v0 + h*k1v/2))
        k3r = h * (f1(t0+ h/2, r0 + h*k2r/2, v0 + h*k2v/2))
        k3v = h * (f2(t0+ h/2, r0 + h*k2r/2, v0 + h*k2v/2))
        k4r = h * (f1(t0+ h, r0 + h*k3r/2, v0 + h*k3v/2))
        k4v = h * (f2(t0+ h, r0 + h*k3r/2, v0 + h*k3v/2))
        kr  = (k1r+2*k2r+2*k3r+k4r)/6
        kv  = (k1v+2*k2v+2*k3v+k4v)/6
        rn = r0 + kr
        vn = v0 + kv
        # yn = y0 + kr
        # print('%.4f\t%.4f\t%.4f'% (x0,y0,yn) )
        print('-------------------------')
        # y0 = yn
        r0 = rn
        v0 = vn
        # x0 = x0+h
    
    print(type(rn), type(vn))
    print('\nAt x=%.4f, y=%.4f' %(rn, vn))
    # print('\nAt x=%.4f, y=%.4f' %(xn,yn))

# Inputs
# print('Enter initial conditions:')
# x0 = float(input('x0 = '))
# y0 = float(input('y0 = '))
# x0 = 0
# y0 = 1

# print('Enter calculation point: ')
# xn = float(input('xn = '))
# xn = 2

print('Enter number of steps:')
step = int(input('Number of steps = '))

# RK4 method call
# rk4(x0,y0,xn,step, u)
rk4(t0,r0,v0,step, f1, f2)






# funktioniert vielleicht, keine Ahnung ;) 



