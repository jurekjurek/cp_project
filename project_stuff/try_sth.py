import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants
from sympy import *


# x_save, y_save, vx_save, vy_save = []

# A, B, C, r = symbols('A,B,C,r')
# a = solve(-A/r**2 -B/r**3 - C , r)

# print(a, type(a))
a = np.zeros((5,8))

for i in range(5):
    a[i,:] = np.ones(8)
    print(a)