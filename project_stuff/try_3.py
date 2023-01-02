import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants

G = constants.G
c = constants.speed_of_light

m = 1.989 * 10**30              # mass of the sun in kg
M = 10**6 * m                   # mass of BH


def eom_r_newton(M, m, r, L):
    return -G*M/r**2 + L**2/(m**2*r**3)

def eom_r_bh(M, m, r, L):
    return -G*M/r**2 + L**2/(m**2*r**3)+3*G*M*L**2/(m*c**2*r**4)

def eom_phi(M, m, r, L):
    mu = m
    return L/(mu*r**2)


