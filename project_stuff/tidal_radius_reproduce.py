import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants
from sympy import *
# import oct2py

# rk4 only works for 1st order DE

G = constants.G
c = constants.speed_of_light

m = 1.989 * 10**30              # mass of the sun in kg
M = 10**6 * m                   # mass of BH

r_star = 696340*10**3

'''
Initial conditions
'''

# tidal radius
r_tidal = r_star * (M/m)**(1/3)

# x0 =    0
# y0 =    1.2*r_tidal #8*10**9     # close to the isco radius

x0 =    -r_tidal*1.2 
y0 =    r_tidal*1.5

# v_x0 =  10**8
v_y0 =  0

# orbit velo for kepler potential
# v_x0 = np.sqrt(G*M/y0)

# orbit velo for BH
# it actually works...
v_x0 = np.sqrt(G*M/y0 + 3*G**2*M**2/(y0**2*c**2))*0.5
v_y0 = v_x0 *0.4
# escape velo, A = G*M
# we see that for 70% e.g. the particle falls back into the Mass
# v_x0 = 0
# v_y0 = np.sqrt(A*2/y0)*0.7

print('vy_0 is ', v_x0/c*100, '% of the speed of light.')

def L(x,y, m):
    '''
    function that calculates the angular momentum for given x, y
    L is constant over time
    Since at point t = 0, we have only a velocity component in x-direction, so we calculate the L(0) = L(t) to be:
    '''
    x = float(x)
    y = float(y)
    r = np.sqrt((x**2+y**2))
    return r*m*v_x0

L = L(x0, y0, m)


A = G*M
B = L**2 / m**2             # 0.005
C = 3*G*M*L**2 / (m**2*c**2)   # 0.8


# things for the Black hole

r_ss   = 2*A/c**2
r_isco = 3*r_ss

# stable orbit for photons
r_photon = 1.5 * r_ss

# tidal radius
r_tidal = r_star * (M/m)**(1/3)





'''
Here this whole grav_layer stuff starts
'''

def r_help(A, B, C, r_i):
    '''
    C is the force on the ith particle
    Solves -A/r^2 - B/r^3 = -C for r
    '''

    r = symbols('r')
    return solve(-A/((r-r_i)**2)-B/((r-r_i)**3)+C, r)

    r_thr =  ((2/3)**(1/3) * A)/(np.sqrt(3) * (27 * B**2 * C**4 - 4 * A**3 * C**3)**(1/2) + 9 * B * C**2)**(1/3) + (np.sqrt(3) * (27 * B**2 * C**4 - 4 * A**3 * C**3)**(1/2) + 9 * B * C**2)**(1/3)/(2**(1/3) * 3**(2/3) * C)

    return r_thr

def r_help_bh(D, r_i):
    A = G*M*m/5
    B = 1/10*L**2 / (m**2) * m/5             # 0.005
    C = 1/10*3*G*M*L**2 / (m**2*c**2) * m/5   # 0.8

    r = symbols('r')
    return solve(-A / (r-r_i)**2 + B / (r-r_i)**3 - C / (r-r_i)**4 + D, r)



# this is the function that kinda works, not assuming constant density!!

def grav_layer(r_star, m_star, split = 5):
    '''
    Calculate the different radii for the different layers
    assuming constant density
    '''
    m = m_star 
    # r1 = r_star / (split**(1/3)) 
    # r1 = r1
    # r2 = (2**(1/3)-1)*r1
    # r3 = (3**(1/3)-2**(1/3))*r1
    # r4 = (4**(1/3)-3**(1/3))*r1
    # r5 = r_star

    # r1 = r_star / (split**(1/3)) 
    # r1 = r1
    # r2 = (2**(1/3)-1)*r1 #+ r1
    # r3 = (3**(1/3)-2**(1/3))*r1 #+ r1 + r2 
    # r4 = (4**(1/3)-3**(1/3))*r1 #+ r1 + r2 +r3
    # r5 = (5**(1/3)-4**(1/3))*r1 #+ r1 + r2 +r3 + r4  # this is the rad of the star
    
    r1 = r_star / 5
    r2 = r_star / 5 
    r3 = r_star / 5 
    r4 = r_star / 5 
    r5 = r_star / 5 
    
    
    print('should reproduce whole star radius...: ', (r1+r2+r3+r4+r5)/100000, r_star/100000)
    
    r1_tot = r1
    r2_tot = r2 + r1
    r3_tot = r3+r2+r1
    r4_tot = r4+r3+r2+r1

    # F5 = G*m**2/(((r4_tot+r_star)/2)**2)
    # F4 = G*(4/5*m**2)/(((r3_tot+r4_tot)/2)**2)

    # F3 = G*(3/5*m**2)/(((r3_tot+r2_tot)/2)**2)
    # F2 = G*(2/5*m**2)/(((r2_tot+r1_tot)/2)**2)
    # F1 = G*(1/5*m**2)/(((r1_tot+0)/2)**2)


    # F_test = G

    F5 = G*2/10*m*(9/10*m)/(((r4_tot+r_star)/2)**2)
    F4 = G*2/10*m*(7/10*m)/(((r3_tot+r4_tot)/2)**2)

    F3 = G*2/10*m*(5/10*m)/(((r3_tot+r2_tot)/2)**2)
    F2 = G*2/10*m*(3/10*m)/(((r2_tot+r1_tot)/2)**2)
    F1 = G*2/10*m*(1/10*m)/(((r1_tot+0)/2)**2)

    print('radii 2,3,4,5: ', r1/10000,r2/10000,r3/10000,r4/10000)

    print('Forces ascending: ',F1, F2, F3, F4, F5)


    A5 = G*M*m*1/5
    A4 = G*M*m*1/5
    A3 = G*M*m*1/5
    A2 = G*M*m*1/5

    B5 = 3*(G*M/c)**2*m*1/5
    B4 = 3*(G*M/c)**2*m*1/5
    B3 = 3*(G*M/c)**2*m*1/5
    B2 = 3*(G*M/c)**2*m*1/5

    # r1_thr = np.array([complex(item) for item in r_help(A2,B2,F1, r1_tot)])[2].real
    # r2_thr = np.array([complex(item) for item in r_help(A2,B2,F2, r2_tot)])[2].real
    # r3_thr = np.array([complex(item) for item in r_help(A3,B3,F3, r3_tot)])[2].real
    # r4_thr = np.array([complex(item) for item in r_help(A4,B4,F4, r4_tot)])[2].real
    # r5_thr = np.array([complex(item) for item in r_help(A5,B5,F5, r_star)])[2].real 

    r1_thr = np.array([complex(item) for item in r_help_bh(F1, r1_tot/2)])[2].real
    r2_thr = np.array([complex(item) for item in r_help_bh(F2, (r2_tot+r1_tot)/2)])[2].real
    r3_thr = np.array([complex(item) for item in r_help_bh(F3, (r3_tot+r2_tot)/2)])[2].real
    r4_thr = np.array([complex(item) for item in r_help_bh(F4, (r4_tot+r3_tot)/2)])[2].real
    r5_thr = np.array([complex(item) for item in r_help_bh(F5, (r_star+r4_tot)/2)])[2].real
 
    print(r1_thr/100000, r2_thr/100000, r3_thr/100000, r4_thr/100000, r5_thr/100000)

    # eq to solve: -A/r^2 - B/r^3 = -F_i for r
    # return F2,F3,F4,F5
    return r1_thr-r1_tot/2, r2_thr-(r2_tot+r1_tot)/2, r3_thr-(r3_tot+r2_tot)/2, r4_thr-(r4_tot+r3_tot)/2, r5_thr-(r_star+r4_tot)/2


'''
Don't assume constant density probably
'''

r1, r2,r3,r4,r5 = grav_layer(r_star, m)

print(r2, r3, r4, r5)
print(r1/r_tidal, r2/r_tidal, r3/r_tidal, r4/r_tidal, r5/r_tidal)





