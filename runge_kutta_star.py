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

m = 1.989 * 10**30              # mass of the sun in kg
M = 10**6 * m                   # mass of BH

r_star = 696340*10**3

'''
Initial conditions
'''

# we're gonna use kepler = 2 potential
# as soon as r<r_isco the orbits become unstable

x0 =    0
y0 =    8*10**9     # close to the isco radius
# v_x0 =  10**8
v_y0 =  0

# orbit velo
# v_x0 = np.sqrt(G*M/y0)

# orbit velo for BH
# it actually works...
v_x0 = np.sqrt(G*M/y0 + 3*G**2*M**2/(y0**2*c**2))

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
    print(x,y)
    x = float(x)
    y = float(y)
    r = np.sqrt((x**2+y**2))
    return r*m*v_x0

L = L(x0, y0, m)


A = G*M
B = L**2 / m**2             # 0.005
C = 3*G*M*L**2 / (m*c**2)   # 0.8


# things for the Black hole

r_ss   = 2*A/c**2
r_isco = 3*r_ss

# stable orbit for photons
r_photon = 1.5 * r_ss

# tidal radius
r_tidal = r_star * (M/m)**(1/3)

##############################################################################################################
#                               Kepler orbit or BH orbit? 
kepler = 2              
# set to True for Newtonian potential, False for BH potential, which does not seem to work at all
# and to 2 for the alternative potential from the paper
##############################################################################################################


# small r, high v are most intersting region for us 

'''
for the kepler potential we see that we still have an orbit for r<r_isco
That makes total sense, there is no isco for keplerian potential
these values give reasonable results for the keplerian potential: 

m = 1.989 * 10**30              # mass of the sun in kg
M = 10**6 * m                   # mass of BH

r_star = 1 

A = G*M


# things for the Black hole

r_ss   = 2*A/c**2
r_isco = 3*r_ss

x0 =    0
y0 =    5*10**9    
v_y0 =  0

# orbit velo
v_x0 = np.sqrt(A/y0)

'''





# escape velocity for newtonian potential
# v0 = np.sqrt(2/r0) 

t = np.linspace(0,100000,10000)

t0 = t[0]

h = len(t)

def r_help(A, B, C):
    '''
    C is the force on the ith particle
    '''

    # print('################################')
    # print((27 * B**2 * C**4 - 4 * A**3 * C**3)**(1/2))
    # print(27 * B**2 * C**4 - 4 * A**3 * C**3)
    # print('################################')

    r_thr =  ((2/3)**(1/3) * A)/(np.sqrt(3) * (27 * B**2 * C**4 - 4 * A**3 * C**3)**(1/2) + 9 * B * C**2)**(1/3) + (np.sqrt(3) * (27 * B**2 * C**4 - 4 * A**3 * C**3)**(1/2) + 9 * B * C**2)**(1/3)/(2**(1/3) * 3**(2/3) * C)

    return r_thr

def r_help_bh(A,B,C,D):
    A = G*M
    B = L**2 / m**2             # 0.005
    C = 3*G*M*L**2 / (m*c**2)   # 0.8

# Try, aber hat nicht funktioniert ... nur nan als number bekommen

# def grav_layer(r_star, m_star, split = 5):
#     '''
#     Calculate the different radii for the different layers
#     assuming constant density
#     '''
#     m = m_star 
#     r1 = r_star / (split**(1/3)) 
#     r1 = r1
#     r2 = (2**(1/3)-1)*r1
#     r3 = (3**(1/3)-2**(1/3))*r1
#     r4 = (4**(1/3)-3**(1/3))*r1
#     r5 = r_star

#     # gravitational pull on i-th layer by star
#     F5 = G*1/5*m*(m-1/5*m)/(r4**2)
#     F4 = G*1/5*m*(m-2/5*m)/(r3**2)
#     F3 = G*1/5*m*(m-3/5*m)/(r2**2)
#     F2 = G*1/5*m*(m-4/5*m)/(r1**2)

   


#     A = G*M*m
#     B = 3*(G*M/c)**2*m

#     print(A, B, F3)
#     r2_thr = r_help(A,B,F2).real
#     r3_thr = r_help(A,B,F3).real
#     r4_thr = r_help(A,B,F4).real
#     r5_thr = r_help(A,B,F5).real
#     # eq to solve: -A/r^2 - B/r^3 = -F_i for r

#     return r2_thr, r3_thr, r4_thr, r5_thr

def grav_layer(r_star, m_star, split = 5):
    '''
    Calculate the different radii for the different layers
    assuming constant density
    '''
    m = m_star 
    r1 = r_star / 5
    r2 = 2*r1
    r3 = 3*r1
    r4 = 4*r1
    r5 = r_star

    # gravitational pull on i-th layer by star
    F5 = G*1/5*m*(m-1/5*m)/((r4)**2)
    F4 = G*1/5*m*(m-2/5*m)/((r3)**2)
    F3 = G*1/5*m*(m-3/5*m)/((r2)**2)
    F2 = G*1/5*m*(m-4/5*m)/((r1)**2)
    F_whole_star = G*m**2/(r_star/2)**2

    A = G*M*m
    B = 3*(G*M/c)**2*m

    print(F2, F3, F4, F5) 

    print(A, B, F3)
    r2_thr = r_help(A,B,F2).real
    r3_thr = r_help(A,B,F3).real
    r4_thr = r_help(A,B,F4).real
    r5_thr = r_help(A,B,F5).real

    r_whole_star = r_help(A,B, F_whole_star).real
    # eq to solve: -A/r^2 - B/r^3 = -F_i for r

    return r2_thr, r3_thr, r4_thr, r5_thr



r2, r3, r4, r5 = grav_layer(r_star, m)

print(r2, r3, r4, r5)
print(r2/r_tidal, r3/r_tidal, r4/r_tidal, r5/r_tidal)

print(r2/r_isco, r3/r_isco, r4/r_isco, r5/r_isco)





'''
What we have to do now:
Look at every timestep -> if F_BH on the ith layer is greater than the constant value of F_i
-> layer turns into n particles
-> these particles are only interacting with the black hole now, not with the star, not with each other 
-  all of these particles have the same initial velocity, but different initial positions
'''


# as proposed in the paper: 
# https://arxiv.org/pdf/2008.04922.pdf

def v_x_(t,x,y,v_x,v_y):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -A*x / r**3 #+ B*x / r**4 - C*x / r**5
    elif kepler == 2:
        return -A*x / r**3 - 3*(G*M/c)**2*x / r**4
    else:
        return -A*x / r**3 + B*x / r**4 - C*x / r**5

def v_y_(t,x,y,v_x,v_y):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -A*y / r**3 #+ B*y / r**4 - C*y / r**5
    elif kepler == 2:
        return -A*y / r**3 - 3*(G*M/c)**2*y / r**4
    else: 
        return -A*y / r**3 + B*y / r**4 - C*y / r**5

# this acts as v = dot(x) and v = dot(y)
def x_(t,x,y,v_x,v_y):
    return v_x

def y_(t,x,y,v_x,v_y):
    return v_y





def my_rk4(x0, v_x0, y0, v_y0, h = 0.2):
    '''
    This function solves the differential equations (the equation of motion)
    for the x- and y-component of a particle. At each step, we evaluate v[i] for x and y which corresponds to their derivatives 
    as well as the derivatives dot(v) for x and y which corresponds to their second order derivatives. 
    Since for the (i+1)-th step, we need the values of x, y, v_x and v_y for the i-th step, for each step we calculate these 4 values
    '''
    x = np.zeros(len(t))
    v_x = np.zeros(len(t))
    y = np.zeros(len(t))
    v_y = np.zeros(len(t))
    x[0] = x0
    v_x[0] = v_x0
    y[0] = y0
    v_y[0] = v_y0
    for i in range(0,len(t)-1):
        k1x =   (x_(    t[i], x[i], y[i], v_x[i], v_y[i]))
        k1v_x = (v_x_(  t[i], x[i], y[i], v_x[i], v_y[i]))
        k1y =   (y_(    t[i], x[i], y[i], v_x[i], v_y[i]))
        k1v_y = (v_y_(  t[i], x[i], y[i], v_x[i], v_y[i]))

        k2x =   (x_(    t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2))
        k2v_x = (v_x_(  t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2))
        k2y =   (y_(    t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2))
        k2v_y = (v_y_(  t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2))

        k3x =   (x_(    t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2))
        k3v_x = (v_x_(  t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2))
        k3y =   (y_(    t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2))
        k3v_y = (v_y_(  t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2))
        
        k4x =   (x_(    t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2))
        k4v_x = (v_x_(  t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2))
        k4y =   (y_(    t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2))
        k4v_y = (v_y_(  t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2))

        x[i+1]  = x[i] + (k1x+2*k2x+2*k3x+k4x)*h/6
        v_x[i+1]  = v_x[i] + (k1v_x+2*k2v_x+2*k3v_x+k4v_x)*h/6
        y[i+1]  = y[i] + (k1y+2*k2y+2*k3y+k4y)*h/6
        v_y[i+1]  = v_y[i] + (k1v_y+2*k2v_y+2*k3v_y+k4v_y)*h/6


        # define black hole radius to be one right now
        if np.sqrt(x[i]**2+y[i]**2) <= 10**(-17):
            x[i] = 0
            y[i] = 0
            x = x[:i]
            y = y[:i]
            v_x = v_x[:i]
            v_y = v_y[:i]
            break
    return x, y, v_x, v_y


# def test():
#     '''
#     function that calculates, when the F_BH on the i-th layer of the star is bigger than the gravitational pull of the star itself
#     '''
#     r = np.sqrt(x**2+y**2)
#     F_bh = -A / r**2 - 3*(G*M/c)**2 / r**3
#     F2, F3, F4, F5 = grav_layer(r_star, m)
#     for i in range(len(t)):
#         if F_bh <= F2:
#             print('Distance for first layer to drop: ', r)
#         if F_bh <= F3:
#             print('Distance for second layer to drop: ', r)
#         if F_bh <= F4:
#             print('Distance for third layer to drop: ', r)
#         if F_bh <= F5:
#             print('Distance for fourth layer to drop: ', r)

# test()


# write a function that creates multiple particles out of a layer

def particles(layer, x_star, y_star, vx_star, vy_star, number_of_particles=5):
    '''
    All particles have the same velocity but different initial positions
    '''
    vx0 = vx_star
    vy0 = vy_star

    # x10 = 
    # x20 = 
    # x30 = 
    # x40 = 
    # x50 = 

    # y10 = 
    # y20 = 
    # y30 = 
    # y40 = 
    # y50 = 


    # for i in range(number_of_particles):

    




