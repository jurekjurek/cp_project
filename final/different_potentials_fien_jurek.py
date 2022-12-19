import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants

'''
This file compares the orbits for different potential. 
Given initial conditions x, y, vx and vy and the masses of the star and the black hole, the file creates 
an animation of said star orbiting the black hole. 
One can, varying the factor kepler between True and False, explore the effects of the different potentials on the 
orbiting behaviour of the star. 
kepler = False corresponds to a purely newtonian potential, 
kepler = True corresponds to an alternative potential as proposed in this paper:  https://arxiv.org/pdf/2008.04922.pdf
Since working with the true black hole potential didn't work out to be succesful in this part of the project (comparing the orbits), 
we tried using the alternativ one, and it reproduces physically realstic results. 

In the file tde_fien_jurek.py we move on to creating different particles in the (this time true) potential of a black hole, when the star (that will disappear, once turned into particles)
crosses a certain threshold (is less than a certain distance away from the black hole - this distance will be closely related to the tidal radius)
'''

##############################################################################################################
#                                       Kepler orbit or BH orbit? 
kepler = False     
# set to True for Newtonian potential, False for BH potential (the alternative potential from the paper)
##############################################################################################################


# Global Constants like G, c, the mass of the star, the mass of the black hole, the radius of the star

G = constants.G
c = constants.speed_of_light

m = 1.989 * 10**30              # mass of the sun in kg
M = 8*10**6 * m                   # mass of BH

r_star = 696340*10**3


# Black hole quantities

# tidal radius
r_tidal = r_star * (M/m)**(1/3)

# Schwarzschild radius
r_ss   = 2*G*M/c**2

# innermost stable circular orbit
r_isco = 3*r_ss

# stable photon orbit
r_photon = 1.5*r_ss

'''
Initial conditions
'''

x0 = 0
y0 = r_tidal#r_isco*0.5#3*10**9



def orbit_velo():
    '''
    returns the orbit velocities depending on potential
    '''
    if kepler == True:
        return np.sqrt(G*M/y0)
    if kepler == False:
        return np.sqrt(G*M/y0 + 3*G**2*M**2/(y0**2*c**2))

v_x0 = orbit_velo()
v_y0 = 0


print('vy_0 is ', v_x0/c*100, '% of the speed of light.')




# define a timescale on which the orbiting happens
# How much time do we want to calculate the trajectory for? 
t = np.linspace(0,100000,20000)

# where the time at time zero is the first element in this period (which is zero)
t0 = t[0]


'''
Since the Runge-Kutta algorithm is only heplful for solving differential equations of first order, 
we split our two differential equations of second order (one for x, one for y) up in two equations of motion of first order each.
dot(x) = v
dot(v) = - grad(Potential)
'''

def v_x_(t,x,y,v_x,v_y):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -G*M*x / r**3 
    if kepler == False:
        return -G*M*x / r**3 - 3*(G*M/c)**2*x / r**4

def v_y_(t,x,y,v_x,v_y):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -G*M*y / r**3 
    if kepler == False:
        return -G*M*y / r**3 - 3*(G*M/c)**2*y / r**4

def x_(t,x,y,v_x,v_y):
    return v_x

def y_(t,x,y,v_x,v_y):
    return v_y





def runge_kutta(x0, v_x0, y0, v_y0, h = 1):
    '''
    This function solves the differential equations (the equation of motion) of second order
    for the x- and y-component of a particle. At each step, we evaluate v[i] for x and y which corresponds to their derivatives 
    as well as the derivatives dot(v) for x and y which corresponds to their second order derivatives. 
    Since the results of the (i+1)-th step depends on the values of x, y, v_x and v_y for the i-th step, for each step we calculate these 4 values
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

        # update x, y, v_x, v_y according to Runge-Kutta algorithm
        x[i+1]      = x[i]      + (k1x+2*k2x+2*k3x+k4x)*h/6             
        v_x[i+1]    = v_x[i]    + (k1v_x+2*k2v_x+2*k3v_x+k4v_x)*h/6
        y[i+1]      = y[i]      + (k1y+2*k2y+2*k3y+k4y)*h/6
        v_y[i+1]    = v_y[i]    + (k1v_y+2*k2v_y+2*k3v_y+k4v_y)*h/6


        # If the star is passing the Schwarzschild radius, it disappears. 
        # The loop breaks and the arrays that contain the values for x, y, vx, vy are shortened to the length array[:i]
        if np.sqrt(x[i]**2+y[i]**2) <= r_ss:
            x[i] = 0
            y[i] = 0
            x = x[:i]
            y = y[:i]
            v_x = v_x[:i]
            v_y = v_y[:i]
            break
    return x, y, v_x, v_y


# solve the eom for the x and y coordinate of the star for the initial conditions given above.
x, y, vx, vy = runge_kutta(x0, v_x0, y0, v_y0)

r = np.sqrt(x**2+ y**2)


def number_rotations():
    '''
    calculates index for that one rotation is done, prints gain in distance after one rotation
    to compare the instability of the orbits to each other
    '''
    for i in range(10,len(x)-1):
        if x[i]<x[i+1] and x[i] > x[i-1] and y[i] > y[i-1] and y[i] > y[i+1]:
            print('index of interest: ', i)
            print('The difference in distance from the BH after one rotation: ', (np.sqrt(x[i]**2 + y[i]**2)-np.sqrt(x[0]**2 + y[0]**2))/10**4, 'e4')
            
    print('There was no full rotation simulated.')
    return None
    
number_rotations()


###################################################################################################################################################
#                                                      Plotting
###################################################################################################################################################

'''
This step creates an animation of the orbit of the star in the potential of the black hole 
(which we also approximate to be a newtonian potential which is unphysical; it serves illustrating purposes)
'''

# create an array for the animation
array_x_y = np.array([x, y])
num = int(len(t)*100)


def animation_(i):
    ax.clear()  
    # Updating Trajectory Line
    ax.plot(array_x_y[0, :(i+1)*100], array_x_y[1, :(i+1)*100], c='red')
    # Updating Point Location 
    ax.scatter(array_x_y[0, (i+1)*100], array_x_y[1, (i+1)*100], 
               c='yellow', marker='o', s = 100)


    bh = plt.Circle((0,0), r_ss, color = 'black', fill = True, lw = 2, label = 'Black hole')
    isco_radius = plt.Circle((0,0), r_isco, color = 'purple', fill = False, lw = 2, label = 'Isco radius')
    tidal_radius = plt.Circle((0,0), r_tidal, color = 'green', fill = False, lw = 2, label = 'Tidal radius')
    ax.add_patch(isco_radius)
    ax.add_patch(bh)
    ax.add_patch(tidal_radius)

    ax.set_xlim([-y0*2, y0*2])
    ax.set_ylim([-y0*2, y0*2])

    ax.set_title('Trajectory of the star')
    ax.set_xlabel('x-coordinate')
    ax.set_ylabel('y-coordinate')
    ax.legend()
    


fig = plt.figure()
ax = plt.axes()

orbit_animation = animation.FuncAnimation(fig, animation_, interval=1, frames=num)


plt.show()

'''
Simple 2-d plot of stars trajectory
'''
# fig, ax = plt.subplots()


# xla = plt.xlabel('X-coordinate of object')
# yla = plt.ylabel('Y-coordinate of object')
# end = plt.scatter(x[-1], y[-1], label = 'End', s = 100, color = 'yellow')

# plt.plot(x, y, label = 'Movement of star', color = 'red')
# bh = plt.Circle((0,0), r_ss, color = 'black', fill = True, lw = 2, label = 'Black hole')
# isco_radius = plt.Circle((0,0), r_isco, color = 'purple', fill = False, lw = 2, label = 'Isco radius')
# tidal_radius = plt.Circle((0,0), r_tidal, color = 'green', fill = False, lw = 2, label = 'Tidal radius')
# ax.add_patch(isco_radius)
# ax.add_patch(bh)
# ax.add_patch(tidal_radius)

# plt.xlim(-y0*4, y0*4)
# plt.ylim(-y0*4, y0*4)
# if kepler == True:
#     plt.title('Star motion in Newtonian potential')
# if kepler == False: 
#     plt.title('Star motion in black hole potential')
# plt.legend(loc = 'upper right')

# plt.show()




