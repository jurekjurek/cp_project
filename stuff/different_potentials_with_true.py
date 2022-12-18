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
One can, varying the factor kepler between True, False and 'alternative', explore the effects of the different potentials on the 
orbiting behaviour of the star. 
kepler = False corresponds to a purely newtonian potential, 
kepler = True corresponds to the true potential of a black hole
and kepler = 'alternative' corresponds to an alternative potential as proposed in this paper:  https://arxiv.org/pdf/2008.04922.pdf
Since working with the true black hole potential didn't work out to be succesful, we tried using the alternativ one, and the alternative one
reproduces physically realstic results. 

In the file .... we move on to creating different particles in the potential of a black hole, when the star (that will disappear, once turned into particles)
crosses a certain threshold (is less than a certain distance away from the black hole - this distance will be closely related to the tidal radius)
'''

##############################################################################################################
#                                       Kepler orbit or BH orbit? 
kepler = True#'alternative' #True#False              
# set to True for Newtonian potential, False for BH potential, which does not seem to work at all
# and to 'alternative' for the alternative potential from the paper
##############################################################################################################


# Global Constants like G, c, the mass of the star, the mass of the black hole, the radius of the star

G = constants.G
c = constants.speed_of_light

m = 1.989 * 10**30              # mass of the sun in kg
M = 10**6 * m                   # mass of BH

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
y0 = 2*r_isco#3*10**9



def orbit_velo():
    '''
    returns the orbit velocities depending on potential
    '''
    if kepler == True:
        return np.sqrt(G*M/y0)
    if kepler == 'alternative' or kepler == False:
        return np.sqrt(G*M/y0 + 3*G**2*M**2/(y0**2*c**2))

v_x0 = orbit_velo()
v_y0 = 0


print('vy_0 is ', v_x0/c*100, '% of the speed of light.')



def A_B_C(x,y, m):
    '''
    function that calculates the angular momentum for given x, y
    L is constant over time
    Since at point t = 0, we have only a velocity component in x-direction, we calculate the L(0) = L(t).
    Then the quantities A, B, and C are calculated that we need for the true potential of a black hole which depends on the angular momentum
    '''
    x = float(x)
    y = float(y)
    r = np.sqrt((x**2+y**2))
    L = r*m*v_x0
    A = G*M
    B = L**2 / m**2                     # 0.005
    C = 3*G*M*L**2 / (m**2*c**2)        # 0.8
    return A, B, C

A, B, C = A_B_C(x0, y0, m)



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
        return -A*x / r**3 
    if kepler == 'alternative':
        return -A*x / r**3 - 3*(G*M/c)**2*x / r**4
    if kepler == False:
        return -A*x / r**3 + B*x / r**4 - C*x / r**5

def v_y_(t,x,y,v_x,v_y):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -A*y / r**3 
    if kepler == 'alternative':
        return -A*y / r**3 - 3*(G*M/c)**2*y / r**4
    if kepler == False: 
        return -A*y / r**3 + B*y / r**4 - C*y / r**5

def x_(t,x,y,v_x,v_y):
    return v_x

def y_(t,x,y,v_x,v_y):
    return v_y





def runge_kutta(x0, v_x0, y0, v_y0, h = 0.08):
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

        # update x, y, v_x, v_y
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

# how much did we orbit away during one rotation?
# print('Is the orbit unstable? The r in the beginning was: ', np.sqrt(x[0]**2 + y[0]**2), 'The r in the end is:', np.sqrt(x[-1]**2 + y[-1]**2))

def number_rotations():
    '''
    returns index for that one rotation is done
    '''
    count = 0
    for i in range(1,len(x)):
        # if x[i] == 0 and count == 0: 
        #     count += 1
        # if x[i] == 0 and count == 1:
        #     return(i)
        if x[i]<x[i+1] and x[i] > x[i-1]:
            return i

i = number_rotations()
# print('The difference between the distance from the BH in the beginning and end is: ', (np.sqrt(x[-1]**2 + y[-1]**2)-np.sqrt(x[0]**2 + y[0]**2))/10**8, 'e8')
print('The difference between the distance from the BH in the beginning and end is: ', (np.sqrt(x[i]**2 + y[i]**2)-np.sqrt(x[0]**2 + y[0]**2))/10**4, 'e4')

###################################################################################################################################################
#                                                      Working Animation of the orbit
###################################################################################################################################################

'''
This step creates an animation of the orbit of the star in the potential of the black hole 
(which we also approximate to be a newtonian potential which is unphyical; it serves illustrating purposes)
'''


dataSet = np.array([x, y])
numDataPoints = int(len(t)*100)


def animate_func(num):
    # if np.sqrt(dataSet[0,num+1]**2+ dataSet[1,num+1]**2) <= r_Earth:
    #     line_ani.event_source.stop()
    #     return None
    ax.clear()  # Clears the figure to update the line, point,   
                # title, and axes
    # ax.plot_su6rface(X_Earth, Y_Earth, Z_Earth, color='green', alpha=0.7)
    # Updating Trajectory Line (num+1 due to Python indexing)
    ax.plot(dataSet[0, :(num+1)*100], dataSet[1, :(num+1)*100], c='black')
    # Updating Point Location 
    ax.scatter(dataSet[0, (num+1)*100], dataSet[1, (num+1)*100], 
               c='black', marker='o', s = 1)


    ax.scatter(0,0, label= 'Black Hole', color = 'black')#, s = 300)
    ax.scatter(0,-r_isco, label = 'isco radius')
    ax.scatter(0,-r_ss, label = 'Schwarzschild radius')
    circle = plt.Circle((0,0), r_tidal, color = 'red', fill = False, lw = 2, label = 'Tidal radius')
    ax.add_patch(circle)

    # Setting Axes Limits
    ax.set_xlim([-y0*2, y0*2])
    ax.set_ylim([-y0*2, y0*2])
    # ax.set_zlim3d([-10000, 10000])


    # fancier way of limiting, but it doesn't work with the animation

    # xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(),      
    #                ax.get_zlim3d()]).T
    # XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])
    # ax.set_xlim3d(XYZlim)
    # ax.set_ylim3d(XYZlim)
    # ax.set_zlim3d(XYZlim * 3/4)

    # Adding Figure Labels
    ax.set_title('Trajectory \nTime = ' + str(np.round(t[num],    
                 decimals=2)) + ' sec')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend()
    


fig = plt.figure()
ax = plt.axes()

line_ani = animation.FuncAnimation(fig, animate_func, interval=1,   
                                   frames=numDataPoints)




plt.show()

###################################################################################################################################################
#               Simple 2d-plot of trajectory
###################################################################################################################################################


fig, ax = plt.subplots()


xla = plt.xlabel('X-coordinate of object')
yla = plt.ylabel('Y-coordinate of object')
beg = plt.scatter(x[0], y[0], label = 'beginning')
end = plt.scatter(x[-1], y[-1], label = 'end')

plt.plot(x, y, label = 'Movement of star')
bh = plt.scatter(0,0, label= 'Black Hole', color = 'black')#, s = 300)
ss = plt.scatter(0,-r_ss, label = 'Schwarzschild radius')
isco = plt.scatter(0,-r_isco, label = 'isco radius')
circle = plt.Circle((0,0), r_tidal, color = 'red', fill = False, lw = 2, label = 'Tidal radius')
ax.add_patch(circle)

plt.xlim(-y0*4, y0*4)
plt.ylim(-y0*4, y0*4)
if kepler == True:
    plt.title('Star motion in Newtonian potential')
else: 
    plt.title('Star motion in black hole potential')
plt.legend(loc = 'upper right')

plt.show()




