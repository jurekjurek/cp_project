import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants

'''
This file creates the trajectory of earths orbit around the sun
'''


kepler = True


# Global Constants like G, c, the mass of the star, the mass of the black hole, the radius of the star

G = constants.G
c = constants.speed_of_light

if kepler == True:
    m = 5.972*10**24                # mass of earth
    M = 1.989 * 10**30              # mass of sun

r_star = 696340*10**3


'''
Initial conditions
'''

x0 = 0

if kepler == True:
    y0 =    147200000               # distance sun-earth


def orbit_velo():
    '''
    returns the orbit velocitie for the earth in suns potential
    '''
    if kepler == True:
        return np.sqrt(G*M/y0)

v_x0 = orbit_velo()
v_y0 = 0

# escape velocity for newtonian potential, given initial conditions x = 0, so velocity has only component in x-direction when orbiting
# v_x0 = 0
# v_y0 = np.sqrt(A*2/y0)*0.7

print('vy_0 is ', v_x0/c*100, '% of the speed of light.')



# define a timescale on which the orbiting happens
# How much time do we want to calculate the trajectory for? 
t = np.linspace(0,100000,40000)

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

def v_y_(t,x,y,v_x,v_y):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -G*M*y / r**3 

def x_(t,x,y,v_x,v_y):
    return v_x

def y_(t,x,y,v_x,v_y):
    return v_y





def runge_kutta(x0, v_x0, y0, v_y0, h = 0.2):
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

    return x, y, v_x, v_y


# solve the eom for the x and y coordinate of the star for the initial conditions given above.
x, y, vx, vy = runge_kutta(x0, v_x0, y0, v_y0)

r = np.sqrt(x**2+ y**2)

'''
animation
'''


array_x_y = np.array([x, y])
num = int(len(t)*100)


def animation_(i):
    ax.clear()
    ax.plot(array_x_y[0, :(i+1)*300], array_x_y[1, :(i+1)*300], c='green', label = 'Earth')
    ax.scatter(array_x_y[0, (i+1)*300], array_x_y[1, (i+1)*300], 
               c='blue', marker='o', s = 100)

    ax.scatter(0,0, label= 'Sun', color = 'yellow', s = 300)

    # Setting Axes Limits
    ax.set_xlim([-y0*2, y0*2])
    ax.set_ylim([-y0*2, y0*2])

    ax.set_title('Movement of Earth in the gravitational potential of the sun\n')
    ax.set_xlabel('x-coordinate')
    ax.set_ylabel('y-coordinate')
    ax.legend()
    


fig = plt.figure()
ax = plt.axes()

line_ani = animation.FuncAnimation(fig, animation_, interval=1, frames=num)

# line_ani.save('earth_h02.gif', fps=10)


plt.show()

'''
2d plot
'''
# fig, ax = plt.subplots()


# xla = plt.xlabel('X-coordinate of object')
# yla = plt.ylabel('Y-coordinate of object')
# # beg = plt.scatter(x[0], y[0], label = 'beginning')
# end = plt.scatter(x[-1], y[-1], label = 'End', color = 'blue', s = 100)

# plt.plot(x, y, label= 'Earth', color = 'green')
# sun = plt.scatter(0,0, label= 'Sun', color = 'yellow', s = 300)
    

# plt.xlim(-y0*2, y0*2)
# plt.ylim(-y0*2, y0*2)
# plt.title('Movement of earth in the gravitational potential of the sun\n')
# plt.legend(loc = 'upper right')
# plt.savefig('earth_h02.png')
# plt.show()




