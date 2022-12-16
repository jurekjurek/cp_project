import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants
from sympy import *

G = constants.G
c = constants.speed_of_light

# We need a supermassive Black hole!!!  

m = 1.989 * 10**30              # mass of the sun in kg
M = 4*10**6 * m                   # mass of BH

r_star = 696340*10**3

A = G*M
# B = L**2 / m**2             # 0.005
# C = 3*G*M*L**2 / (m*c**2)   # 0.8


# things for the Black hole

r_ss   = 2*A/c**2
r_isco = 3*r_ss

# stable orbit for photons
r_photon = 1.5 * r_ss

# tidal radius
r_tidal = r_star * (M/m)**(1/3)


'''
Initial conditions
'''
x0 =    -r_tidal*1.2 
y0 =    r_tidal*1.5

# v_x0 =  10**8
v_y0 =  0

# orbit velo for kepler potential
# v_x0 = np.sqrt(G*M/y0)

# orbit velo for BH
# it actually works...
v_x0 = np.sqrt(G*M/y0 + 3*G**2*M**2/(y0**2*c**2))*0.7

t = np.linspace(0,1000000000,20000)





x = np.load('x.npy')
y = np.load('y.npy')
vx = np.load('vx.npy')
vy = np.load('vy.npy')
x_p = np.load('x_p.npy')
y_p = np.load('y_p.npy')
vx_p = np.load('vx_p.npy')
vy_p = np.load('vy_p.npy')


'''
Plotting 
'''

# fig, ax = plt.subplots()

# plt.plot(x,y, label = 'Star', color = 'black')
# plt.plot(x_p[4,0], y_p[4,0], label = 'innermost layer')
# plt.plot(x_p[0,0], y_p[0,0], label = 'outermost layer')
# plt.plot(x_p[1,0], y_p[1,0], label = '4th layer')
# plt.plot(x_p[2,0], y_p[2,0], label = '3rd layer')
# plt.plot(x_p[3,0], y_p[3,0], label = '2nd layer')
# xla = plt.xlabel('X-coordinate of the particle')
# yla = plt.ylabel('Y-coordinate of the particle')
# # beg = plt.scatter(x[0], y[0], label = 'beginning')
# # end = plt.scatter(x[-1], y[-1], label = 'end')          # last element of x and y
# bh = plt.scatter(0,0, label= 'Black Hole', color = 'black')#, s = 300)
# ss = plt.scatter(0,-r_ss, label = 'Schwarzschild radius', s = 5)
# isco = plt.scatter(0,-r_isco, label = 'isco radius', s = 5)
# # tidal = plt.scatter(0,-r_tidal, label = 'Tidal radius for star')
# # if we maybe want to do it more fancy:
# #
# circle = plt.Circle((0,0), r_tidal, color = 'red', fill = False, lw = 2, label = 'Tidal radius')

# ax.add_patch(circle)

# plt.xlim(-y0*2, y0*2)
# plt.ylim(-y0*2, y0*2)
# plt.title('Motion of the particle in the Kepler field - solved with RK4')

# plt.legend(loc = 'upper right')

# plt.show()






'''
Animating the star turning into particles
'''

###################################################################################################################################################
#               Working Animation of the orbit
###################################################################################################################################################

data_star = np.array([x, y])
data_layer1 = np.array([x_p[0,0], y_p[0,0]])
data_layer2 = np.array([x_p[1,1], y_p[1,1]])
data_layer3 = np.array([x_p[2,2], y_p[2,2]])
data_layer4 = np.array([x_p[3,3], y_p[3,3]])
data_layer5 = np.array([x_p[4,4], y_p[4,4]])
numDataPoints = int(len(t)*100)

# time_offset = np.zeros(len(data_star))
data_layer1 = np.hstack((data_star, data_layer1))
data_layer2 = np.hstack((data_star, data_layer2))
data_layer3 = np.hstack((data_star, data_layer3))
data_layer4 = np.hstack((data_star, data_layer4))
data_layer5 = np.hstack((data_star, data_layer5))

def animate_func(num):
    if num >= 100:
        line_ani.event_source.stop()
        return None
    ax.clear()  # Clears the figure to update the line, point,   
                # title, and axes
    # Updating Trajectory Line (num+1 due to Python indexing)

    # Animating the star movement:
    # if np.sqrt(x[num]**2+y[num]**2) <= r_tidal:
    # ax.plot(data_star[0, :(num+1)*100], data_star[1, :(num+1)*100], label = 'star', c='black')
    # # Updating Point Location 
    # ax.scatter(data_star[0, num], data_star[1, num],  c='black', marker='o', s = 1)

    if True: 
        ax.plot(data_layer1[0, :(num+1)*100], data_layer1[1, :(num+1)*100], label = 'outermost')#, c='black')
        # Updating Point Location 
        ax.scatter(data_layer1[0, num], data_star[1, num], 
                c='black', marker='o', s = 1)

        ax.plot(data_layer2[0, :(num+1)*100], data_layer2[1, :(num+1)*100], label = '4')#, c='black')
        # Updating Point Location 
        ax.scatter(data_layer2[0, num], data_star[1, num], 
                c='black', marker='o', s = 1)

        ax.plot(data_layer3[0, :(num+1)*100], data_layer3[1, :(num+1)*100], label = '3')#, c='black')
        # Updating Point Location 
        ax.scatter(data_layer3[0, num], data_star[1, num], 
                c='black', marker='o', s = 1)

        ax.plot(data_layer4[0, :(num+1)*100], data_layer4[1, :(num+1)*100], label = '2')#, c='black')
        # Updating Point Location 
        ax.scatter(data_layer4[0, num], data_star[1, num], 
                c='black', marker='o', s = 1)

        ax.plot(data_layer5[0, :(num+1)*100], data_layer5[1, :(num+1)*100], label = 'innermost')#, c='black')
        # Updating Point Location 
        ax.scatter(data_layer5[0, num], data_star[1, num], 
                c='black', marker='o', s = 1)

    ax.scatter(0,0, label= 'Black Hole', color = 'black')#, s = 300)
    ax.scatter(0,-r_isco, label = 'isco radius')
    ax.scatter(0,-r_ss, label = 'Schwarzschild radius')

    circle = plt.Circle((0,0), r_tidal, color = 'red', fill = False, lw = 2, label = 'Tidal radius')
    ax.add_patch(circle)
    # Setting Axes Limits
    ax.set_xlim([-y0*2, y0*2])
    ax.set_ylim([-y0*2, y0*2])


    # Adding Figure Labels
    ax.set_title('Trajectory \nTime = ' + str(np.round(t[num]/(3600*24),    
                 decimals=2)) + ' days')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend()
    


fig = plt.figure()
ax = plt.axes()

line_ani = animation.FuncAnimation(fig, animate_func, interval=1,   
                                   frames=numDataPoints)


# line_ani.save('particles_anim.gif', fps=10)

# f = r"particles_animation.gif" 
# writergif = animation.PillowWriter(fps=30) 
# line_ani.save(f, writer='imagemagick', fps = 30)


plt.show()


