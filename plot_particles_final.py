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
M = 8*10**6 * m                   # mass of BH

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
escape = np.load('escape.npy')

# This many particles escaped:
print('This many particles escaped:' ,np.sum(escape))


'''
Animating the star turning into particles
'''
Animation = True


###################################################################################################################################################
#               Working Animation of the orbit
###################################################################################################################################################

if Animation == True:
    data_star = np.array([x, y])
    data_layer1 = np.array([x_p[0,0], y_p[0,0]])
    data_layer2 = np.array([x_p[0,1], y_p[0,1]])
    data_layer3 = np.array([x_p[0,2], y_p[0,2]])
    data_layer4 = np.array([x_p[0,3], y_p[0,3]])
    data_layer5 = np.array([x_p[0,4], y_p[0,4]])
    data_layer6 = np.array([x_p[0,5], y_p[0,5]])
    data_layer7 = np.array([x_p[0,6], y_p[0,6]])
    data_layer8 = np.array([x_p[0,7], y_p[0,7]])
    numDataPoints = int(len(t)*100)

    # time_offset = np.zeros(len(data_star))
    data_layer1 = np.hstack((data_star, data_layer1))
    data_layer2 = np.hstack((data_star, data_layer2))
    data_layer3 = np.hstack((data_star, data_layer3))
    data_layer4 = np.hstack((data_star, data_layer4))
    data_layer5 = np.hstack((data_star, data_layer5))
    data_layer6 = np.hstack((data_star, data_layer6))
    data_layer7 = np.hstack((data_star, data_layer7))
    data_layer8 = np.hstack((data_star, data_layer8))


    def animate_func(num):
        # when do we want the animation to stop?
        if num >= 100:
            line_ani.event_source.stop()
            return None
        ax.clear()   # each iteration, there's a new plot
        speed = 30

        ax.plot(data_layer1[0, :(num+1)*speed], data_layer1[1, :(num+1)*speed], label = 'P 1')

        ax.plot(data_layer2[0, :(num+1)*speed], data_layer2[1, :(num+1)*speed], label = 'P 2')

        ax.plot(data_layer3[0, :(num+1)*speed], data_layer3[1, :(num+1)*speed], label = 'P 3')

        ax.plot(data_layer4[0, :(num+1)*speed], data_layer4[1, :(num+1)*speed], label = 'P 4')

        ax.plot(data_layer5[0, :(num+1)*speed], data_layer5[1, :(num+1)*speed], label = 'P 5')

        # ax.plot(data_layer6[0, :(num+1)*speed], data_layer6[1, :(num+1)*speed], label = 'P 6')

        # ax.plot(data_layer7[0, :(num+1)*speed], data_layer7[1, :(num+1)*speed], label = 'P 7')

        # ax.plot(data_layer8[0, :(num+1)*speed], data_layer8[1, :(num+1)*speed], label = 'P 8')

        bh = plt.Circle((0,0), r_ss, color = 'black', fill = True, lw = 2, label = 'Black hole')
        isco_radius = plt.Circle((0,0), r_isco, color = 'purple', fill = False, lw = 2, label = 'Isco radius')
        tidal_radius = plt.Circle((0,0), r_tidal, color = 'green', fill = False, lw = 2, label = 'Tidal radius')
        ax.add_patch(isco_radius)
        ax.add_patch(bh)
        ax.add_patch(tidal_radius)

        # Setting Axes Limits
        ax.set_xlim([-y0*2, y0*2])
        ax.set_ylim([-y0*2, y0*2])


        # Adding Figure Labels
        ax.set_title('Tidal disruption event')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.legend()
        


    fig = plt.figure()
    ax = plt.axes()

    line_ani = animation.FuncAnimation(fig, animate_func, interval=1,   
                                    frames=numDataPoints)


#     line_ani.save('particles_anim_1.gif', fps=10)

    # f = r"particles_animation.gif" 
    # writergif = animation.PillowWriter(fps=30) 
    # line_ani.save(f, writer='imagemagick', fps = 30)


    plt.show()

if Animation == False:
    fig, ax = plt.subplots()

    # plt.plot(x,y, label = 'Star', color = 'black')
    # plt.plot((x_p[4,0]),(y_p[4,0]),label = 'innermost layer')
    # plt.plot((x_p[0,0]),(y_p[0,0]),label = 'outermost layer')
    # plt.plot((x_p[1,0]),(y_p[1,0]),label = '4th layer')
    # plt.plot((x_p[2,0]),(y_p[2,0]),label = '3rd layer')
    # plt.plot((x_p[3,0]),(y_p[3,0]),label = '2nd layer')

    # outermost layer
    plt.plot((x_p[0,0]),(y_p[0,0]),label = 'L 1, P 1, E '+str(escape[0,0]))
    plt.plot((x_p[0,1]),(y_p[0,1]),label = 'L 1, P 2, E '+str(escape[0,1]))
    plt.plot((x_p[0,2]),(y_p[0,2]),label = 'L 1, P 3, E '+str(escape[0,2]))
    plt.plot((x_p[0,3]),(y_p[0,3]),label = 'L 1, P 4, E '+str(escape[0,3]))
    plt.plot((x_p[0,4]),(y_p[0,4]),label = 'L 1, P 5, E '+str(escape[0,4]))
    plt.plot((x_p[0,5]),(y_p[0,5]),label = 'L 1, P 6, E '+str(escape[0,5]))
    plt.plot((x_p[0,6]),(y_p[0,6]),label = 'L 1, P 7, E '+str(escape[0,6]))
    plt.plot((x_p[0,7]),(y_p[0,7]),label = 'L 1, P 8, E '+str(escape[0,7]))

    # 2nd layer
#     plt.plot((x_p[1,0]),(y_p[1,0]),label = 'L 2, P 1'+str(escape[0,0]))
#     plt.plot((x_p[1,1]),(y_p[1,1]),label = 'L 2, P 2'+str(escape[0,0]))
#     plt.plot((x_p[2,2]),(y_p[1,2]),label = 'L 2, P 3'+str(escape[0,0]))
#     plt.plot((x_p[1,3]),(y_p[1,3]),label = 'L 2, P 4'+str(escape[0,0]))
#     plt.plot((x_p[1,4]),(y_p[1,4]),label = 'L 2, P 5, E '+str(escape[1,4]))
#     plt.plot((x_p[1,5]),(y_p[1,5]),label = 'L 2, P 6'+str(escape[0,0]))
#     plt.plot((x_p[1,6]),(y_p[1,6]),label = 'L 2, P 7'+str(escape[0,0]))
#     plt.plot((x_p[1,7]),(y_p[1,7]),label = 'L 2, P 8'+str(escape[0,0]))

    # 3rd layer
#     plt.plot((x_p[2,0]),(y_p[2,0]),label = 'L 3, P 1'+str(escape[0,0]))
#     plt.plot((x_p[2,1]),(y_p[2,1]),label = 'L 3, P 2'+str(escape[0,0]))
#     plt.plot((x_p[2,2]),(y_p[2,2]),label = 'L 3, P 3'+str(escape[0,0]))
    # plt.plot((x_p[2,3]),(y_p[2,3]),label = 'L 3, P 4'+str(escape[0,0]))
#     plt.plot((x_p[2,4]),(y_p[2,4]),label = 'L 3, P 5, E'+str(escape[2,4]))
    # plt.plot((x_p[2,5]),(y_p[2,5]),label = 'L 3, P 6'+str(escape[0,0]))
    # plt.plot((x_p[2,6]),(y_p[2,6]),label = 'L 3, P 7'+str(escape[0,0]))
    # plt.plot((x_p[2,7]),(y_p[2,7]),label = 'L 3, P 8+str(escape[0,0])')

    # 4th layer
#     plt.plot((x_p[3,0]),(y_p[3,0]),label = 'L 4, P 1'+str(escape[0,0]))
#     plt.plot((x_p[3,1]),(y_p[3,1]),label = 'L 4, P 2'+str(escape[0,0]))
#     plt.plot((x_p[3,2]),(y_p[3,2]),label = 'L 4, P 3'+str(escape[0,0]))
    # plt.plot((x_p[3,3]),(y_p[3,3]),label = 'L 4, P 4'+str(escape[0,0]))
#     plt.plot((x_p[3,4]),(y_p[3,4]),label = 'L 4, P 5, E'+str(escape[3,4]))
    # plt.plot((x_p[3,5]),(y_p[3,5]),label = 'L 4, P 6'+str(escape[0,0]))
    # plt.plot((x_p[3,6]),(y_p[3,6]),label = 'L 4, P 7'+str(escape[0,0]))
    # plt.plot((x_p[3,7]),(y_p[3,7]),label = 'L 4, P 8'+str(escape[0,0]))
    # 5th layer
#     plt.plot((x_p[4,0]),(y_p[4,0]),label = 'L 5, P 1'+str(escape[0,0]))
#     plt.plot((x_p[4,1]),(y_p[4,1]),label = 'L 5, P 2'+str(escape[0,0]))
#     plt.plot((x_p[4,2]),(y_p[4,2]),label = 'L 5, P 3'+str(escape[0,0]))
    # plt.plot((x_p[4,3]),(y_p[4,3]),label = 'L 5, P 4'+str(escape[0,0]))
#     plt.plot((x_p[4,4]),(y_p[4,4]),label = 'L 5, P 5, E'+str(escape[4,4]))
    # plt.plot((x_p[4,5]),(y_p[4,5]),label = 'L 5, P 6'+str(escape[0,0]))
    # plt.plot((x_p[4,6]),(y_p[4,6]),label = 'L 5, P 7'+str(escape[0,0]))
    # plt.plot((x_p[4,7]),(y_p[4,7]),label = 'L 5, P 8'+str(escape[0,0]))

    xla = plt.xlabel('X-coordinate of object')
    yla = plt.ylabel('Y-coordinate of object')
    end = plt.scatter(x[-1], y[-1], label = 'End', s = 100, color = 'yellow')

    bh = plt.Circle((0,0), r_ss, color = 'black', fill = True, lw = 2, label = 'Black hole')
    isco_radius = plt.Circle((0,0), r_isco, color = 'purple', fill = False, lw = 2, label = 'Isco radius')
    tidal_radius = plt.Circle((0,0), r_tidal, color = 'green', fill = False, lw = 2, label = 'Tidal radius')
    ax.add_patch(isco_radius)
    ax.add_patch(bh)
    ax.add_patch(tidal_radius)

    plt.xlim(-y0*2, y0*2)
    plt.ylim(-y0*2, y0*2)

    plt.title('Star motion in Black hole potential')
    plt.legend(loc = 'upper right')

    plt.show()
