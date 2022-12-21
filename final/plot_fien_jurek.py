import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants

'''
This file simply animates or plots the trajectory of the star and the particles it turns into. 
So this file illustrates the TDE. 
'''



G = constants.G
c = constants.speed_of_light


m = 1.989 * 10**30              # mass of the sun in kg
M = 8*10**6 * m                   # mass of BH

r_star = 696340*10**3

A = G*M


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
Animating (=True) or plotting (=False) the star turning into particles
'''
Animation = True

# Animation

if Animation == True:
    layer = 1
    star = np.array([x, y])
    layer1 = np.array([x_p[layer,0], y_p[layer,0]])
    layer2 = np.array([x_p[layer,1], y_p[layer,1]])
    layer3 = np.array([x_p[layer,2], y_p[layer,2]])
    layer4 = np.array([x_p[layer,3], y_p[layer,3]])
    layer5 = np.array([x_p[layer,4], y_p[layer,4]])
    layer6 = np.array([x_p[layer,5], y_p[layer,5]])
    layer7 = np.array([x_p[layer,6], y_p[layer,6]])
    layer8 = np.array([x_p[layer,7], y_p[layer,7]])
    num = int(len(t)*100)

    # time_offset = np.zeros(len(star))
    layer1 = np.hstack((star, layer1))
    layer2 = np.hstack((star, layer2))
    layer3 = np.hstack((star, layer3))
    layer4 = np.hstack((star, layer4))
    layer5 = np.hstack((star, layer5))
    layer6 = np.hstack((star, layer6))
    layer7 = np.hstack((star, layer7))
    layer8 = np.hstack((star, layer8))


    def animation_(i):
        # when do we want the animation to stop?
        if i >= 100:
            orbit_animation.event_source.stop()
            return None
        ax.clear()   # each iteration, there's a new plot
        speed = 30      # how fast is the animation? 

        ax.plot(layer1[0, :(i+1)*speed], layer1[1, :(i+1)*speed], label = 'P 1, E ' + str(escape[layer,0]))

        ax.plot(layer2[0, :(i+1)*speed], layer2[1, :(i+1)*speed], label = 'P 2, E ' + str(escape[layer,1]))

        ax.plot(layer3[0, :(i+1)*speed], layer3[1, :(i+1)*speed], label = 'P 3, E ' + str(escape[layer,2]))

        ax.plot(layer4[0, :(i+1)*speed], layer4[1, :(i+1)*speed], label = 'P 4, E ' + str(escape[layer,3]))

        ax.plot(layer5[0, :(i+1)*speed], layer5[1, :(i+1)*speed], label = 'P 5, E ' + str(escape[layer,4]))

        ax.plot(layer6[0, :(i+1)*speed], layer6[1, :(i+1)*speed], label = 'P 6, E ' + str(escape[layer,5]))

        ax.plot(layer7[0, :(i+1)*speed], layer7[1, :(i+1)*speed], label = 'P 7, E ' + str(escape[layer,6]))

        ax.plot(layer8[0, :(i+1)*speed], layer8[1, :(i+1)*speed], label = 'P 8, E ' + str(escape[layer,7]))

        bh                  = plt.Circle((0,0), r_ss, color = 'black', fill = True, lw = 2, label = 'Black hole')
        isco_radius         = plt.Circle((0,0), r_isco, color = 'purple', fill = False, lw = 2, label = 'Isco radius')
        tidal_radius        = plt.Circle((0,0), r_tidal, color = 'green', fill = False, lw = 2, label = 'Tidal radius')
        ax.add_patch(isco_radius)
        ax.add_patch(bh)
        ax.add_patch(tidal_radius)

        ax.set_xlim([-r_tidal*2, r_tidal*2])
        ax.set_ylim([-r_tidal*2, r_tidal*2])

        ax.set_title('Tidal disruption event')
        ax.set_xlabel('x-coordinate')
        ax.set_ylabel('y-coordinate')
        ax.legend()
        


    fig = plt.figure()
    ax = plt.axes()

    orbit_animation = animation.FuncAnimation(fig, animation_, interval=1, frames=num)

    # save animation as .gif
    # orbit_animation.save('name.gif', fps=10)

    plt.show()

# Plot

if Animation == False:
    fig, ax = plt.subplots()

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
#     plt.plot((x_p[1,0]),(y_p[1,0]),label = 'L 2, P 1, E '+str(escape[1,0]))
#     plt.plot((x_p[1,1]),(y_p[1,1]),label = 'L 2, P 2, E '+str(escape[1,1]))
#     plt.plot((x_p[2,2]),(y_p[1,2]),label = 'L 2, P 3, E '+str(escape[1,2]))
#     plt.plot((x_p[1,3]),(y_p[1,3]),label = 'L 2, P 4, E '+str(escape[1,3]))
#     plt.plot((x_p[1,4]),(y_p[1,4]),label = 'L 2, P 5, E '+str(escape[1,4]))
#     plt.plot((x_p[1,5]),(y_p[1,5]),label = 'L 2, P 6, E '+str(escape[1,5]))
#     plt.plot((x_p[1,6]),(y_p[1,6]),label = 'L 2, P 7, E '+str(escape[1,6]))
#     plt.plot((x_p[1,7]),(y_p[1,7]),label = 'L 2, P 8, E '+str(escape[1,7]))

    # 3rd layer
#     plt.plot((x_p[2,0]),(y_p[2,0]),label = 'L 3, P 1, E '+str(escape[2,0]))
#     plt.plot((x_p[2,1]),(y_p[2,1]),label = 'L 3, P 2, E '+str(escape[2,1]))
#     plt.plot((x_p[2,2]),(y_p[2,2]),label = 'L 3, P 3, E '+str(escape[2,2]))
    # plt.plot((x_p[2,3]),(y_p[2,3]),label = 'L 3, P 4, E '+str(escape[2,3]))
#     plt.plot((x_p[2,4]),(y_p[2,4]),label = 'L 3, P 5, E '+str(escape[2,4]))
    # plt.plot((x_p[2,5]),(y_p[2,5]),label = 'L 3, P 6, E '+str(escape[2,5]))
    # plt.plot((x_p[2,6]),(y_p[2,6]),label = 'L 3, P 7, E '+str(escape[2,6]))
    # plt.plot((x_p[2,7]),(y_p[2,7]),label = 'L 3, P 8, E '+str(escape[2,7]))

    # 4th layer
#     plt.plot((x_p[3,0]),(y_p[3,0]),label = 'L 4, P 1, E '+str(escape[3,0]))
#     plt.plot((x_p[3,1]),(y_p[3,1]),label = 'L 4, P 2, E '+str(escape[3,1]))
#     plt.plot((x_p[3,2]),(y_p[3,2]),label = 'L 4, P 3, E '+str(escape[3,2]))
    # plt.plot((x_p[3,3]),(y_p[3,3]),label = 'L 4, P 4, E '+str(escape[3,3]))
#     plt.plot((x_p[3,4]),(y_p[3,4]),label = 'L 4, P 5, E '+str(escape[3,4]))
    # plt.plot((x_p[3,5]),(y_p[3,5]),label = 'L 4, P 6, E '+str(escape[3,5]))
    # plt.plot((x_p[3,6]),(y_p[3,6]),label = 'L 4, P 7, E '+str(escape[3,6]))
    # plt.plot((x_p[3,7]),(y_p[3,7]),label = 'L 4, P 8, E '+str(escape[3,7]))
    # 5th layer
#     plt.plot((x_p[4,0]),(y_p[4,0]),label = 'L 5, P 1, E '+str(escape[4,0]))
#     plt.plot((x_p[4,1]),(y_p[4,1]),label = 'L 5, P 2, E '+str(escape[4,1]))
#     plt.plot((x_p[4,2]),(y_p[4,2]),label = 'L 5, P 3, E '+str(escape[4,2]))
    # plt.plot((x_p[4,3]),(y_p[4,3]),label = 'L 5, P 4, E '+str(escape[4,3]))
#     plt.plot((x_p[4,4]),(y_p[4,4]),label = 'L 5, P 5, E '+str(escape[4,4]))
    # plt.plot((x_p[4,5]),(y_p[4,5]),label = 'L 5, P 6, E '+str(escape[4,5]))
    # plt.plot((x_p[4,6]),(y_p[4,6]),label = 'L 5, P 7, E '+str(escape[4,6]))
    # plt.plot((x_p[4,7]),(y_p[4,7]),label = 'L 5, P 8, E '+str(escape[4,7]))

    xla = plt.xlabel('X-coordinate of object')
    yla = plt.ylabel('Y-coordinate of object')
    end = plt.scatter(x[-1], y[-1], label = 'End', s = 100, color = 'yellow')

    bh                  = plt.Circle((0,0), r_ss, color = 'black', fill = True, lw = 2, label = 'Black hole')
    isco_radius         = plt.Circle((0,0), r_isco, color = 'purple', fill = False, lw = 2, label = 'Isco radius')
    tidal_radius        = plt.Circle((0,0), r_tidal, color = 'green', fill = False, lw = 2, label = 'Tidal radius')
    ax.add_patch(isco_radius)
    ax.add_patch(bh)
    ax.add_patch(tidal_radius)

    plt.xlim(-r_tidal*2, r_tidal*2)
    plt.ylim(-r_tidal*2, r_tidal*2)

    plt.title('Star motion in Black hole potential')
    plt.legend(loc = 'upper right')

    plt.show()
