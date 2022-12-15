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
v_x0 = np.sqrt(G*M/y0 + 3*G**2*M**2/(y0**2*c**2))*0.7

# escape velo, A = G*M
# we see that for 70% e.g. the particle falls back into the Mass
# v_x0 = 0
# v_y0 = np.sqrt(A*2/y0)*0.7

print('vy_0 is ', v_x0/c*100, '% of the speed of light.')

# def L(x,y, m):
#     '''
#     function that calculates the angular momentum for given x, y
#     L is constant over time
#     Since at point t = 0, we have only a velocity component in x-direction, so we calculate the L(0) = L(t) to be:
#     '''
#     x = float(x)
#     y = float(y)
#     r = np.sqrt((x**2+y**2))
#     return r*m*v_x0

# L = L(x0, y0, m)


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

##############################################################################################################
#                               Kepler orbit or BH orbit? 
kepler = 2              
# set to True for Newtonian potential, False for BH potential, which does not seem to work at all
# and to 2 for the alternative potential from the paper
##############################################################################################################






# escape velocity for newtonian potential
# v0 = np.sqrt(2/r0) 

# t = np.linspace(0,100000,10000)
t = np.linspace(0,1000000000,20000)
t0 = t[0]

h = len(t)


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




# the higher h, the more happens / the less detailed
# which makes sense
def my_rk4(x0, v_x0, y0, v_y0, star, h = 10):
    '''
    This function solves the differential equations (the equation of motion)
    for the x- and y-component of a particle. At each step, we evaluate v[i] for x and y which corresponds to their derivatives 
    as well as the derivatives dot(v) for x and y which corresponds to their second order derivatives. 
    Since for the (i+1)-th step, we need the values of x, y, v_x and v_y for the i-th step, for each step we calculate these 4 values
    star is a variable that tells us what kind of object we are solving the eom for, 
    star = True: The star that further splits into particles
    star = False: One particle the star turned into
    '''
    # F1,F2,F3,F4 = grav_layer(r_star, m)
    print('#')
    x = np.zeros(len(t))
    v_x = np.zeros(len(t))
    y = np.zeros(len(t))
    v_y = np.zeros(len(t))
    x[0] = x0
    v_x[0] = v_x0
    y[0] = y0
    v_y[0] = v_y0

    # lists to safe x,y,vx,vy for when the star crosses certain thresholds
    a_secure = 0            # to make sure of sth, gleich erklären

    # create arrays to store the position and velocity of the star for when each layer separates
    x_save =    np.zeros(5)
    y_save =    np.zeros(5)
    vx_save =   np.zeros(5) 
    vy_save =   np.zeros(5)

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

        # if the star (or the particles) passes the ss_radius, it's gone and we cannot describe its motion anymore
        if np.sqrt(x[i]**2+y[i]**2) <= r_ss:
            x[i] = 0
            y[i] = 0
            print('######### \n Star passed Schwarzschild radius! \n ########')
            # sets the valus after this time to zero, since it's not going to evolve anymore
            # x = x[:i]
            # y = y[:i]
            # v_x = v_x[:i]
            # v_y = v_y[:i]
            break

        # if we want the eom for the particles, we do not want to divide them up further into more particles
        # if we want the eom for the star, we will divide it into particles at this point
        if star == True: 
            

            # (a_secure makes sure that we only divide each layer into particles once and makes sure we can only split the i-1 th layer if we already split up the ith layer)
            # outermost (5th) shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*1.2 and a_secure == 0:
                x_save[0] = (x[i])
                y_save[0] = (y[i])
                vx_save[0] = (v_x[i])
                vy_save[0] = (v_y[i])
                print('a_secure = ', a_secure)
                print(x[i], y[i])
                a_secure = 1

            # 4th shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*1.1 and a_secure == 1:
                x_save[1] = (x[i])
                y_save[1] = (y[i])
                vx_save[1] = (v_x[i])
                vy_save[1] = (v_y[i])
                
                print('a_secure = ', a_secure)
                print(x[i], y[i])
                a_secure = 2

            # 3rd shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*1.00 and a_secure == 2:
                x_save[2] = (x[i])
                y_save[2] = (y[i])
                vx_save[2] = (v_x[i])
                vy_save[2] = (v_y[i])
                
                print('a_secure = ', a_secure)
                print(x[i], y[i])
                a_secure = 3

            # 2nd shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*0.975 and a_secure == 3:
                x_save[3] = (x[i])
                y_save[3] = (y[i])
                vx_save[3] = (v_x[i])
                vy_save[3] = (v_y[i])
                
                print('a_secure = ', a_secure)
                print(x[i], y[i])
                a_secure = 4

            # innermost (1st) shell:
            # if the innermost shell turns into particles, the star does not exist anymore as a star, so this runge kutta for-loop stops
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*0.95 and a_secure == 4:
                x_save[4] = (x[i])
                y_save[4] = (y[i])
                vx_save[4] = (v_x[i])
                vy_save[4] = (v_y[i]) 
                
                # the star is no longer flying at this point, so we shorten its array to the length it has now
                x[i] = 0
                y[i] = 0
                # sets the valus after this time to zero, since it's not going to evolve anymore
                x = x[:i]
                y = y[:i]
                v_x = v_x[:i]
                v_y = v_y[:i]

                print('a_secure = ', a_secure)  
                # print(x[i], y[i])
                a_secure = 5 
                break

     

    if star == True:
        if len(x_save) != 5 or len(y_save) != 5 or len(vx_save) != 5 or len(vy_save) != 5:
            print('Not every layer has turned into a star. Break.')
            # return x, y, v_x, v_y
        return x, y, v_x, v_y, x_save, y_save, vx_save, vy_save
    if star == False:
        return x, y, v_x, v_y



# creates 8 particles for each layer in the star 

def particles(r_layer, x_star, y_star, number_of_particles=8, scale = 20):
    '''
    This function creates initial positions for the particles once the shell splits up
    This initial position is dependent on the position of the star (x_star, y_star)
    All particles have the same velocity but different initial positions
    '''

    x_is = np.zeros(number_of_particles)
    y_is = np.zeros(number_of_particles)

    # for every one of the 8 particles in the j-th shell define the initial positions
    x_is[0] = 0
    y_is[0] = 1

    x_is[1] = 1/(np.sqrt(2))
    y_is[1] = 1/(np.sqrt(2))

    x_is[2] = 1
    y_is[2] = 0

    x_is[3] = 1/(np.sqrt(2))
    y_is[4] = -1/(np.sqrt(2))

    x_is[4] = 0
    y_is[4] = -1

    x_is[5] = -1/(np.sqrt(2))
    y_is[5] = -1/(np.sqrt(2))

    x_is[6] = -1
    y_is[6] = 0

    x_is[7] = -1/(np.sqrt(2))
    y_is[7] = 1/(np.sqrt(2))

    x_is = x_is * r_layer*scale
    y_is = y_is * r_layer*scale

    x_is = x_is + x_star    # adding the constant value x_star to the array x_is 
                            # every particle has an initial position relative to the star, then we add the star position to get the global initial position
    y_is = y_is + y_star    

    return x_is, y_is



# print('######### test particle function #########')
# x = [10,10,10,10,10]
# x_is, y_is = particles(1, 10, 10)
# print(x_is, y_is)

# print('######### test particle function #########')


def execute():

    # for now:
    x_star =  x0
    y_star =  y0 
    vx_star = v_x0
    vy_star = v_y0

    number_of_layers = 5
    number_of_particles = 8



    # we get x_star, y_star, vx_star ... from calling rk 4 for the star 
    # but! These will also be 1-dim np.arrays with length 5 (for every layer)
    # so: the velocity of the star when the outermost layer separated was: vx_star[0]
    # my_x, ... describes the actual movement of the star
    my_x, my_y, my_vx, my_vy, x_star, y_star, vx_star, vy_star = my_rk4(x0, v_x0,y0, v_y0, star=True)

    # Wir bekommen die richtigen x_save und y_save von rk4 wieder! 

    print('Das hier sind alle shapes... ', np.shape(x_star), np.shape(y_star), np.shape(vx_star), np.shape(vy_star))
    print('Und das hier sind die arrays an sich: ', x_star, y_star, vx_star, vy_star)
    print('####################################')

    '''
    Create the matrices with the initial conditions 
    '''
    # x_is and y_is will consist of the initial conditions (x and y) for every particle in every layer
    x_is = np.zeros((number_of_layers, number_of_particles))
    y_is = np.zeros((number_of_layers, number_of_particles))

    # if every star is divided in 5 equal radii, r_layer_list contains the radius for every shell, r[0] for 1st shell, r[1] for 2nd and so on...
    a = r_star/5
    # r_layer_list = np.array([a/2, 3/2*a, 5/2*a, 7/2*a, 9/2*a])
    # first element in every array considered corresponds to outermost layer
    r_layer_list = np.array([9/2*a, 7/2*a, 5/2*a, 3/2*a, a/2])
    for i in range(number_of_layers):
        # in every layer we have 8 particles created 
        # the positions where the particles are created are dependent on the radius of the shell and 
        # the position of the star at the point where the shell splits up
        x_is_, y_is_ = particles(r_layer_list[i], x_star[i], y_star[i])         # initial positions for particles in i-th layer 
        x_is[i,:] = x_is_
        y_is[i,:] = y_is_

    '''
    particles function works
    '''
    # print('#####')
    # print('x_is = ', x_is)
    # print('y_is = ', y_is)
    # print('#####')
    '''
    Eom for the particles are solved with the respective initial conditions 
    '''

    # here I solve the eom for the particles (for all 40)
    # create empty arrays for the particle eoms (for x and y only for now, we are only interested in the trajectory)
    x_final = np.zeros((5,8,len(t)))
    y_final = np.zeros((5,8,len(t)))
    vx_final = np.zeros((5,8,len(t)))
    vy_final = np.zeros((5,8,len(t)))
    
    print(np.shape(my_x))
    for i in range(number_of_layers):
        for j in range(number_of_particles):
            # solve the eom and save the trajectories x and y in x_final and y_final where x_is and y_is are the initial conditions for the [i,j]th particle
            x_final[i,j,:], y_final[i,j,:], vx_final[i,j,:], vy_final[i,j,:] = my_rk4(x_is[i,j], vx_star[i], y_is[i,j], vy_star[i], star=False)
            
    # for only one, test it: (0,0)th
    # print(x_is[0,0], y_is[0,0], x_star[0], y_star[0], vy_star[0], vx_star[0])
    # x_test, y_test, vx_test, vy_test = my_rk4(x_is[0,0], vx_star[0], y_is[0,0], vy_star[0], star=False)

    # print('initial conditions for rk that does not seem to work: ', x_is[0,0], y_is[0,0], vx_star[0], vy_star[0])

    # so now, we should have: 
    # 40 eoms for 40 particles, all movement in x_final and y_final
    # and 1 eom of the star. The star disappears at some point. 
    
    return my_x, my_y, my_vx, my_vy, x_final, y_final, vx_final, vy_final
    # return my_x, my_y, my_vx, my_vy, x_test, y_test, vx_test, vy_test

x, y, vx, vy, x_p, y_p, vx_p, vy_p = execute()
print(np.shape(x), np.shape(y), np.shape(vx), np.shape(vy), np.shape(x_p), np.shape(y_p), np.shape(vx_p), np.shape(vy_p))

'''
Saving files 
'''
np.save('x.npy', x)
np.save('y.npy', y)
np.save('vx.npy', vx)
np.save('vy.npy', vy)
np.save('x_p.npy', x_p)
np.save('y_p.npy', y_p)
np.save('vx_p.npy', vx_p)
np.save('vy_p.npy', vy_p)


'''
Plotting
'''







'''
Remember: x[0] is outermost, x[4] is innermost layer
'''

print('The tidal radius for our problem is: ', r_tidal)

fig, ax = plt.subplots()

plt.plot(x,y, label = 'Star', color = 'black')
plt.plot(x_p[4,0], y_p[4,0], label = 'innermost layer')
plt.plot(x_p[0,0], y_p[0,0], label = 'outermost layer')
plt.plot(x_p[1,0], y_p[1,0], label = '4th layer')
plt.plot(x_p[2,0], y_p[2,0], label = '3rd layer')
plt.plot(x_p[3,0], y_p[3,0], label = '2nd layer')
xla = plt.xlabel('X-coordinate of the particle')
yla = plt.ylabel('Y-coordinate of the particle')
# beg = plt.scatter(x[0], y[0], label = 'beginning')
# end = plt.scatter(x[-1], y[-1], label = 'end')          # last element of x and y
bh = plt.scatter(0,0, label= 'Black Hole', color = 'black')#, s = 300)
ss = plt.scatter(0,-r_ss, label = 'Schwarzschild radius', s = 5)
isco = plt.scatter(0,-r_isco, label = 'isco radius', s = 5)
# tidal = plt.scatter(0,-r_tidal, label = 'Tidal radius for star')
# if we maybe want to do it more fancy:
#
circle = plt.Circle((0,0), r_tidal, color = 'red', fill = False, lw = 2, label = 'Tidal radius')

ax.add_patch(circle)

plt.xlim(-y0*2, y0*2)
plt.ylim(-y0*2, y0*2)
plt.title('Motion of the particle in the Kepler field - solved with RK4')

plt.legend(loc = 'upper right')

plt.show()
