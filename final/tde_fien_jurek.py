import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from scipy import constants

'''
This file simulates a Tidal disruption event as well as it is possible given the number of assumptions we have.
A star is split up into, in this case, 40 particles. The star is divided into 5 layers, each layer splits up into 8 particles. 
The layers split up when the star is closer a certain distance from the star. This distance is greatest for the outer layer and 
smallest for the inner layer.  
The interactions between particles are being neglected. 
The potential of the black hole is the actual true Potential, where we assume the angular momentum to be constant over time (for one particle). 
The angular momentum of the star is going to change when it splits up into particles, and these particles will have only a fraction of the initial angular momentum.
'''

# we use the true black hole potential here
kepler = False            


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

x0 = -1.5*r_tidal
y0 = r_tidal #* 3#3*10**9



def orbit_velo():
    '''
    returns the orbit velocities depending on potential
    '''
    if kepler == True:
        return np.sqrt(G*M/y0)
    if kepler == False:
        return np.sqrt(G*M/y0 + 3*G**2*M**2/(y0**2*c**2))

v_x0 = orbit_velo()*0.7
v_y0 = -0.2*v_x0


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

# a scales the Angular momentum, dependent on the mass
# A particle of mass 1/40 m will have a smaller angular momentum than the whole star

def v_x_(t,x,y,v_x,v_y, a):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -G*M*x / r**3 
    if kepler == False:
        return -A*x / r**3 + a**2*B*x / r**4 - a**2*C*x / r**5

def v_y_(t,x,y,v_x,v_y, a):
    r = np.sqrt((x**2+y**2))
    if kepler == True:
        return -G*M*y / r**3 
    if kepler == False:
        return -A*y / r**3 + a**2*B*y / r**4 - a**2*C*y / r**5

def x_(t,x,y,v_x,v_y):
    return v_x

def y_(t,x,y,v_x,v_y):
    return v_y



def runge_kutta(x0, v_x0, y0, v_y0, star, a, h = 20, r_threshold = 20*r_tidal, escape = False):
    '''
    star is a variable that tells us what kind of object we are solving the eom for, 
    star = True:    The star that further splits into particles
    star = False:   One particle the star turned into
    if a particle is farther away from the black hole than r_threshold, it disappears 
    escape will help us determine how many particles escaped the gravitational potential of the black hole
    '''
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
    a_secure = 0            # to make sure that the layers split in the right order and never in the same timestep

    # create arrays to store the position and velocity of the star for when each layer separates
    x_save =    np.zeros(5)
    y_save =    np.zeros(5)
    vx_save =   np.zeros(5) 
    vy_save =   np.zeros(5)

    for i in range(0,len(t)-1):
        k1x =   (x_(    t[i], x[i], y[i], v_x[i], v_y[i]))
        k1v_x = (v_x_(  t[i], x[i], y[i], v_x[i], v_y[i], a))
        k1y =   (y_(    t[i], x[i], y[i], v_x[i], v_y[i]))
        k1v_y = (v_y_(  t[i], x[i], y[i], v_x[i], v_y[i], a))

        k2x =   (x_(    t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2))
        k2v_x = (v_x_(  t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2, a))
        k2y =   (y_(    t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2))
        k2v_y = (v_y_(  t[i]+ h/2, x[i]+ h*k1x/2, y[i]+ h*k1y/2, v_x[i]+ h*k1v_x/2, v_y[i]+ h*k1v_y/2, a))

        k3x =   (x_(    t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2))
        k3v_x = (v_x_(  t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2, a))
        k3y =   (y_(    t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2))
        k3v_y = (v_y_(  t[i]+ h/2, x[i]+ h*k2x/2, y[i]+ h*k2y/2, v_x[i]+ h*k2v_x/2, v_y[i]+ h*k2v_y/2, a))
        
        k4x =   (x_(    t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2))
        k4v_x = (v_x_(  t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2, a))
        k4y =   (y_(    t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2))
        k4v_y = (v_y_(  t[i]+ h, x[i]+ h*k3x/2, y[i]+ h*k3y/2, v_x[i]+ h*k3v_x/2, v_y[i]+ h*k3v_y/2, a))


        # update the values of x,y,vx and vy
        x[i+1]  = x[i] + (k1x+2*k2x+2*k3x+k4x)*h/6
        v_x[i+1]  = v_x[i] + (k1v_x+2*k2v_x+2*k3v_x+k4v_x)*h/6
        y[i+1]  = y[i] + (k1y+2*k2y+2*k3y+k4y)*h/6
        v_y[i+1]  = v_y[i] + (k1v_y+2*k2v_y+2*k3v_y+k4v_y)*h/6


        if np.sqrt(x[i]**2+y[i]**2) > r_threshold:
            escape = True
            x[i:] = x[i]
            y[i:] = y[i]
            print('The particle escapes!')
            break 

        # if the star (or the particles) passes the ss_radius, it's gone and we cannot describe its motion anymore
        if np.sqrt(x[i]**2+y[i]**2) <= r_ss:
            x[i] = 0
            y[i] = 0
            print('######### \n Object passed Schwarzschild radius! \n ########')
            break

        # if we want the eom for the particles, we do not want to divide them up further into more particles
        # if we want the eom for the star, we will divide it into particles at this point
        if star == True: 
            

            # (a_secure makes sure that we only divide each layer into particles once and makes sure we can only split the i-1 th layer if we already split up the ith layer)
            # outermost (5th) shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*1.1 and a_secure == 0:
                x_save[0] = (x[i])
                y_save[0] = (y[i])
                vx_save[0] = (v_x[i])
                vy_save[0] = (v_y[i])
                print('Outermost layer splits')
                a = 4/5                             # the star loses 1/5 of his mass
                a_secure = 1                        # update a so that next layer can split

            # 4th shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*1.05 and a_secure == 1:
                x_save[1] = (x[i])
                y_save[1] = (y[i])
                vx_save[1] = (v_x[i])
                vy_save[1] = (v_y[i])
                print('4th layer splits')
                a = 3/5                             # the star loses 1/5 of his mass
                a_secure = 2                        # update a so that next layer can split

            # 3rd shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*1.00 and a_secure == 2:
                x_save[2] = (x[i])
                y_save[2] = (y[i])
                vx_save[2] = (v_x[i])
                vy_save[2] = (v_y[i])
                print('3rd layer splits')
                a = 2/5                             # the star loses 1/5 of his mass
                a_secure = 3                        # update a so that next layer can split

            # 2nd shell:
            if np.sqrt(x[i]**2+y[i]**2) <= r_tidal*0.975 and a_secure == 3:
                x_save[3] = (x[i])
                y_save[3] = (y[i])
                vx_save[3] = (v_x[i])
                vy_save[3] = (v_y[i])
                print('2nd layer splits')
                a = 1/5                             # the star loses 1/5 of his mass
                a_secure = 4                        # update a so that next layer can split

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
                print('Innermost layer splits')
                a_secure = 5 
                break

     

    if star == True:
        if len(x_save) != 5 or len(y_save) != 5 or len(vx_save) != 5 or len(vy_save) != 5:
            print('Not every layer has turned into a star.')
        return x, y, v_x, v_y, x_save, y_save, vx_save, vy_save
    if star == False:
        return x, y, v_x, v_y, escape



# creates 8 particles for each layer in the star 

def particles(r_layer, x_star, y_star, number_of_particles=8, scale = 10):
    '''
    This function creates initial positions for the particles once the shell splits up
    This initial position is dependent on the position of the star (x_star, y_star) but not on its velocity
    All particles have the same velocity but different initial positions
    The variable scale indicates how far away from the center of the scale the particle will be, scale=1 corresponds
    to the radius of the corresponding shell, shell = 2 corresponds to two times this radius... 
    '''

    x_is = np.zeros(number_of_particles)
    y_is = np.zeros(number_of_particles)

    # for every one of the 8 particles in the j-th shell define the initial positions relative to the star
    # the positions are in a circular shape around the star
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

    x_is = x_is * r_layer*scale                 # distance from star of each particle is set to the distance of its shell *plus*
    y_is = y_is * r_layer*scale                 # a scale due to the 'explosion' of the star

    x_is = x_is + x_star    # adding the constant value x_star to the array x_is 
                            # every particle has an initial position relative to the star, then we add the star position to get the global initial position
    y_is = y_is + y_star    

    return x_is, y_is




def execute():

    # for now:
    x_star =  x0
    y_star =  y0 
    vx_star = v_x0
    vy_star = v_y0

    number_of_layers = 5
    number_of_particles = 8

    x_star, y_star, vx_star, vy_star, x_initial, y_initial, vx_initial, vy_initial = runge_kutta(x0, v_x0,y0, v_y0, star=True, a=1)

    
    # Create the matrices with the initial conditions 
    

    # x_is and y_is will consist of the initial conditions (x and y) for every particle in every layer
    x_is = np.zeros((number_of_layers, number_of_particles))
    y_is = np.zeros((number_of_layers, number_of_particles))

    # if every star is divided in 5 equal radii, r_layer_list contains the radius for every shell, r[0] for 1st shell, r[1] for 2nd and so on...
    a = r_star/5
    # first element in every array corresponds to outermost layer
    r_layer_list = np.array([9/2*a, 7/2*a, 5/2*a, 3/2*a, a/2])
    for i in range(number_of_layers):
        # in every layer we have 8 particles created 
        # the positions where the particles are created are dependent on the radius of the shell and 
        # the position of the star at the point where the shell splits up
        x_is_, y_is_ = particles(r_layer_list[i], x_initial[i], y_initial[i])         # initial positions for particles in i-th layer 
        x_is[i,:] = x_is_
        y_is[i,:] = y_is_


    '''
    Eom for the particles are solved with the respective initial conditions 
    '''

    # here I solve the eom for the particles (for all 40)
    # create empty arrays for the particle eoms (for x and y only for now, we are only interested in the trajectory)
    x_final = np.zeros((5,8,len(t)))
    y_final = np.zeros((5,8,len(t)))
    vx_final = np.zeros((5,8,len(t)))
    vy_final = np.zeros((5,8,len(t)))
    escape = np.zeros((5,8))
    
    for i in range(number_of_layers):
        for j in range(number_of_particles):
            # solve the eom and save the trajectories x and y in x_final and y_final where x_is and y_is are the initial conditions for the [i,j]th particle
            # we set a = 1/40, since every particle has only 1/40 of the initial angular momentum of the star
            x_final[i,j,:], y_final[i,j,:], vx_final[i,j,:], vy_final[i,j,:], escape[i,j] = runge_kutta(x_is[i,j], vx_initial[i], y_is[i,j], vy_initial[i], star=False, a=1/40)

    
    return x_star, y_star, vx_star, vy_star, x_final, y_final, vx_final, vy_final, escape


def main():
    '''
    This function executes the program and creates the necessary plots
    '''
    x, y, vx, vy, x_p, y_p, vx_p, vy_p, escape = execute()

    # This many particles escaped:
    print('This many particles escaped:' ,np.sum(escape), 'out of 40.')

    '''
    Saving numpy arrays to plot in program plot_fien_jurek.py
    '''
    np.save('x.npy', x)
    np.save('y.npy', y)
    np.save('vx.npy', vx)
    np.save('vy.npy', vy)
    np.save('x_p.npy', x_p)
    np.save('y_p.npy', y_p)
    np.save('vx_p.npy', vx_p)
    np.save('vy_p.npy', vy_p)
    np.save('escape.npy', escape)


    '''
    Plotting
    '''
    # x[0] is outermost, x[4] is innermost layer


    print('The tidal radius for our problem is: ', r_tidal)

    fig, ax = plt.subplots()

    plt.plot(x,y, label = 'Star', color = 'black')
    plt.plot(x_p[4,0], y_p[4,0], label = 'innermost layer')
    plt.plot(x_p[0,0], y_p[0,0], label = 'outermost layer')
    plt.plot(x_p[1,0], y_p[1,0], label = '4th layer')
    plt.plot(x_p[2,0], y_p[2,0], label = '3rd layer')
    plt.plot(x_p[3,0], y_p[3,0], label = '2nd layer')
    xla = plt.xlabel('X-coordinate of object')
    yla = plt.ylabel('Y-coordinate of object')
    # end = plt.scatter(x[-1], y[-1], label = 'End', s = 100, color = 'yellow')

    bh = plt.Circle((0,0), r_ss, color = 'black', fill = True, lw = 2, label = 'Black hole')
    isco_radius = plt.Circle((0,0), r_isco, color = 'purple', fill = False, lw = 2, label = 'Isco radius')
    tidal_radius = plt.Circle((0,0), r_tidal, color = 'green', fill = False, lw = 2, label = 'Tidal radius')
    ax.add_patch(isco_radius)
    ax.add_patch(bh)
    ax.add_patch(tidal_radius)

    plt.xlim(-y0*4, y0*4)
    plt.ylim(-y0*4, y0*4)

    plt.title('Star motion in Black hole potential')
    plt.legend(loc = 'upper right')

    plt.show()

main()
