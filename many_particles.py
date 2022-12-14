import runge_kutta_star
from runge_kutta_star import grav_layer
'''
write a function that creates N particles out of i-th shell of star
'''


def particles(layer, x_star, y_star, vx_star, vy_star, number_of_particles=5):
    '''
    All particles have the same velocity but different initial positions
    '''
    vx0 = vx_star
    vy0 = vy_star

    # the radii of the different particles relative to the center of the star 
    r4,r3,r2,r1 = grav_layer(r_star, m_star)

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



