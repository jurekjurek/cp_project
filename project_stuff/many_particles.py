import runge_kutta_star
from runge_kutta_star import grav_layer
'''
write a function that creates N particles out of i-th shell of star
'''


def particles(r_layer, x_star, y_star, vx_star, vy_star, number_of_particles=8):
    '''
    All particles have the same velocity but different initial positions
    '''
    vx0 = vx_star
    vy0 = vy_star

    # the radii of the different particles relative to the center of the star 
    r4,r3,r2,r1 = grav_layer(r_star, m_star)

    x_is = np.zeros(number_of_particles)
    y_is = np.zeros(number_of_particles)

    # for every one of the 9 particles in the j-th shell define the initial positions
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

    x_is = x_is * r_layer
    y_is = y_is * r_layer

    for i in range(number_of_particles):
        x_is[i] = x_star + x_is[i]
        y_is[i] = y_star + y_is[i]
        

    # for i in range(number_of_particles):



