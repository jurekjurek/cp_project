import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation


# Earth model

def model_2BP(state, t):
    mu = 3.986004418E+05  # Earth's gravitational parameter  
                          # [km^3/s^2]
    x = state[0]
    y = state[1]
    # z = state[2]
    x_dot = state[2]
    y_dot = state[3]
    # z_dot = state[5]
    x_ddot = -mu * x / (x ** 2 + y ** 2) ** (3 / 2)
    y_ddot = -mu * y / (x ** 2 + y ** 2) ** (3 / 2)
    # z_ddot = -mu * z / (x ** 2 + y ** 2 + z ** 2) ** (3 / 2)
    dstate_dt = [x_dot, y_dot, x_ddot, y_ddot]
    return dstate_dt



# Initial Conditions
X_0 = -2500+10000  # [km]
Y_0 = -5500  # [km]
# Y_0 = 3400  # [km]
VX_0 = 7.5  # [km/s]
# VY_0 = 0.0  # [km/s]
VY_0 = 5.0  # [km/s]
state_0 = [X_0, Y_0, VX_0, VY_0]

# Time Array
t = np.linspace(0, 6*3600, 200)  # Simulates for a time period of 6
                                 # hours [s]


# Solving ODE
sol = odeint(model_2BP, state_0, t)
X_Sat = sol[:, 0]  # X-coord [km] of satellite over time interval 
Y_Sat = sol[:, 1]  # Y-coord [km] of satellite over time interval
# Z_Sat = sol[:, 2]  # Z-coord [km] of satellite over time interval


# Setting up Spherical Earth to Plot
N = 50
phi = np.linspace(0, 2 * np.pi, N)
theta = np.linspace(0, np.pi, N)
theta, phi = np.meshgrid(theta, phi)

r_Earth = 6378.14  # Average radius of Earth [km]
X_Earth = r_Earth * np.cos(phi) * np.sin(theta)
Y_Earth = r_Earth * np.sin(phi) * np.sin(theta)
# Z_Earth = r_Earth * np.cos(theta)






# Animation

dataSet = np.array([X_Sat, Y_Sat])
numDataPoints = len(t)


def animate_func(num):
    if np.sqrt(dataSet[0,num+1]**2+ dataSet[1,num+1]**2) <= r_Earth:
        line_ani.event_source.stop()
        return None
    # ax.clear()  # Clears the figure to update the line, point,   
                # title, and axes
    # ax.plot_su6rface(X_Earth, Y_Earth, Z_Earth, color='green', alpha=0.7)
    # Updating Trajectory Line (num+1 due to Python indexing)
    ax.plot(dataSet[0, :num+1], dataSet[1, :num+1], c='black')
    # Updating Point Location 
    ax.scatter(dataSet[0, num], dataSet[1, num], 
               c='black', marker='o', s = 1)

    # Setting Axes Limits
    ax.set_xlim([-50000, 50000])
    ax.set_ylim([-50000, 50000])
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
    # ax.set_zlabel('z')
    





# Plotting Earth and Orbit
fig = plt.figure()
ax = plt.axes()
# ax.plot3D(X_Sat, Y_Sat, Z_Sat, 'black')
# ax.view_init(30, 145)  # Changing viewing angle (adjust as needed)
# plt.title('Two-Body Orbit')
# ax.set_xlabel('X [km]')
# ax.set_ylabel('Y [km]')
# ax.set_zlabel('Z [km]')


# Make axes limits
# xyzlim = np.array([ax.get_xlim3d(), ax.get_ylim3d(),      
#                    ax.get_zlim3d()]).T
# XYZlim = np.asarray([min(xyzlim[0]), max(xyzlim[1])])
# ax.set_xlim3d(XYZlim)
# ax.set_ylim3d(XYZlim)
# ax.set_zlim3d(XYZlim * 3/4)

# Adding Constant Origin
# ax.plot3D(dataSet[0, 0], dataSet[1, 0], dataSet[2, 0],     
#                c='black', marker='o')

line_ani = animation.FuncAnimation(fig, animate_func, interval=100,   
                                   frames=numDataPoints)

# ax.plot_surface(X_Earth, Y_Earth, color='blue', alpha=0.4)
earth_plot = plt.Circle((0,0),  r_Earth, color = 'blue')
ax.add_patch(earth_plot)
plt.show()







