import numpy as np
import matplotlib.pyplot as plt

# function that returns convolution of 2 arrays
def convolution(x,y):
    convol_xy = np.zeros(len(x))
    for i in range(len(x)):
        for j in range(len(x)):
            if i <= j:
                convol_xy[i] += x[i]*y[j-i]
            # if (i+j)<len(convol_xy):
            #     convol_xy[i+j] = convol_xy[i+j]+x[i]*y[j]
    return convol_xy

x_input = np.array([1,-1,-0.5,0,0])
y_input = np.array([1,1,-0.5,0,0])
convolution_xy = convolution(x_input, y_input)
print(convolution(x_input, y_input))

y = np.linspace(0,4,5)

plt.figure()
plt.scatter(y, x_input, label = 'X')
plt.scatter(y, y_input, label = 'Y')
plt.plot(y, convolution_xy, label = 'Convolution')
plt.legend()
plt.show()