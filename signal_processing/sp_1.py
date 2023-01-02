import numpy as np
import matplotlib.pyplot as plt

# 1st exercise 

c = 360/(2*np.pi)

def sine(x):
    return np.sin(10/0.1592*x)      # 10 Hz in rad 

def sample(sample_f,f, amplitude, t_max):
    sample_f = int(sample_f)
    x = np.linspace(0,t_max,sample_f)
    y_ = np.zeros(sample_f)
    for i in range(sample_f):
        y_[i] = f(x[i], amplitude)
    return x, y_

# x1 and x2 have same range, but different sampling
# we don't have to convert it from rad to Hz
# x_1, y_1 = sample(1000, sine)
# x_2, y_2 = sample(12, sine)
# x_3, y_3 = sample(41, sine)


# plt.figure()
# plt.plot(x_1,y_1, label = '1000 Hz')
# plt.scatter(x_1, y_1, label = '1000 Hz')
# plt.plot(x_2,y_2, label = '11 Hz')
# plt.scatter(x_2, y_2, label = '11 Hz')
# plt.plot(x_3,y_3, label = '21 Hz')
# plt.scatter(x_3, y_3, label = '21 Hz')
# plt.legend()
# plt.show()

# 2nd exercise

# from time to f domaine
# a = np.fft.rfft()

# n samples above, n/2 samples below, since complex -> more information

# from f to t
# b = np.fft.irfft()

# normalization 2pi * # of samples is normalization factor for digital fouoriertrafo

# what are the correct unit for the x-axis?
# what are correct frequencies?

def f(x, amplitude):
    return amplitude*np.sin(8/0.1592*x)+np.sin(12/0.1592*x)

# now sampling at 1000, 24, 20 Hz
# then (plot in f-domaine): plot amplitude of fft

x_1, y_1 = sample(1001, f, 1, 10)
x_2, y_2 = sample(24, f, 10, 10)
x_3, y_3 = sample(20, f, 3, 10)

# plt.figure()
# plt.plot(x_1,y_1, label = '1000 Hz')
# # plt.scatter(x_1, y_1, label = '1000 Hz')
# plt.plot(x_2,y_2, label = '25 Hz')
# plt.scatter(x_2, y_2, label = '25 Hz')
# plt.plot(x_3,y_3, label = '20 Hz')
# plt.scatter(x_3, y_3, label = '20 Hz')
# plt.legend()
# plt.show()

# print(len(x_1), len(x_2), len(x_3))

y_1_f = np.fft.rfft(y_1)
y_2_f = np.fft.rfft(y_2)
y_3_f = np.fft.rfft(y_3)



x_1_f = np.fft.rfftfreq(1)
x_2_f = np.fft.rfftfreq(1)
x_3_f = np.fft.rfftfreq(1)

print(len(y_1_f))

n_1 = len(y_1)
n_2 = len(y_2)
n_3 = len(y_3)

def asdf():
    if n_1 % 2 == 0:
        x_1_f = np.linspace(0,(1/2)+1,len(y_1_f))
    else:
        x_1_f = np.linspace(0,(1+1)/2,len(y_1_f))

    if n_2 % 2 == 0:
        x_2_f = np.linspace(0,1,len(y_2_f))
    else:
        x_2_f = np.linspace(0,1,len(y_2_f))

    if n_3 % 2 == 0:
        x_3_f = np.linspace(0,1,len(y_3_f))
    else:
        x_3_f = np.linspace(0,1,len(y_3_f))
    return x_1_f, x_2_f, x_3_f

# x_1_f, x_2_f, x_3_f = asdf()
x_1_f = np.linspace(0,1, len(y_1_f))
x_2_f = np.linspace(0,2, len(y_2_f))
x_3_f = np.linspace(0,2, len(y_3_f))

plt.figure()
plt.plot(x_1_f, y_1_f, label = '1000 Hz')
plt.plot(x_2_f, y_2_f, label = '25 Hz')
plt.plot(x_3_f, y_3_f, label = '20 Hz')
plt.ylabel('F(w)')
plt.xlabel('w')
plt.legend()
plt.show()

