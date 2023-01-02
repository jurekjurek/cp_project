import numpy as np
import matplotlib.pyplot as plt

# USE np.convolution!!!


f=200*10**6
T=100*10**(-6)
f_cutoff=48*10**6
n=int(T*f)

t=np.linspace(0,T,n)
signal=np.sin(2*np.pi*f*t)
for i in range(n):
    signal[i]+=np.random.normal(0,1)





# Generate random noise with a sampling frequency of 200 MHz and a total length of 100 microseconds.

def sample(sample_f,f, amplitude, t_max):
    sample_f = int(sample_f)
    x = np.linspace(0,t_max,sample_f)
    y_ = np.zeros(sample_f)
    for i in range(sample_f):
        y_[i] = f(x[i], amplitude)
    return x, y_

# generate random values between -1 and 1
def random_values(sampling_f, total_length):
    length = int(sampling_f*total_length)
    noise = np.random.random(length)*2-1
    return noise

noise = random_values(200*10**6, 0.1)
print(np.shape(noise))

# check if noise makes sense
# works on Macos, the plot does not seem to work on windows, because it's too large
# x = np.linspace(0, 1000, len(noise))

# plt.figure()
# plt.plot(noise[::5])
# plt.show()


# Implement a low-pass filter of length M=50 with a cutoff frequency of 48 MHz. Calculate the impulse response of the filter.
# For a low pass filter, I have to find a h(n) whose Fourier-trafo approximates this desired Freuquency response below(with a cutoff at 48 MHz)


# A will be a step function that only lets values smaller than 48MHz pass (amplitude function)
# 26 = 50/2 + 1 bc of np indexing
Amp = np.zeros(len(noise))
Amp = np.zeros(26)

interest = int((48)/(200) * len(Amp))


Amp[0:interest] = 1

print('#####')
print(Amp)

# check if the amplitude function looks right

# plt.figure()
# plt.plot(Amp)
# plt.title('Amplitude function; responsible for filtering')
# plt.show()


def h_of_n(M=50, A=Amp):
    N=M+1
    # the first element in h_sum will be zero, because in this case we only add A[0], which is 1
    h_sum = np.zeros(26, dtype='complex_')
    z = complex(0,1)
    print(z)
    h_array = np.zeros(26)
    for i in range(len(h_array)):         # compute ith element of h
        

        # if we vary the range of the loop we get completely different stuff
        for j in range(1,26):
            h_sum[i] += A[j]*np.exp(-z*np.pi*j*((M-2*i)/N))   # hier war vorher das j im exp nicht da 
    
    

    for i in range(len(h_array)):
        h_array[i] = (1/(N) * (A[0] + 2* h_sum[i]))
    print(h_array)
    return h_array

def h_of_n_whole(M=50):
    N=M+1
    A = np.zeros(N) 
    A[0:int((48)/(200) * len(Amp))]=1

    # the first element in h_sum will be zero, because in this case we only add A[0], which is 1
    h_sum = np.zeros(N)#, dtype='complex_')
    z = complex(0,1)
    print(z)
    h_array = np.zeros(N)
    for i in range(len(h_array)):         # compute ith element of h
        

        # if we vary the range of the loop we get completely different stuff
        for j in range(1,N):
            h_sum[i] += A[j]*np.exp(-z*np.pi*j*((M-2*i)/N)).real  # hier war vorher das j im exp nicht da 
    
    

    for i in range(len(h_array)):
        h_array[i] = (1/(N) * (A[0] + 2* h_sum[i]))
    print(h_array)
    return h_array

# Apply the impulse response to your generated noise.
def apply_to_noise():
    h = h_of_n_whole()
    plt.figure()
    plt.title('Left half of the Impulse response function')
    plt.xlabel('n')
    plt.ylabel('h(n)')
    plt.plot(h)
    plt.show()
    return np.convolve(h, noise)
    # new_array = np.zeros(len(noise))
    # number = len(noise) / len(h)
    # for i in range(len(h)):
    #     if 
    #     new_array[i] = 
    # return new_array

# transform into frequency domain
# A_back = np.fft.rfft(Amp)


convolution = apply_to_noise()

# plt.figure()
# plt.plot(np.fft.rfft(noise))
# plt.show()

# convolution = convolution[::4]

# Make a Fourier transform of your filtered noise and produce a plot of the frequency spectrum.
convolution_f = np.fft.rfft(convolution)

convolution_f = convolution_f[::4]
print(np.shape(convolution), np.shape(convolution_f))
plt.figure()
plt.plot(convolution_f)
# plt.xlim(0, 100*10**6)
plt.title('Frequency spectrum of the filtered noise')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.show()

# Make a Fourier transform of your filtered noise and produce a plot of the frequency spectrum.

# Make the same plot for a filter of length M=200 and discuss the results.