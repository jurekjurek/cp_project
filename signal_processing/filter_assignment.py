import numpy as np
import matplotlib.pyplot as plt

# USE np.convolution!!!

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

# x = np.linspace(0,1000, len(noise))

# plt.figure()
# plt.plot(noise)
# plt.show()


# Implement a low-pass filter of length M=50 with a cutoff frequency of 48 MHz. Calculate the impulse response of the filter.

# this is the required fequency response
noise_f = np.fft.rfft(noise)
x = np.linspace(0,200*10**6, len(noise_f))
indices_ = np.where(x > 48*10**6)
# noise_f[indices_] = 0

# plt.figure()
# plt.plot(x, noise_f)
# plt.show()

noise_back = np.fft.irfft(noise_f)

# plt.figure()
# plt.plot(noise_back)
# plt.show()

# low pass filter:
# I have to find a h(n) whose Fourier-trafo approximates this desired Freuquency response above

# for a Filter or length M:
# N = M+1
# x(n) = 1/M * X (exp())

# A will be a step function, that only lets values smaller than 48MHz pass (amplitude function)
Amp = np.zeros(len(noise))
Amp = np.zeros(26)
# Amp[0:int(48*10**6*0.1)] = 1
Amp[0:int(len(Amp)/4)] = 1

def h_of_n(M=50, A=Amp):
    N=M+1
    h_sum = np.zeros(50, dtype='complex_')
    z = complex(0,1)
    print(z)
    h_array = np.zeros(50)
    for i in range(len(h_array)):         # compute ith element of h
        

        # if we vary the range of the loop we get completely different stuff
        for j in range(0,26):
            h_sum[i] += A[j]*np.exp(-z*np.pi*j*((M-2*i)/N))   # hier war vorher das j im exp nicht da 
    
    for i in range(len(h_array)):
        h_array[i] = (1/(N) * (A[0] + 2* h_sum[i]))
    return h_array




# Apply the impulse response to your generated noise.
def apply_to_noise():
    h = h_of_n()
    new_array = np.zeros(len(noise))
    for i in range(len(h)):
        new_array
    return new_array

# transform into frequency domain
# A_back = np.fft.rfft(Amp)


convolution = apply_to_noise()
convolution_f = np.fft.rfft(convolution)
plt.figure()
plt.plot(convolution_f)
plt.title('Frequency spectrum of the filtered noise')
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.show()

# Make a Fourier transform of your filtered noise and produce a plot of the frequency spectrum.

# Make the same plot for a filter of length M=200 and discuss the results.