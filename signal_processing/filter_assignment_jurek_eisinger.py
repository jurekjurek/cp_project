import numpy as np
import matplotlib.pyplot as plt

f=200*10**6
T=100*10**(-6)
f_cutoff=48*10**6
n=int(T*f)

t=np.linspace(0,T,n)
signal=np.sin(2*np.pi*f*t)
for i in range(n):
    signal[i]+=np.random.normal(0,1)





# Generate random noise with a sampling frequency of 200 MHz and a total length of 100 microseconds.


def random_values(sampling_f, total_length):
    length = int(sampling_f*total_length)
    noise = np.random.random(length)*2-1        # generate random values between -1 and 1
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

def impulse_response(M=50):
    '''
    this function calculates the impulse resonse for a given filter length M accroding to the formula derived in class
    '''
    N=M+1
    # define the step function responsible for letting nothing over a given frequency pass
    A = np.zeros(N) 
    A[0:int((48)/(200) * len(A))]=1

    # the first element in h_sum will be zero, because in this case we only add A[0], which is 1
    h_sum = np.zeros(N)#, dtype='complex_')
    z = complex(0,1)
    impulse_response = np.zeros(N)
    for i in range(N):                                                  # compute ith element of h
        for j in range(1,N):
            h_sum[i] += A[j]*np.exp(-z*np.pi*j*((M-2*i)/N)).real        # ith element of h is sum over all 
        impulse_response[i] = (1/(N) * (A[0] + 3*h_sum[i]))
        
    # print(impulse_response)
    return impulse_response/np.sum(impulse_response)

# Apply the impulse response to your generated noise.
def main(signal = noise):
    h = impulse_response()

    # plot the impulse response 
    plt.figure()
    plt.title('Impulse response function of a filter of length 50')
    plt.xlabel('n')
    plt.ylabel('h(n)')
    plt.plot(h)
    plt.show()
    
    convolution =  np.convolve(h, signal)
    convolution_f = np.fft.rfft(convolution)

    # plot the noise in the frequency domain after applying the filter (convolving the impulse response with the noise and fourier transforming)
    plt.figure()
    plt.plot(convolution_f)
    plt.title('Frequency spectrum of the filtered noise')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.show()

    return None 

main()

# Make a Fourier transform of your filtered noise and produce a plot of the frequency spectrum

# Make the same plot for a filter of length M=200 and discuss the results.