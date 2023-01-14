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
print(noise[-1])

# check if noise makes sense
# works on Macos, the plot does not seem to work on windows, because it's too large
# x = np.linspace(0, 1000, len(noise))

plt.figure()
plt.plot(noise)
plt.show()


# Implement a low-pass filter of length M=50 with a cutoff frequency of 48 MHz. Calculate the impulse response of the filter.
# For a low pass filter, I have to find a h(n) whose Fourier-trafo approximates this desired Freuquency response below(with a cutoff at 48 MHz)

def impulse_response(M):
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
        impulse_response[i] = (1/(N) * (A[0] + h_sum[i]))
        
    # print(impulse_response)
    return impulse_response/np.sum(impulse_response)

# Apply the impulse response to your generated noise.
def main(signal = noise):
    '''
    This function generates the plots of the impulse functions calculated above. 
    The fouriertransformation of the with the impulse response convolved signal is being carried out and plotted 
    for both filter lengths; 50 and 200. 
    We see that the filter of length 200 produces a frequency sprectrum that oscillates more rapidly, thus converging 
    (in the low-pass-filtered, so the high frequency region) faster to zero. 
    '''
    h_50 = impulse_response(M=50)
    h_200 = impulse_response(M=200)
    # plot the impulse response 
    plt.figure()
    plt.title('Impulse response function of a filter of length 50')
    plt.xlabel('n')
    plt.ylabel('h(n)')
    plt.plot(h_50, label = 'M = 50')
    plt.plot(h_200, label = 'M = 200')
    plt.legend()
    plt.show()
    
    # return None 

    # convolve impulse function for filter of length 50 with noise
    convolution_50 =  np.convolve(h_50, signal)
    convolution_f_50 = np.fft.rfft(convolution_50)

    # convolve impulse function for filter of length 200 with noise
    convolution_200 =  np.convolve(h_200, signal)
    convolution_f_200 = np.fft.rfft(convolution_200)

    # plot the noise in the frequency domain after applying the filter (convolving the impulse response with the noise and fourier transforming)
    plt.figure()
    plt.plot(convolution_f_50, label = 'M = 50', alpha = 0.6, color = 'green')
    plt.plot(convolution_f_200, label = 'M = 200', alpha = 0.6, color = 'blue')
    plt.title('Frequency spectrum of the filtered noise for two filter lenghts')
    plt.xlabel('Frequency')
    plt.ylabel('Amplitude')
    plt.legend() 
    plt.show()

    return None 

main()


'''
np.fft.rfftfreq() is the fuction you need für die Grenzen
'''





# Make a Fourier transform of your filtered noise and produce a plot of the frequency spectrum

# Make the same plot for a filter of length M=200 and discuss the results.