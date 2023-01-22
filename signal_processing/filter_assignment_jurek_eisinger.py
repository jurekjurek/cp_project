import numpy as np
import matplotlib.pyplot as plt

'''
In this last assignment we implement lowpassfilters of two different lengths, convolve the corresponding impulse response functions 
with random noise and fourier transform the signal to see the effects of the lowpass filter on the frequency distribution of the signal.
'''


def random_values(sampling_f, total_length):
    '''
    this function generates random values between -1 and 1 with a sampling frequency of sampling_f, 
    and a total length of total_length
    '''
    length = int(sampling_f*total_length)
    noise = np.random.random(length)*2-1        # generate random values between -1 and 1
    return noise

# generate random noise with given sampling frequency (200 MHz) and total length (0.1s)
noise = random_values(200*10**6, 0.1)


def impulse_response(M):
    '''
    this function calculates the impulse resonse for a given filter length M according to the formula derived in class
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
    return impulse_response/np.max(impulse_response)                    # return the normalized impulse response

def main(signal = noise):
    '''
    This function generates the plots of the impulse functions calculated above. 
    The fouriertransformation of the with the impulse response convolved signal is being carried out and plotted 
    for both filter lengths; 50 and 200. 
    We see that the filter of length 200 produces a frequency sprectrum that oscillates more rapidly, thus converging 
    (in the low-pass-filtered, so the high frequency region) faster to zero. Furthermore, the cutoff at the frequency of 48 MHz is 
    significantly sharper for the filter of length 200. 
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

    # convolve impulse response function for filter of length 50 with noise
    convolution_50 =  np.convolve(h_50, signal)
    # and fourier transfrom the convoluted signal to the frequency domain
    convolution_f_50 = np.fft.rfft(convolution_50)

    # convolve impulse response function for filter of length 200 with noise and fourier transform
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
For an ideal filter, we would expect a sharp cut at the cutoff frequency of 48MHz. 
For the filter of length 200, we see a more or less sharp cut at this cutoff frequency, and some further oscillation that decreases in amplitude. 
The filtered signal where we used the filter of length 50 shows a rather broad descent starting at around 44MHz. The following converging oscillation
has a significantly lower frequency compared to the signal filtered with a filter of length 200, thus converging slower and we see a higher amount of frequencies
larger than the cutoff frequency still being part of the signal. 
'''