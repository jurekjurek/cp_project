import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('/Users/jurekeisinger/Documents/Uni/cp/signal_processing/findpulse.txt')


# due to f-trafo we get rid of the noise bc it's gaussian 
data_f = np.fft.rfft(data)

threshold = 100
pulse_ind = np.where(data_f > threshold)

# delete the frequencies of the two radio pulses
# thats all the pulses that lay above a certain threshold, here 100
data_f[pulse_ind] = 0

# do an inverse fourier-trafo to get in the actual time domain again 
data_without_pulse = np.fft.irfft(data_f)
pulse_of_interest = np.max(data_without_pulse)
# PLOT

# in this plot, the amplitudes without the two radio frequencies, we see the peak of the actual pulse now
plt.figure()
plt.plot(data_without_pulse)
# plt.hlines(pulse_of_interest,0,1,color = 'green', label = 'pulse')
# plt.plot(data_f)
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()
plt.show()
