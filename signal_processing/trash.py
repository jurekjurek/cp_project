
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






# for a Filter or length M:
# N = M+1
# x(n) = 1/M * X (exp())