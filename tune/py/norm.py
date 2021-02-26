import numpy as np
import scipy.signal

import matplotlib.pyplot as plt
sr = 11245.

def fill_data(N):
    vals = []
    inc = 2 * np.pi * 83.21 # 83.21 magic freq
    for i in range(N):
        vals.append( np.sin( (i / sr) * inc ) )
    return np.asarray(vals)

def magnitudes(x, N):
    tau = np.pi * 2. / N
    mags = []

    for i in range(N):
        sum_real = 0
        sum_imag = 0

        for j in range(N):
            angle = tau * i * j
            sum_real += x[j] * np.cos( angle ) / N 
            sum_imag += x[j] * -np.sin( angle ) / N 
        
        mags.append( np.sqrt(sum_real**2 + sum_imag**2) )
    return np.asarray(mags)

if __name__ == '__main__':
    N = 1024*8
    m = 2
    # Signal with window size = 1024
    y1 = fill_data(N)#np.sin(44*2*np.pi*x1)
    ffty1 = np.fft.fft(y1)

    # signal with window size = 2048
    y2 = fill_data(N*m) #np.sin(44*2*np.pi*x2)
    ffty2 = np.fft.fft(y2)

    ffty1 = ffty1 / len(y1)
    ffty2 = (ffty2 / len(y2))[::m]

    # Magnitudes
    mag1 = np.abs(ffty1[0:int(len(ffty1)/2+1)])
    mag2 = np.abs(ffty2[0:int(len(ffty2)/2+1)])

    mag3 = magnitudes(y1, N)
    print(len(mag1), len(mag2), len(mag3), np.max(mag3))

    plt.figure()
    #plt.plot(mag2,label='N=2048 (downsampled)')
    plt.plot(mag3[0:int(len(mag3)/2+1)],label='dft')
    #plt.plot(mag1,label='N=1024')



    plt.legend()
    plt.show()