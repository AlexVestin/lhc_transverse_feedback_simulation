import numpy as np
import scipy.signal

import matplotlib.pyplot as plt


def abs2(x):
    return x.real**2 + x.imag**2

if __name__ == '__main__':
    framelength=1.0
    N=1024
    x=np.linspace(0,framelength,N,endpoint=False)
    y=np.sin(44*2*np.pi*x)
    #y=y-np.mean(y)
    ffty = np.fft.fft(y)
    #power spectrum, after real2complex transfrom (factor )
    scale=2.0/(len(y)*len(y))
    power=scale*abs2(ffty)
    freq=np.fft.fftfreq(len(y) , framelength/len(y) )

    # power spectrum, via scipy welch. 'boxcar' means no window, nperseg=len(y) so that fft computed on the whole signal.
    freq2,power2=scipy.signal.welch(y, fs=len(y)/framelength,window='boxcar',nperseg=len(y),scaling='spectrum', axis=-1, average='mean')

    for i in range(len(freq2)):
        print (i, freq2[i], power2[i], freq[i], power[i])
    print (np.sum(power2))


    plt.figure()
    plt.plot(freq[0:int(len(y)/2+1)],power[0:int(len(y)/2+1)],label='np.fft.fft()')
    plt.plot(freq2,power2,label='scipy.signal.welch()')
    plt.legend()
    plt.xlim(0,np.max(freq[0:int(len(y)/2+1)]))


    plt.show()