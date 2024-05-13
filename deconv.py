import numpy as np
from scipy.signal import welch, convolve
from scipy.signal import convolve, butter, filtfilt
import pandas as pd
from numpy.fft import fft, ifft
from matplotlib import pyplot as plt

###### User-defined functions
# functions for gACh and rACh impulse responses in time domain 
def gach_impulse_response(time):
    falltime = 0.580
    return 1/(falltime)*np.exp(-(time-np.min(time))/falltime) 

def rach_impulse_response(time):
    falltime = 0.245
    return 1/(falltime)*np.exp(-(time-np.min(time))/falltime) 





###### Helper functions
# convolves two signals and trims the ends to remove edge effects
def trim_convolve(x, y, time, trim_samples):
    start = trim_samples
    end = len(time) - trim_samples
    return convolve(x, y)[start:end]

def trim(x, trim_samples = 100):
    return x[trim_samples:-trim_samples]

# trimmed reconstruction of signal given deconvolved activity, time, and impulse response
def trim_reconst(deconv, time, impulse_response = gach_impulse_response):
    return trim_convolve(deconv, impulse_response(time), time, trim_samples = 100) 

def abs_err(x, y):
    n = len(x)
    return np.sum(np.abs(x-y))/n

# get coefficients for default 1Hz lowpass filter
def get_lp(fs, cutoff = 1, order=3):
    return butter(order, cutoff, fs=fs, btype='lowpass', analog=False)

# filter data
def lp(data, fs, cutoff = 1, order=3):
    b, a = get_lp(fs, cutoff = cutoff, order=order)
    #print(data)
    y = filtfilt(b, a, np.array(data))
    return y

# normalize signal to minimum 0; maximum 1
def mmnorm(x): 
    return (x - np.min(x))/(np.max(x) - np.min(x))

# get array element nearest to value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



##### Main functions

def wiener_deconv_noncausal(x, time, W = 0.1, impulse_response = gach_impulse_response):
    """ Deconvolves signal x, with timepoints time, and parameter W
    impulse_response is passed as a function of time (defaults to gACh impulse response from Neyhart et al)
    Returns deconvolved trace, filter H, and parameter W
    """
    k = impulse_response(time)
    K = fft(k)
    S = fft(x)
    H = np.conj(K)/(K*np.conj(K) + W)
    deconvolved = ifft(S*H)
    return deconvolved, H, W

def param_search(x, time, fs, lowpass = False, tol = 0.1, impulse_response = gach_impulse_response, params = np.linspace(0, 20, 200)):
    """ Returns Wiener deconvolution parameter best matching an error tolerance. 
    Uses trimmed signal and reconstruction (100 samples from both ends) to compute error.
    (default case: tol=0.1 means the filter will be chosen at 10% of the "worst-case" reconstruction error --- 
    this assumes the computed errors approximately cover the range of possible reconstruction errors
    which will be the case if min(params) = 0 and max(params) is large enough)
    """
    err = []
    if lowpass:
        x = lp(x, fs)
    stdx = (x - np.mean(x))/np.std(x)
    for W in params:
        
        deconv, _, _ = wiener_deconv_noncausal(stdx, time, W, impulse_response)
        trim_rec = trim_reconst(deconv, time, impulse_response)
        # only want to compute error on trimmed trace/trimmed reconstruction
        err.append(abs_err(trim(stdx, 100), trim_rec))
    err = mmnorm(np.real(err))
    idx = find_nearest(err, tol)
    return params[idx], 



def deconv_reconst(trace, time, impulse_response, W):
    """ Returns deconvolved trace, reconstructed trace, filter, and filter parameter.
    Impulse response is passed as a function of time.
    Trimming edge effects is left to the user here.
    """
    print(f'Deconvolving (noncausal method; |Wn| = {W})...')
    deconvolved, H, W = wiener_deconv_noncausal(trace, time, W, impulse_response)
    reconst = convolve(impulse_response(time), deconvolved)[0:len(time)]
    return deconvolved, reconst, H, W