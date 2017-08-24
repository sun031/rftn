#!/usr/bin/env python

"""
Utilities for RFTN code
"""

from scipy import signal
from obspy import read
import numpy as np
import math
import os
import matplotlib.pyplot as plt

def ricker(fp=1.0, dt=0.025, tlen=3, zerophase=True):
    """
    Generate waveform of a Ricker wavelet.


    Paramters
    =========
    fp : float
        peak frequency

    dt : float
        sampling interval

    tlen : float
        length of wavelet

    zerophase : bool
        the maximum of wavelet at zero if True

    Returns
    =======
    t : array_like
        time series
    rick : array_like
        Ricker wavelet

    Notes
    =====
    For details, please refer to  wiki.seg.org.

    """

    nt = np.int(tlen/dt+1)
    t = np.arange(nt)*dt
    pi2 = math.pi * math.pi
    fp2 = fp * fp

    if zerophase:
        t2 = t * t
    else:
        t2 = (t-1/fp) * (t-1/fp)

    x = pi2*fp2*t2

    print len(x), len(t)

    rick = (1 - 2.0*x) * np.exp(-x)

    return t, rick

def conv(sacfile, fp, savepath, zerophase=False):
    """
    Convolution between the reflectivity and the Ricker wavelet

    Parameters
    ==========
    sacfile : string
        filename of sac trace

    fp : float
        peak frequency of wavelet

    savepath : string
        path name to save the convoluted signal

    zerophase : bool
        the maximum of wavelet at zero if True


    Returns
    -------


    """
    try:
        os.makedirs(savepath)
    except:
        pass

    tr = read(sacfile)[0]
    dt = tr.stats.delta
    tlen = 4.0/fp
    s1 = tr.data
    t, s2 = ricker(fp=fp, dt=dt, tlen=tlen, zerophase=False)

    # c = signal.convolve(in1=s1, in2=s2, mode="full")
    c = signal.fftconvolve(in1=s1, in2=s2, mode="full")

    if not zerophase:
        troll = 1./fp
        ntroll = np.int(troll/dt)
        c = np.roll(c, -ntroll)

    bn = os.path.basename(sacfile)
    fn = "/".join([savepath, bn])
    tr.data = c
    tr.write(fn, format="SAC")


def autocorrelate(sacfile, oneside=True, savepath="None"):
    """
    Autocorrelation of a SAC trace


    Parameters
    ==========
    sacfile : string
        filename of sac trace

    oneside : bool
        output the right part of autocorrelation if True

    savepath : string
        path name to save the convoluted signal
        If `None`, save the same directory of sacfile.
        The suffix is ".ac".


    Returns
    -------


    """

    tr = read(sacfile)[0]
    sig = tr.data
    ac = signal.fftconvolve(sig, sig[::-1], mode="full")

    if oneside:
        h = np.int(len(ac)/2)
        ac = ac[h:]

    tr.data = ac
    tr.normalize()

    if savepath=="None":
        tr.write(filename=sacfile+".ac", format="SAC")
    else:
        bn = os.path.basename(sacfile)
        fn = "/".join([savepath, bn])
        try:
            os.makedirs(savepath)
        except:
            pass
        tr.write(filename=fn+".ac", format="SAC")

if __name__=="__main__":

    # t, w1 = ricker(fp=1.0, dt=0.025, tlen=3)
    # t, w2 = ricker(fp=1.0, dt=0.025, tlen=3, zerophase=False)
    # plt.plot(t, w1, 'r', t, w2)
    # plt.show()

    # sacfile = "s1.mod_sp.z"
    # conv(sacfile=sacfile, fp=1.0, savepath="conv")

    autocorrelate(sacfile="./conv/s1.mod_sp.z")