# Copyright 2019 Pascal Audet
#
# This file is part of RfPy.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
Functions to calculate piercing points from a velocity model and slowness
values, bin them, and produce CCP stacks using receiver functions.

"""

import numpy as np
import scipy as sp
from scipy.signal import hilbert

# Calculate horizontal distance for interval dz and velocity vs


def ppoint_distance(tr, dz, vs):

    # Get horizontal slowness
    slow = tr.stats.sac.user0

    # Calculate distance
    dx = dz*np.tan(np.arcsin(slow*vs))

    return dx

# Determine geographic location of piercing point


def ppoint(tr, dist):

    # Conversion factors
    lat2km = 111.
    lon2km = 90.

    # Get lat and lon of station location
    slat = tr.stats.sac.stla
    slon = tr.stats.sac.stlo

    # Back-azimuth of event
    baz = tr.stats.sac.baz*np.pi/180.

    # location of piercing point on geographical grid
    plat = dist*np.sin(-baz+np.pi/2.)/lat2km + slat
    plon = dist*np.cos(-baz+np.pi/2.)/lon2km + slon

    return plon, plat

# Calculate travel time for interval dz and velocities vp and vs


def ttime(tr, dz, vp, vs, phase='Ps'):

    # Get horizontal slowness
    slow = tr.stats.sac.user0

    # Calculate travel time for phase
    if phase == 'Ps':
        tt = dz*(np.sqrt((1./vs)**2 - slow**2) -
                 np.sqrt((1./vp)**2 - slow**2))
    elif phase == 'Pps':
        tt = dz*(np.sqrt((1./vs)**2 - slow**2) +
                 np.sqrt((1./vp)**2 - slow**2))
    elif phase == 'Pss':
        tt = 2.*dz*(np.sqrt((1./vs)**2 - slow**2))
    else:
        print('Error - unrecognized phase, ', phase)
        print('Returning tt = 0')
        tt = 0.

    return tt

# Shift a trace by a travel time tt and take amplitude at zero


def timeshift(tr, tt):

    # Define frequencies
    nt = int(tr.stats.npts)
    dt = tr.stats.delta
    freq = np.fft.fftfreq(int(nt), d=dt)
    trace = tr.data
    hilb = hilbert(tr.data)
    hilb_index = np.rint(tt/dt)
    hilb_tt = hilb[int(hilb_index)]
    hilb_tt_phase = np.arctan2(hilb_tt.imag, hilb_tt.real)

    # Fourier transform
    ftr = np.fft.fft(tr.data)

    # Shift using Fourier transform
    for i in range(len(freq)):
        # Fourier timeshift theorem
        ftr[i] = ftr[i]*np.exp(2.*np.pi*1j*freq[i]*tt)

    # Back to time domain (inverse Fourier transform)
    rtr = np.fft.ifft(ftr)

    # Take first sample from trace (convert to real value)
    amp = np.real(rtr[0])

    return amp, hilb_tt_phase


def raypath(tr, nz=50, dep=None, vp=None, vs=None):

    # Define arrays with zeros
    plat = np.zeros(nz)
    plon = np.zeros(nz)
    ttps = np.zeros(nz)
    ttpps = np.zeros(nz)
    ttpss = np.zeros(nz)

    # Default velocity model - can be updated later
    if (dep is None) and (vp is None) and (vs is None):
        dep = np.array([0., 4., 8., 14., 25.9, 35.7, 45., 110., 200.])
        vp = np.array([4.0, 5.9, 6.2, 6.3, 6.8, 7.2, 8.0, 8.1, 8.2])
        vs = vp/1.73

    # Get regular depth array
    idep = np.linspace(dep.min(), dep.max(), nz)

    # Interpolate Vp and Vs models on depth grid
    ivp = sp.interpolate.interp1d(dep, vp, kind='linear')(idep)
    ivs = sp.interpolate.interp1d(dep, vs, kind='linear')(idep)

    # Get exact depth interval
    dz = idep[1]-idep[0]

    # Now loop through all depths
    for iz in range(nz):

        # Initialize travel time and distance counters
        dtps = 0.
        dtpps = 0.
        dtpss = 0.
        dx = 0.

        # Sum over depths from 0 to iz
        for i in range(iz):
            dtps = dtps + ttime(tr, dz, ivp[i], ivs[i], 'Ps')
            dtpps = dtpps + ttime(tr, dz, ivp[i], ivs[i], 'Pps')
            dtpss = dtpss + ttime(tr, dz, ivp[i], ivs[i], 'Pss')
            dx = dx + ppoint_distance(tr, dz, ivs[i])

        # Get piercing point from distance
        plo, pla = ppoint(tr, dx)

        # Assign values to arrays
        ttps[iz] = dtps
        ttpps[iz] = dtpps
        ttpss[iz] = dtpss
        plon[iz] = plo
        plat[iz] = pla

    return ttps, ttpps, ttpss, plon, plat, idep
