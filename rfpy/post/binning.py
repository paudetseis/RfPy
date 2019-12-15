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
Functions to bin receiver functions along back-azimuth or slowness
bins, or both, or for all data, regardless of the x axis (time or depth).

"""

# Import modules and functions
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
from obspy.core import read, Stream, Trace, AttribDict
from matplotlib import pyplot as plt
from scipy.signal import hilbert


# Function to stack receiver functions into (baz or slow) bins
def bins(rfV, rfH, typ='baz', dd=36+1, pws=True):

    if not (typ == 'baz' or typ == 'slow' or typ == 'dist'):
        print('type has to be "baz" or "slow" or "dist"')
        return

    print()
    print('Stacking receiver functions by '+typ)

    # Define empty streams
    rfVbin = Stream()
    rfHbin = Stream()
    statsbin = rfV[0].stats

    if typ == 'baz':
        bmin = 0
        bmax = 360
        stat = [rfV[i].stats.sac.baz for i in range(len(rfV))]
    elif typ == 'slow':
        bmin = 0.04
        bmax = 0.08
        stat = [rfV[i].stats.sac.user0 for i in range(len(rfV))]
    elif typ == 'dist':
        bmin = 30.
        bmax = 90.
        stat = [rfV[i].stats.sac.user0 for i in range(len(rfV))]

    # Define bins
    bins = np.linspace(bmin, bmax, dd)

    # Digitize stat
    ind = np.digitize(stat, bins)

    # Loop through bins
    for i in range(len(bins)):

        nbin = 0
        Vtmp = np.zeros(len(rfV[0].data))
        Htmp = np.zeros(len(rfH[0].data))
        Vweight = np.zeros(len(rfV[0].data), dtype=complex)
        Hweight = np.zeros(len(rfH[0].data), dtype=complex)

        # Loop through stat
        for j in range(len(stat)):

            # If index of bins is equal to ind
            if i == ind[j]:

                nbin += 1
                Vtmp += rfV[j].data
                Htmp += rfH[j].data
                Vhilb = hilbert(rfV[j].data)
                Hhilb = hilbert(rfH[j].data)
                Vphase = np.arctan2(Vhilb.imag, Vhilb.real)
                Hphase = np.arctan2(Hhilb.imag, Hhilb.real)
                Vweight += np.exp(1j*Vphase)
                Hweight += np.exp(1j*Hphase)
                break

        if nbin > 0:

            # Average and update stats
            Vtmp = Vtmp/np.float(nbin)
            Htmp = Htmp/np.float(nbin)
            Vweight = Vweight/np.float(nbin)
            Hweight = Hweight/np.float(nbin)
            Vweight = np.real(abs(Vweight))
            Hweight = np.real(abs(Hweight))
            rfVt = Trace(header=statsbin)
            rfHt = Trace(header=statsbin)
            rfVt.stats.sac = AttribDict()
            rfHt.stats.sac = AttribDict()
            rfVt.stats.sac.user1 = nbin
            rfHt.stats.sac.user1 = nbin
            if typ == 'baz':
                rfVt.stats.sac.baz = bins[i]
                rfHt.stats.sac.baz = bins[i]
            elif typ == 'slow' or typ == 'dist':
                rfVt.stats.sac.user0 = bins[i]
                rfHt.stats.sac.user0 = bins[i]
            if pws is False:
                Vweight = np.ones(len(rfV[0].data))
                Hweight = np.ones(len(rfH[0].data))
            rfVt.data = Vweight*Vtmp
            rfHt.data = Hweight*Htmp
            rfVbin.append(rfVt)
            rfHbin.append(rfHt)

    return rfVbin, rfHbin


# Function to stack receiver functions into back-azimuth and slowness bins
def baz_slow_bins(rfV, rfH, dbaz=36+1, dslow=20+1, pws=True):

    print()
    print('Stacking receiver functions by baz and slow')

    # Define empty streams
    rfVbin = Stream()
    rfHbin = Stream()
    stats = rfV[0].stats
    stlo = stats.sac.stlo
    stla = stats.sac.stla

    # Define back-azimuth and slowness bins
    baz_bins = np.linspace(0, 360, dbaz)
    slow_bins = np.linspace(0.04, 0.08, dslow)

    # Extract baz and slowness
    baz = [rfV[i].stats.sac.baz for i in range(len(rfV))]
    slow = [rfV[i].stats.sac.user0 for i in range(len(rfV))]

    # Digitize baz and slowness
    ibaz = np.digitize(baz, baz_bins)
    islow = np.digitize(slow, slow_bins)

    # Loop through baz_bins
    for i in range(len(baz_bins)):
        for j in range(len(slow_bins)):

            nbin = 0
            Vtmp = np.zeros(len(rfV[0].data))
            Htmp = np.zeros(len(rfH[0].data))
            Vweight = np.zeros(len(rfV[0].data), dtype=complex)
            Hweight = np.zeros(len(rfH[0].data), dtype=complex)

            # Loop through baz
            for k in range(len(baz)):

                # If index of baz_bins is equal to ibaz
                if i == ibaz[k] and j == islow[k]:

                    nbin += 1
                    Vtmp += rfV[k].data
                    Htmp += rfH[k].data
                    Vhilb = hilbert(rfV[k].data)
                    Hhilb = hilbert(rfH[k].data)
                    Vphase = np.arctan2(Vhilb.imag, Vhilb.real)
                    Hphase = np.arctan2(Hhilb.imag, Hhilb.real)
                    Vweight += np.exp(1j*Vphase)
                    Hweight += np.exp(1j*Hphase)
                    break

            if nbin > 0:

                # Average and update stats
                Vtmp = Vtmp/np.float(nbin)
                Htmp = Htmp/np.float(nbin)
                Vweight = Vweight/np.float(nbin)
                Hweight = Hweight/np.float(nbin)
                Vweight = np.real(abs(Vweight))
                Hweight = np.real(abs(Hweight))
                rfVt = Trace(header=stats)
                rfHt = Trace(header=stats)
                rfVt.stats.sac = AttribDict()
                rfHt.stats.sac = AttribDict()
                rfVt.stats.sac.user1 = nbin
                rfHt.stats.sac.user1 = nbin
                rfVt.stats.sac.baz = baz_bins[i]
                rfHt.stats.sac.baz = baz_bins[i]
                rfVt.stats.sac.user0 = slow_bins[j]
                rfHt.stats.sac.user0 = slow_bins[j]
                rfVt.stats.sac.user1 = nbin
                rfHt.stats.sac.user1 = nbin
                rfVt.stats.sac.stlo = stlo
                rfHt.stats.sac.stlo = stlo
                rfVt.stats.sac.stla = stla
                rfHt.stats.sac.stla = stla

                if pws is False:
                    Vweight = np.ones(len(rfV[0].data))
                    Hweight = np.ones(len(rfH[0].data))
                rfVt.data = Vweight*Vtmp
                rfHt.data = Hweight*Htmp
                rfVbin.append(rfVt)
                rfHbin.append(rfHt)

    return rfVbin, rfHbin

# Stack all traces


def stack_all(rfV, rfH, pws=True):

    print()
    print('Stacking ALL receiver functions')

    # Copy stats from stream
    str_stats = rfV[0].stats

    # Initialize arrays
    Vtmp = np.zeros(len(rfV[0].data))
    Htmp = np.zeros(len(rfH[0].data))
    Vweight = np.zeros(len(rfV[0].data), dtype=complex)
    Hweight = np.zeros(len(rfH[0].data), dtype=complex)

    # Stack all traces
    for tr in rfV:
        Vtmp += tr.data
        Vhilb = hilbert(tr.data)
        Vphase = np.arctan2(Vhilb.imag, Vhilb.real)
        Vweight += np.exp(1j*Vphase)

    for tr in rfH:
        Htmp += tr.data
        Hhilb = hilbert(tr.data)
        Hphase = np.arctan2(Hhilb.imag, Hhilb.real)
        Hweight += np.exp(1j*Hphase)

    # Normalize
    Vtmp = Vtmp/np.float(len(rfV))
    Htmp = Htmp/np.float(len(rfH))
    Vweight = Vweight/np.float(len(rfV))
    Hweight = Hweight/np.float(len(rfH))
    Vweight = np.real(abs(Vweight))
    Hweight = np.real(abs(Hweight))

    # Fill up traces
    if pws is False:
        Vweight = np.ones(len(rfV[0].data))
        Hweight = np.ones(len(rfV[0].data))

    # Put back into traces
    rfVstack = Trace(data=Vweight*Vtmp, header=str_stats)
    rfHstack = Trace(data=Hweight*Htmp, header=str_stats)

    return rfVstack, rfHstack

# Migrate receiver functions


def migrate(rfV, rfH, nz=50, dep=None, vp=None):

    print()
    print('Migrating receiver functions to depth')

    # Initialize streams
    rfV_z = Stream()
    rfH_z = Stream()

    if (dep is None) and (vp is None):
        dep = np.array([0., 4., 8., 14., 25.9, 35.7, 45., 110.])
        vp = np.array([4.0, 5.9, 6.2, 6.3, 6.8, 7.2, 8.0, 8.1])

    # Get interpolated velocity profile
    idep = np.linspace(dep.min(), dep.max(), nz)
    ivp = sp.interpolate.interp1d(dep, vp, kind='linear')(idep)
    dz = idep[1]-idep[0]

    # Uniform Vp/Vs ratio
    R = 1.77

    # Calculate traveltime integral and stack
    # Loop over all traces
    for ib in range(len(rfV)):

        # Get slowness and back-azimuth
        slow = rfV[ib].stats.sac.user0
        baz = rfV[ib].stats.sac.baz
        stats = rfV[ib].stats

        # Initialize arrays
        Vtmp = np.zeros(nz)
        Htmp = np.zeros(nz)

        # Loop over interpolated depth range
        for iz in range(nz):

            # Initialize travel time
            tt = 0.

            # integrate travel times over depth
            for i in range(iz):
                tt = tt+dz*(np.sqrt((R/ivp[i])**2-slow**2) -
                            np.sqrt(1./ivp[i]**2-slow**2))

            # Shift traces to zero time
            Vtmp[iz] = timeshift(rfV[ib], tt)
            Htmp[iz] = timeshift(rfH[ib], tt)

        # Copy into new traces and update stats
        rfVt = Trace(data=Vtmp, header=stats)
        rfHt = Trace(data=Htmp, header=stats)
        rfVt.stats.delta = dz
        rfHt.stats.delta = dz
        rfVt.stats.npts = nz
        rfHt.stats.npts = nz

        rfV_z.append(rfVt)
        rfH_z.append(rfHt)

    return rfV_z, rfH_z


def timeshift(trace, tt):

    # Define frequencies
    nt = trace.stats.npts
    dt = trace.stats.delta
    freq = np.fft.fftfreq(nt, d=dt)

    # Fourier transform
    ftrace = np.fft.fft(trace.data)

    # Shift
    for i in range(len(freq)):
        ftrace[i] = ftrace[i]*np.exp(2.*np.pi*1j*freq[i]*tt)

    # Back Fourier transform
    rtrace = np.fft.ifft(ftrace)

    # Take first sample from trace
    val = np.real(rtrace[0])

    return val


def find_moho(rfV):

    nz = rfV.stats.npts
    dz = rfV.stats.delta
    depth = np.arange(nz)*dz

    tr = rfV.copy()
    tr.data[depth < 20.] = 0.
    tr.data[depth > 50.] = 0.
    moho = depth[np.where(tr.data == tr.data.max())]
    print('Moho is :', moho, 'km')

    return moho
