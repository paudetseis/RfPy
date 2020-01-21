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

import os
import sys
import pickle
import numpy as np
import scipy as sp
from scipy.signal import hilbert
from rfpy import binning
import matplotlib.pyplot as plt
from matplotlib import cm


class CCPimage(object):

    def __init__(self, coord_start=[None, None], coord_end=[None, None],
                 weights=[1., 1., -1.],
                 dep=np.array([0., 4., 8., 14., 30., 35., 45., 110.]),
                 vp=np.array([4.0, 5.9, 6.2, 6.3, 6.8, 7.2, 8.0, 8.1]),
                 vpvs=1.73):

        self.radialRF = []
        self.dep = dep
        self.vp = vp
        self.vs = vp/vpvs
        self.weights = weights
        self.xs_lat1 = coord_start[0]
        self.xs_lon1 = coord_start[1]
        self.xs_lat2 = coord_end[0]
        self.xs_lon2 = coord_end[1]
        self.is_ready_for_prep = False
        self.is_ready_for_prestack = False
        self.is_ready_for_ccp = False
        self.is_ready_for_gccp = False

    def add_rfstream(self, rfstream):

        if len(rfstream)>0:
            self.radialRF.append(rfstream)
            self.is_ready_for_prep = True

    def prep_data(self, f1=0.05, f2ps=0.5, f2pps=0.3, f2pss=0.3, n_depth=120,
                  nbaz=36+1, nslow=40+1, save=None):

        if not self.is_ready_for_prep:
            raise(Exception("CCPimage not ready for pre-prep"))

        i_key = 0

        # Process streams one at a time
        for RF in self.radialRF:

            # Bin RFs into back-azimuth and slowness bins to speed up
            # calculations
            RFbin = binning.bin_baz_slow(
                RF, nbaz=nbaz, nslow=nslow, pws=True)[0]
            n_traces = len(RFbin)
            amp_ps_tr = np.empty([n_traces, n_depth])
            amp_pps_tr = np.empty([n_traces, n_depth])
            amp_pss_tr = np.empty([n_traces, n_depth])
            lon_tr = np.empty([n_traces, n_depth])
            lat_tr = np.empty([n_traces, n_depth])

            st_ps = RFbin.copy()
            st_pps = RFbin.copy()
            st_pss = RFbin.copy()

            # Filter Ps, Pps and Pss
            st_ps.filter(
                'bandpass', freqmin=f1, freqmax=f2ps,
                corners=4, zerophase=True)
            st_pps.filter(
                'bandpass', freqmin=f1, freqmax=f2pps,
                corners=4, zerophase=True)
            st_pss.filter(
                'bandpass', freqmin=f1, freqmax=f2pss,
                corners=4, zerophase=True)
            del RFbin

            print("Station: "+st_ps[0].stats.station)
            for itr in _progressbar(range(len(st_ps)), '', 25):

                # Get raypath and travel time for all phases
                tt_ps, tt_pps, tt_pss, plon, plat, idep = \
                    raypath(st_ps[itr], nz=n_depth,
                            dep=self.dep, vp=self.vp, vs=self.vs)

                # Now get amplitude of RF at corresponding travel
                # time along the raypath
                depth_array = np.asarray(idep)
                lon_tr[itr, :] = plon
                lat_tr[itr, :] = plat

                amp_ps = []
                amp_pps = []
                amp_pss = []

                # Loop through travel times and shift RFs to get amplitudes
                for tt in tt_ps:
                    a, phase = timeshift(st_ps[itr], tt)
                    amp_ps.append(self.weights[0]*a)
                amp_ps_tr[itr, :] = amp_ps

                # Loop through travel times and shift RFs to get amplitudes
                for tt in tt_pps:
                    a, phase = timeshift(st_pps[itr], tt)
                    amp_pps.append(self.weights[1]*a)
                amp_pps_tr[itr, :] = amp_pps

                # Loop through travel times and shift RFs to get amplitudes
                for tt in tt_pss:
                    a, phase = timeshift(st_pss[itr], tt)
                    amp_pss.append(self.weights[2]*a)
                amp_pss_tr[itr, :] = amp_pss

            if i_key == 0:
                amp_ps_depth = amp_ps_tr.transpose()
                amp_pps_depth = amp_pps_tr.transpose()
                amp_pss_depth = amp_pss_tr.transpose()
                lon_depth = lon_tr.transpose()
                lat_depth = lat_tr.transpose()

            elif i_key > 0:
                amp_ps_depth = np.concatenate(
                    (amp_ps_depth, amp_ps_tr.transpose()), axis=1)
                amp_pps_depth = np.concatenate(
                    (amp_pps_depth, amp_pps_tr.transpose()), axis=1)
                amp_pss_depth = np.concatenate(
                    (amp_pss_depth, amp_pss_tr.transpose()), axis=1)
                lon_depth = np.concatenate(
                    (lon_depth, lon_tr.transpose()), axis=1)
                lat_depth = np.concatenate(
                    (lat_depth, lat_tr.transpose()), axis=1)

            i_key += 1

        self.amp_ps_depth = amp_ps_depth
        self.amp_pps_depth = amp_pps_depth
        self.amp_pss_depth = amp_pss_depth
        self.lon_depth = lon_depth
        self.lat_depth = lat_depth
        self.depth_array = depth_array
        self.is_ready_for_prestack = True

        del self.radialRF

    def prestack(self, cell_length=1.):

        if not self.is_ready_for_prestack:
            raise(Exception("CCPimage not ready for prestack"))

        (n_depth, n_traces) = self.lon_depth.shape

        # Get total length of grid from end points
        xs_length = haversine(self.xs_lat1, self.xs_lon1,
                              self.xs_lat2, self.xs_lon2)

        # number of cells laterally for specified cell_length (rounded)
        n_lateral = int(np.rint(xs_length/cell_length))

        xs_latitudes = np.asarray(
            np.linspace(self.xs_lat1, self.xs_lat2, n_lateral))
        xs_longitudes = np.asarray(
            np.linspace(self.xs_lon1, self.xs_lon2, n_lateral))
        lateral_distances = np.arange(n_lateral)*cell_length

        xs_amps_ps = np.zeros((n_depth, n_lateral, n_traces))
        xs_amps_pps = np.zeros((n_depth, n_lateral, n_traces))
        xs_amps_pss = np.zeros((n_depth, n_lateral, n_traces))

        for i_depth in _progressbar(range(n_depth), '', 25):

            for i_coor in range(n_traces):

                lat_tr = self.lat_depth[i_depth, i_coor]
                lon_tr = self.lon_depth[i_depth, i_coor]
                distance_tests = np.empty(n_lateral)

                for i_xs in range(n_lateral):
                    lat_xs = xs_latitudes[i_xs]
                    lon_xs = xs_longitudes[i_xs]
                    distance_tests[i_xs] = haversine(
                        lat_xs, lon_xs, lat_tr, lon_tr)

                minimum_distance = np.amin(distance_tests)
                i_cell = np.where(distance_tests ==
                                  np.amin(distance_tests))[0][0]

                nonzero_count = np.count_nonzero(
                    xs_amps_ps[i_depth, i_cell, :])
                new_amp_ps = self.amp_ps_depth[i_depth, i_coor]
                if xs_amps_ps[i_depth, i_cell, 0] == 0.:
                    xs_amps_ps[i_depth, i_cell, 0] = new_amp_ps
                else:
                    xs_amps_ps[i_depth, i_cell, nonzero_count] = new_amp_ps

                nonzero_count = np.count_nonzero(
                    xs_amps_pps[i_depth, i_cell, :])
                new_amp_pps = self.amp_pps_depth[i_depth, i_coor]
                if xs_amps_pps[i_depth, i_cell, 0] == 0.:
                    xs_amps_pps[i_depth, i_cell, 0] = new_amp_pps
                else:
                    xs_amps_pps[i_depth, i_cell, nonzero_count] = new_amp_pps

                nonzero_count = np.count_nonzero(
                    xs_amps_pss[i_depth, i_cell, :])
                new_amp_pss = self.amp_pss_depth[i_depth, i_coor]
                if xs_amps_pss[i_depth, i_cell, 0] == 0.:
                    xs_amps_pss[i_depth, i_cell, 0] = new_amp_pss
                else:
                    xs_amps_pss[i_depth, i_cell, nonzero_count] = new_amp_pss

        self.xs_amps_ps = xs_amps_ps
        self.xs_amps_pps = xs_amps_pps
        self.xs_amps_pss = xs_amps_pss
        self.lateral_distances = lateral_distances
        self.n_lateral = n_lateral
        self.n_depth = n_depth
        self.is_ready_for_ccp = True
        self.is_ready_for_gccp = True

    def ccp(self):

        if not self.is_ready_for_ccp:
            raise(Exception("CCPimage not ready for ccp"))

        xs_ps_avg = np.zeros((self.n_depth, self.n_lateral))
        xs_pps_avg = np.zeros((self.n_depth, self.n_lateral))
        xs_pss_avg = np.zeros((self.n_depth, self.n_lateral))

        for i_depth in _progressbar(range(self.n_depth), '', 25):

            for i_cell in range(self.n_lateral):

                nonzero_count = np.count_nonzero(
                    self.xs_amps_ps[i_depth, i_cell, :])
                if nonzero_count != 0:
                    amps_ps = self.xs_amps_ps[i_depth, i_cell, 0:nonzero_count]
                    xs_ps_avg[i_depth, i_cell] = np.mean(amps_ps)

                nonzero_count = np.count_nonzero(
                    self.xs_amps_pps[i_depth, i_cell, :])
                if nonzero_count != 0:
                    amps_pps = self.xs_amps_pps[i_depth,
                                                i_cell, 0:nonzero_count]
                    xs_pps_avg[i_depth, i_cell] = np.mean(amps_pps)

                nonzero_count = np.count_nonzero(
                    self.xs_amps_pss[i_depth, i_cell, :])
                if nonzero_count != 0:
                    amps_pss = self.xs_amps_pss[i_depth,
                                                i_cell, 0:nonzero_count]
                    xs_pss_avg[i_depth, i_cell] = np.mean(amps_pss)

        self.xs_ps_avg = xs_ps_avg
        self.xs_pps_avg = xs_pps_avg
        self.xs_pss_avg = xs_pss_avg

    def gccp(self, wlen=15.):

        if not self.is_ready_for_gccp:
            raise(Exception("CCPimage not ready for gccp"))

        dlat = max(self.lateral_distances)/self.n_lateral

        import scipy.ndimage as ndimage

        self.xs_gauss_ps = ndimage.filters.gaussian_filter(
            self.xs_ps_avg, sigma=(0, wlen/dlat))
        self.xs_gauss_pps = ndimage.filters.gaussian_filter(
            self.xs_pps_avg, sigma=(0, wlen/dlat))
        self.xs_gauss_pss = ndimage.filters.gaussian_filter(
            self.xs_pss_avg, sigma=(0, wlen/dlat))

    def stack_ccp(self):

        tot_trace = np.zeros((self.n_depth, self.n_lateral))

        for i_cell in range(self.n_lateral):
            ps_trace = self.xs_ps_avg[:, i_cell]
            pps_trace = self.xs_pps_avg[:, i_cell]
            pss_trace = self.xs_pss_avg[:, i_cell]
            tot_trace[:, i_cell] = (ps_trace + pps_trace + pss_trace)

        self.tot_trace_ccp = tot_trace

    def pws_gccp(self):

        tot_trace = np.zeros((self.n_depth, self.n_lateral))

        for i_cell in range(self.n_lateral):
            ps_trace = self.xs_gauss_ps[:, i_cell]
            pps_trace = self.xs_gauss_pps[:, i_cell]
            pss_trace = self.xs_gauss_pss[:, i_cell]

            weight = np.zeros(len(ps_trace), dtype=complex)
            ps_hilb = hilbert(ps_trace)
            ps_phase = np.arctan2(ps_hilb.imag, ps_hilb.real)
            weight += np.exp(1j*ps_phase)

            pps_hilb = hilbert(pps_trace)
            pps_phase = np.arctan2(pps_hilb.imag, pps_hilb.real)
            weight += np.exp(1j*pps_phase)

            pss_hilb = hilbert(pss_trace)
            pss_phase = np.arctan2(pss_hilb.imag, pss_hilb.real)
            weight += np.exp(1j*pss_phase)

            weight = np.abs(weight)/3.

            tot_trace[:, i_cell] = (ps_trace + pps_trace + pss_trace)*weight**2

        self.tot_trace_gccp = tot_trace


    def save(self, title):

        if title is None:
            title = "CCP_image.pkl"

        if not ".pkl" in title:
            file = open(title+".pkl", "wb")
        else:
            file = open(title, "wb")
        pickle.dump(self, file)
        file.close()


    def plot_ccp(self, vmin=-0.05, vmax=0.05, save=False, form='png'):

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(
            4, 1, figsize=(8.5, 8))

        # plt.pcolormesh(lateral_distances,depth_array,xs_ps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im1 = ax1.pcolormesh(self.lateral_distances, self.depth_array,
                             self.xs_ps_avg, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im1, ax=ax1)
        ax1.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax1.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax1.set_xlabel('Lateral Distance (km)', size=10)
        ax1.set_ylabel('Depth (km)', size=10)
        ax1.set_title('Ps CCP image', size=10)
        ax1.invert_yaxis()

        # plt.pcolormesh(lateral_distances,depth_array,xs_pps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im2 = ax2.pcolormesh(self.lateral_distances, self.depth_array,
                             self.xs_pps_avg*3., cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im2, ax=ax2)
        ax2.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax2.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax2.set_xlabel('Lateral Distance (km)', size=10)
        ax2.set_ylabel('Depth (km)', size=10)
        ax2.set_title('Pps CCP image', size=10)
        ax2.invert_yaxis()

        im3 = ax3.pcolormesh(self.lateral_distances, self.depth_array,
                             self.xs_pss_avg*3., cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im3, ax=ax3)
        ax3.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax3.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax3.set_xlabel('Lateral Distance (km)', size=10)
        ax3.set_ylabel('Depth (km)', size=10)
        ax3.set_title('Pss CCP image', size=10)
        ax3.invert_yaxis()

        im4 = ax4.pcolormesh(self.lateral_distances, self.depth_array,
                             self.tot_trace_ccp, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im4, ax=ax4)
        ax4.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax4.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        ax4.set_xlabel('Lateral Distance (km)', size=10)
        ax4.set_ylabel('Depth (km)', size=10)
        ax4.set_title('Weighted CCP image', size=10)
        ax4.invert_yaxis()

        if save:
            if not os.path.isdir("FIGURES"):
                os.makedirs("FIGURES")
            plt.savefig('FIGURES/ccp.' + form)

        plt.tight_layout()
        plt.show()

    def plot_gccp(self, vmin=-0.015, vmax=0.015, save=False, form='png'):

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(
            4, 1, figsize=(8.5, 8))

        # plt.pcolormesh(lateral_distances,depth_array,xs_ps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im1 = ax1.pcolormesh(self.lateral_distances, self.depth_array,
                             self.xs_gauss_ps, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im1, ax=ax1)
        ax1.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax1.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax1.set_xlabel('Lateral Distance (km)', size=10)
        ax1.set_ylabel('Depth (km)', size=10)
        ax1.set_title('Ps GCCP image', size=10)
        ax1.invert_yaxis()

        # plt.pcolormesh(lateral_distances,depth_array,xs_pps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im2 = ax2.pcolormesh(self.lateral_distances, self.depth_array,
                             self.xs_gauss_pps*3, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im2, ax=ax2)
        ax2.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax2.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax2.set_xlabel('Lateral Distance (km)', size=10)
        ax2.set_ylabel('Depth (km)', size=10)
        ax2.set_title('Pps GCCP image', size=10)
        ax2.invert_yaxis()

        im3 = ax3.pcolormesh(self.lateral_distances, self.depth_array,
                             self.xs_gauss_pss*3, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im3, ax=ax3)
        ax3.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax3.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax3.set_xlabel('Lateral Distance (km)', size=10)
        ax3.set_ylabel('Depth (km)', size=10)
        ax3.set_title('Pss GCCP image', size=10)
        ax3.invert_yaxis()

        import scipy.ndimage as ndimage
        self.tot_trace_gccp = ndimage.filters.gaussian_filter(
            self.tot_trace_gccp, sigma=(3, 0), order=0)

        im4 = ax4.pcolormesh(self.lateral_distances, self.depth_array,
                             self.tot_trace_gccp, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im4, ax=ax4)
        ax4.set_xlim((min(self.lateral_distances)),
                     (max(self.lateral_distances)))
        ax4.set_ylim((min(self.depth_array)), (max(self.depth_array)))
        bar.ax.set_ylabel('Amplitude', size=10)
        ax4.set_xlabel('Lateral Distance (km)', size=10)
        ax4.set_ylabel('Depth (km)', size=10)
        ax4.set_title('Phase-weighted GCCP image', size=10)
        ax4.invert_yaxis()

        if save:
            if not os.path.isdir("FIGURES"):
                os.makedirs("FIGURES")
            plt.savefig('FIGURES/gccp.' + form)

        plt.tight_layout()
        plt.show()


def ppoint_distance(tr, dz, vs):
    """
    Calculate horizontal distance for interval dz and velocity vs

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Single trace object to migrate to depth
    dz : float
        Vertical sampling distance
    vs : float
        S-wave velocity (km/s)

    """

    # Calculate distance
    dx = dz*np.tan(np.arcsin(tr.stats.slow*vs))

    return dx


def ppoint(tr, dist):
    """
    Determine geographic location of piercing point

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Single trace object to migrate to depth
    dist : float
        Horizontal istance from the station (km)

    """

    # Conversion factors
    lat2km = 111.
    lon2km = 90.

    # Get lat and lon of station location
    slat = tr.stats.stla
    slon = tr.stats.stlo

    # Back-azimuth of event
    baz = tr.stats.baz*np.pi/180.

    # location of piercing point on geographical grid
    plat = dist*np.sin(-baz+np.pi/2.)/lat2km + slat
    plon = dist*np.cos(-baz+np.pi/2.)/lon2km + slon

    return plon, plat


def ttime(tr, dz, vp, vs, phase=None):
    """
    Calculate travel time for interval dz and velocities vp and vs

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Single trace object to migrate to depth
    dz : float
        Vertical sampling distance (km)
    vp : float
        P-wave velocity
    vs : float
        S-wave velocity
    phase : str
        Phase of interest
    """

    # Get horizontal slowness
    slow = tr.stats.slow

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


def timeshift(tr, tt):
    """
    Shift a trace by a travel time tt and take amplitude at zero

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Single trace object to migrate to depth
    tt : float
        Travel time (sec)

    """

    # Define frequencies
    nt = int(tr.stats.npts)
    dt = tr.stats.delta
    freq = np.fft.fftfreq(int(nt), d=dt)

    # Hilbert transform and instantaneous phase
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
    """
    Calculate travel times through velocity model for all phases of interest

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Single trace object to migrate to depth
    nz : int
        Number of layers in interpolation
    dep : :class:`~numpy.ndarray`
        Depth array for velocity model
    vp : :class:`~numpy.ndarray`
        P-wave velocity array for velocity model
    vs : :class:`~numpy.ndarray`
        S-wave velocity array for velocity model

    """

    # Define arrays with zeros
    plat = np.zeros(nz)
    plon = np.zeros(nz)
    ttps = np.zeros(nz)
    ttpps = np.zeros(nz)
    ttpss = np.zeros(nz)

    # Get regular depth array
    idep = np.linspace(dep.min(), dep.max(), nz)

    # Interpolate Vp and Vs models on depth grid
    ivp = sp.interpolate.interp1d(dep, vp, kind='linear')(idep)
    ivs = sp.interpolate.interp1d(dep, vs, kind='linear')(idep)

    # Get exact depth interval
    dz = idep[1] - idep[0]

    # Now loop through all depths
    for iz in range(nz):

        # Initialize travel time and distance counters
        dtps = 0.
        dtpps = 0.
        dtpss = 0.
        dx = 0.

        # Sum over depths from 0 to iz
        for i in range(iz):
            dtps += ttime(tr, dz, ivp[i], ivs[i], 'Ps')
            dtpps += ttime(tr, dz, ivp[i], ivs[i], 'Pps')
            dtpss += ttime(tr, dz, ivp[i], ivs[i], 'Pss')
            dx += ppoint_distance(tr, dz, ivs[i])

        # Get piercing point from distance
        plo, pla = ppoint(tr, dx)

        # Assign values to arrays
        ttps[iz] = dtps
        ttpps[iz] = dtpps
        ttpss[iz] = dtpss
        plon[iz] = plo
        plat[iz] = pla

    return ttps, ttpps, ttpss, plon, plat, idep


def haversine(lat, lon, xs_lat, xs_lon):  # great-circle distance (kilometres)

    earth_radius = 6371.  # kilometres

    lat = np.radians(lat)
    lon = np.radians(lon)
    xs_lat = np.radians(xs_lat)
    xs_lon = np.radians(xs_lon)
    dlat = lat - xs_lat
    dlon = lon - xs_lon
    a = ((np.sin(dlat/2.))**2.) + \
        (np.cos(xs_lat)*np.cos(lat)*((np.sin(dlon/2.))**2.))
    distance = np.abs(2.*earth_radius*np.arcsin(np.sqrt(a)), dtype=float)

    return np.abs(distance)


def _progressbar(it, prefix="", size=60, file=sys.stdout):
    """
    Show progress bar while looping in for loop

    """

    count = len(it)

    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" %
                   (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()
