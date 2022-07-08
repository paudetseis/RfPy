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
    """
    A CCPimage object contains attributes and methods to produce
    Common Conversion Point (CCP) stacks for each of the main
    three main phases (Ps, Pps and Pss) using radial-component
    receiver functions. The object is used to project the stacks along
    a linear profile, specified by start and end geographical coordinate 
    locations, which are subsequently averaged to produce a final 
    CCP image. The averaging can be done using a linear weighted sum, 
    or a phase-weighted sum. Methods should be used in the appropriate
    sequence (see ``rfpy_ccp.py`` for details).

    Note
    ----
    By default, the object initializes with un-defined coordinate locations for 
    the profile. If not specified during initialization, make sure
    they are specified later, before the other methods are used, e.g.
    ``ccpimage.xs_lat1 = 10.; ccpimage.xs_lon1 = 110.``, etc. Note also that the 
    default 1D velocity model may not be applicable to your region of 
    interest and a different model can be implemented during initialization
    or later during processing. 

    Parameters
    ----------
    coord_start : list
        List of two floats corresponding to the (latitude, longitude)
        pair for the start point of the profile 
    coord_end : list
        List of two floats corresponding to the (latitude, longitude)
        pair for the end point of the profile 
    weights : list
        List of three floats with corresponding weights for the Ps, Pps
        and Pss phases used during linear, weighted averaging
    dep : :class:`~numpy.ndarray`
        Array of depth values defining the 1D background seismic velocity model.
        Note that the maximum depth defined here sets the maximum depth 
        in each of the CCP stacks and the final CCP image.
    vp : :class:`~numpy.ndarray`
        Array of Vp values defining the 1D background seismic velocity model
    vpvs : float
        Constant Vp/Vs ratio for the 1D model. 

    Other Parameters
    ----------------
    radialRF : list
        List of :class:`~obspy.core.Stream` objects containing the radial
        receiver functions along the line. Each item in the list contains the
        streams for one single station.
    vs : :class:`~numpy.ndarray`
        Array of Vp values defining the 1D background seismic velocity model
    xs_lat1 : float
        Latitude of start point defining the linear profile.
    xs_lon1 : float
        Longitude of start point defining the linear profile.
    xs_lat2 : float
        Latitude of end point defining the linear profile.
    xs_lon2 : float
        Longitude of end point defining the linear profile.
    is_ready_for_prep : boolean
        Whether or not the object is ready for the method ``prep_data``
    is_ready_for_prestack : boolean
        Whether or not the object is ready for the method ``prestack``
    is_ready_for_ccp : boolean
        Whether or not the object is ready for the method ``ccp``
    is_ready_for_gccp : boolean
        Whether or not the object is ready for the method ``gccp``

    """

    def __init__(self, coord_start=[None, None], coord_end=[None, None],
                 weights=[1., 3., -3.],
                 dep=np.array([0., 4., 8., 14., 30., 35., 45., 110.]),
                 vp=np.array([4.0, 5.9, 6.2, 6.3, 6.8, 7.2, 8.0, 8.1]),
                 vs=None, vpvs=1.73, dx=2.5, dz=1.):

        self.radialRF = []
        self.dep = dep
        self.weights = weights
        self.xs_lat1 = coord_start[0]
        self.xs_lon1 = coord_start[1]
        self.xs_lat2 = coord_end[0]
        self.xs_lon2 = coord_end[1]
        self.is_ready_for_prep = False
        self.is_ready_for_prestack = False
        self.is_ready_for_ccp = False
        self.is_ready_for_gccp = False

        # Define grid parameters
        self.dz = dz
        self.dx = dx

        # Get total length of grid from end points
        xlength = haversine(self.xs_lat1, self.xs_lon1,
                              self.xs_lat2, self.xs_lon2)

        # number of cells laterally and vertically
        self.nx = int(np.rint(xlength/self.dx))
        self.nz = int(self.dep[-1]/self.dz)

        self.xarray = np.arange(self.nx)*self.dx
        self.zarray = np.arange(self.nz)*self.dz

        # Interpolate Vp and Vs models on depth grid
        if vs is None:
            vs = vp/vpvs
        else:
            if not vp.shape == vs.shape:
                raise(Exception("vp and vs arrays have a different shape"))

        self.vp = sp.interpolate.interp1d(dep, vp, kind='linear')(self.zarray)
        self.vs = sp.interpolate.interp1d(dep, vs, kind='linear')(self.zarray)

    def add_rfstream(self, rfstream):
        """
        Method to add a :class:`~obspy.core.Stream` object to the list
        ``radialRF``. With at least one stream in ``radialRF``, the object
        is ready for the ``prep_data`` method, and the corresponding flag is
        updated. 

        Parameters
        ----------
        rfstream : :class:`~obspy.core.Stream`
            Stream object containing radial receiver functions for one station

        """
        if len(rfstream) > 0:
            # fftshift if the time axis starts at negative lags 
            if rfstream[0].stats.taxis[0]<0.:
                for tr in rfstream:
                    tr.data = np.fft.fftshift(tr.data)

            self.radialRF.append(rfstream)
            self.is_ready_for_prep = True

    def prep_data(self, f1=0.05, f2ps=0.5, f2pps=0.25, f2pss=0.2,
                  nbaz=36+1, nslow=40+1):
        """
        Method to pre-process the data and calculate the CCP points for each 
        of the receiver functions. Pre-processing includes the binning to
        back-azimuth and slowness bins to reduce processing time and generate
        cleaner stacks, as well as filtering to emphasize the energy of the
        various phases as a function of depth. As a general rule, the high 
        frequency corner of the 2 reverberated phases (Pps and Pss) should be 
        approximately half of the high frequency corner of the direct 
        (Ps) phase. At the end of this step, the object is updated with
        the amplitude at the lat, lon and depth values corresponding to the
        raypath for each receiver function. The object is now ready for the
        method ``prestack``, with the corresponding flag updated.

        Parameters
        ----------
        f1 : float
            Low-frequency corner of the bandpass filter used for all phases (Hz)
        f2ps : float
            High-frequency corner of the bandpass filter used for the Ps phase (Hz)
        f2pps : float
            High-frequency corner of the bandpass filter used for the Pps phase (Hz)
        f2pss : float
            High-frequency corner of the bandpass filter used for the Pss phase (Hz)
        nbaz : int
            Number of increments in the back-azimuth bins
        nslow : int
            Number of increments in the slowness bins

        The following attributes are added to the object:

        Other Parameters
        ----------------
        amp_ps_depth : :class:`numpy.ndarray`
            2D array of amplitudes as a function of depth for the Ps phase
        amp_pps_depth : :class:`numpy.ndarray`
            2D array of amplitudes as a function of depth for the Pps phase
        amp_pss_depth : :class:`numpy.ndarray`
            2D array of amplitudes as a function of depth for the Pss phase
        lon_depth : :class:`numpy.ndarray`
            2D array of longitude as a function of depth (i.e., piercing points)
        lat_depth : :class:`numpy.ndarray`
            2D array of latitude as a function of depth (i.e., piercing points)
        is_ready_for_presstack : boolean
            Flag specifying that the object is ready for the presstack() method
        n_traces : float
            Total number of traces used in creating the CCP image.

        """

        if not self.is_ready_for_prep:
            raise(Exception("CCPimage not ready for pre-prep"))

        ikey = 0
        total_traces = 0

        # Process streams one at a time
        for RF in self.radialRF:

            # Bin RFs into back-azimuth and slowness bins to speed up
            # calculations
            RFbin = binning.bin_baz_slow(
                RF, nbaz=nbaz, nslow=nslow)[0]
            n_traces = len(RFbin)
            total_traces += n_traces

            amp_ps_tr = np.empty([n_traces, self.nz])
            amp_pps_tr = np.empty([n_traces, self.nz])
            amp_pss_tr = np.empty([n_traces, self.nz])
            lon_tr = np.empty([n_traces, self.nz])
            lat_tr = np.empty([n_traces, self.nz])

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
                tt_ps, tt_pps, tt_pss, plon, plat = \
                    raypath(st_ps[itr], dep=self.zarray, vp=self.vp, vs=self.vs)

                # Now get amplitude of RF at corresponding travel
                # time along the raypath
                lon_tr[itr, :] = plon
                lat_tr[itr, :] = plat

                amp_ps = []
                amp_pps = []
                amp_pss = []

                # Loop through travel times and shift RFs to get amplitudes
                for tt in tt_ps:
                    a, phase = timeshift(st_ps[itr], tt)
                    amp_ps.append(a)
                amp_ps_tr[itr, :] = amp_ps

                # Loop through travel times and shift RFs to get amplitudes
                for tt in tt_pps:
                    a, phase = timeshift(st_pps[itr], tt)
                    amp_pps.append(a)
                amp_pps_tr[itr, :] = amp_pps

                # Loop through travel times and shift RFs to get amplitudes
                for tt in tt_pss:
                    a, phase = timeshift(st_pss[itr], tt)
                    amp_pss.append(a)
                amp_pss_tr[itr, :] = amp_pss

            if ikey == 0:
                amp_ps_depth = amp_ps_tr.transpose()
                amp_pps_depth = amp_pps_tr.transpose()
                amp_pss_depth = amp_pss_tr.transpose()
                lon_depth = lon_tr.transpose()
                lat_depth = lat_tr.transpose()

            elif ikey > 0:
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

            ikey += 1

        self.amp_ps_depth = amp_ps_depth
        self.amp_pps_depth = amp_pps_depth
        self.amp_pss_depth = amp_pss_depth
        self.lon_depth = lon_depth
        self.lat_depth = lat_depth
        self.is_ready_for_prestack = True
        self.n_traces = total_traces

        del self.radialRF

    def prestack(self):
        """
        Method to project the raypaths onto the 2D profile for each of the three
        phases. The final grid is defined here, using the parameter ``dx`` in km.
        The horizontal extent is pre-determined from the start and end points of 
        the profile. At the end of this step, the object contains the set of 
        amplitudes at each of the 2D grid points, for each of the three phases.
        The object is now ready for the methods ``ccp`` and/or ``gccp``, with 
        the corresponding flag updated.

        The following attributes are added to the object:

        Other Parameters
        ----------------
        xs_amps_ps : :class:`numpy.ndarray`
            3D array of amplitudes (1D array of amplitudes at each grid cell)
            for the Ps phase
        xs_amps_pps : :class:`numpy.ndarray`
            3D array of amplitudes (1D array of amplitudes at each grid cell)
            for the Pps phase
        xs_amps_pss : :class:`numpy.ndarray`
            3D array of amplitudes (1D array of amplitudes at each grid cell)
            for the Pss phase
        is_ready_for_ccp : boolean
            Flag specifying that the object is ready for the ccp() method
        is_ready_for_gccp : boolean
            Flag specifying that the object is ready for the gccp() method

        """

        if not self.is_ready_for_prestack:
            raise(Exception("CCPimage not ready for prestack"))

        xs_latitudes = np.asarray(
            np.linspace(self.xs_lat1, self.xs_lat2, self.nx))
        xs_longitudes = np.asarray(
            np.linspace(self.xs_lon1, self.xs_lon2, self.nx))

        xs_amps_ps = np.zeros((self.nz, self.nx, self.n_traces))
        xs_amps_pps = np.zeros((self.nz, self.nx, self.n_traces))
        xs_amps_pss = np.zeros((self.nz, self.nx, self.n_traces))

        for iz in _progressbar(range(self.nz), '', 25):

            for i_coor in range(self.n_traces):

                lat_tr = self.lat_depth[iz, i_coor]
                lon_tr = self.lon_depth[iz, i_coor]
                distance_tests = np.empty(self.nx)

                for i_xs in range(self.nx):
                    lat_xs = xs_latitudes[i_xs]
                    lon_xs = xs_longitudes[i_xs]
                    distance_tests[i_xs] = haversine(
                        lat_xs, lon_xs, lat_tr, lon_tr)

                minimum_distance = np.amin(distance_tests)
                ix = np.where(distance_tests ==
                                  np.amin(distance_tests))[0][0]

                nonzero_count = np.count_nonzero(
                    xs_amps_ps[iz, ix, :])
                new_amp_ps = self.amp_ps_depth[iz, i_coor]
                if xs_amps_ps[iz, ix, 0] == 0.:
                    xs_amps_ps[iz, ix, 0] = new_amp_ps
                else:
                    xs_amps_ps[iz, ix, nonzero_count] = new_amp_ps

                nonzero_count = np.count_nonzero(
                    xs_amps_pps[iz, ix, :])
                new_amp_pps = self.amp_pps_depth[iz, i_coor]
                if xs_amps_pps[iz, ix, 0] == 0.:
                    xs_amps_pps[iz, ix, 0] = new_amp_pps
                else:
                    xs_amps_pps[iz, ix, nonzero_count] = new_amp_pps

                nonzero_count = np.count_nonzero(
                    xs_amps_pss[iz, ix, :])
                new_amp_pss = self.amp_pss_depth[iz, i_coor]
                if xs_amps_pss[iz, ix, 0] == 0.:
                    xs_amps_pss[iz, ix, 0] = new_amp_pss
                else:
                    xs_amps_pss[iz, ix, nonzero_count] = new_amp_pss

        self.xs_amps_ps = xs_amps_ps
        self.xs_amps_pps = xs_amps_pps
        self.xs_amps_pss = xs_amps_pss
        self.is_ready_for_ccp = True
        self.is_ready_for_gccp = True

        del self.amp_ps_depth
        del self.amp_pps_depth
        del self.amp_pss_depth
        del self.lon_depth
        del self.lat_depth

    def ccp(self):
        """
        Method to average the amplitudes at each grid point to produce 2D images
        for each of the three phases. At the end of this step, the object
        contains the three 2D arrays that can be further averaged into a single
        final image. 

        The following attributes are added to the object:

        Other Parameters
        ----------------
        xs_ps_avg : :class:`numpy.ndarray`
            2D array of stacked amplitudes for the Ps phase
        xs_pps_avg : :class:`numpy.ndarray`
            2D array of stacked amplitudes for the Pps phase
        xs_pss_avg : :class:`numpy.ndarray`
            2D array of stacked amplitudes for the Pss phase

        """

        if not self.is_ready_for_ccp:
            raise(Exception("CCPimage not ready for ccp"))

        xs_ps_avg = np.zeros((self.nz, self.nx))
        xs_pps_avg = np.zeros((self.nz, self.nx))
        xs_pss_avg = np.zeros((self.nz, self.nx))

        for iz in _progressbar(range(self.nz), '', 25):

            for ix in range(self.nx):

                nonzero_count = np.count_nonzero(self.xs_amps_ps[iz, ix, :])
                if nonzero_count != 0:
                    amps_ps = self.xs_amps_ps[iz, ix, 0:nonzero_count]
                    xs_ps_avg[iz, ix] = np.mean(amps_ps)

                nonzero_count = np.count_nonzero(self.xs_amps_pps[iz, ix, :])
                if nonzero_count != 0:
                    amps_pps = self.xs_amps_pps[iz, ix, 0:nonzero_count]
                    xs_pps_avg[iz, ix] = np.mean(amps_pps)

                nonzero_count = np.count_nonzero(self.xs_amps_pss[iz, ix, :])
                if nonzero_count != 0:
                    amps_pss = self.xs_amps_pss[iz, ix, 0:nonzero_count]
                    xs_pss_avg[iz, ix] = np.mean(amps_pss)

        self.xs_ps_avg = xs_ps_avg
        self.xs_pps_avg = xs_pps_avg
        self.xs_pss_avg = xs_pss_avg

        del self.xs_amps_ps
        del self.xs_amps_pps
        del self.xs_amps_pss

    def gccp(self, wlen=15.):
        """
        Method to average the amplitudes at each grid point to produce 2D images
        for each of the three phases. In this method, the grid points are further
        smoothed in the horizontal direction using a Gaussian function to simulate
        P-wave sensitivity kernels. At the end of this step, the object
        contains the three 2D smoothed arrays that can be further averaged into a 
        single final image. 

        Parameters
        ----------
        wlen : float
            Wavelength of the P-wave for smoothing (km).

        The following attributes are added to the object:

        Other Parameters
        ----------------
        xs_gauss_ps : :class:`numpy.ndarray`
            2D array of stacked and Gaussian-filtered amplitudes for the Ps phase
        xs_gauss_pps : :class:`numpy.ndarray`
            2D array of stacked and Gaussian-filtered amplitudes for the Pps phase
        xs_gauss_pss : :class:`numpy.ndarray`
            2D array of stacked and Gaussian-filtered amplitudes for the Pss phase

        """

        if not self.is_ready_for_gccp:
            raise(Exception("CCPimage not ready for gccp"))
        if not hasattr(self, 'xs_ps_avg'):
            self.ccp()

        dlat = max(self.xarray)/self.nx

        import scipy.ndimage as ndimage

        self.xs_gauss_ps = ndimage.filters.gaussian_filter(
            self.xs_ps_avg, sigma=(0, wlen/dlat))
        self.xs_gauss_pps = ndimage.filters.gaussian_filter(
            self.xs_pps_avg, sigma=(0, wlen/dlat))
        self.xs_gauss_pss = ndimage.filters.gaussian_filter(
            self.xs_pss_avg, sigma=(0, wlen/dlat))

    def linear_stack(self, typ='ccp'):
        """
        Method to average the three 2D images into a final, weighted CCP image
        using the weights defined in the attribute.

        Parameters
        ----------
        typ : str
            Type of phase stacks to use (either `ccp` or `gccp`)

        The following attributes are added to the object:

        Other Parameters
        ----------------
        tot_trace : :class:`numpy.ndarray`
            2D array of amplitudes for the linearly combined Ps, Pps and Pss phases

        """

        tot_trace = np.zeros((self.nz, self.nx))

        if typ=='ccp':
            if not hasattr(self, "xs_ps_avg"):
                self.ccp()
            xs_ps = self.xs_ps_avg*self.weights[0]
            xs_pps = self.xs_pps_avg*self.weights[1]
            xs_pss = self.xs_pss_avg*self.weights[2]
        elif typ=='gccp':
            if not hasattr(self, 'xs_gauss_ps'):
                self.gccp()
            xs_ps = self.xs_gauss_ps*self.weights[0]
            xs_pps = self.xs_gauss_pps*self.weights[1]
            xs_pss = self.xs_gauss_pss*self.weights[2]

        for ix in range(self.nx):
            ps_trace = xs_ps[:, ix]
            pps_trace = xs_pps[:, ix]
            pss_trace = xs_pss[:, ix]
            tot_trace[:, ix] = (ps_trace + pps_trace + pss_trace)

        self.tot_trace = tot_trace

    def phase_weighted_stack(self, typ='gccp'):
        """
        Method to average the three 2D smoothed images into a final, 
        phase-weighted CCP image.

        Parameters
        ----------
        typ : str
            Type of phase stacks to use (either `ccp` or `gccp`)

        The following attributes are added to the object:

        Other Parameters
        ----------------
        tot_trace : :class:`numpy.ndarray`
            2D array of amplitudes for the phase-weighted, combined 
            Ps, Pps and Pss phases

        """

        tot_trace = np.zeros((self.nz, self.nx))

        if typ=='ccp':
            if not hasattr(self, "xs_ps_avg"):
                self.ccp()
            xs_ps = self.xs_ps_avg*self.weights[0]
            xs_pps = self.xs_pps_avg*self.weights[1]
            xs_pss = self.xs_pss_avg*self.weights[2]
        elif typ=='gccp':
            if not hasattr(self, 'xs_gauss_ps'):
                self.gccp()
            xs_ps = self.xs_gauss_ps*self.weights[0]
            xs_pps = self.xs_gauss_pps*self.weights[1]
            xs_pss = self.xs_gauss_pss*self.weights[2]

        for ix in range(self.nx):
            ps_trace = xs_ps[:, ix]
            pps_trace = xs_pps[:, ix]
            pss_trace = xs_pss[:, ix]

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

            tot_trace[:, ix] = (ps_trace + pps_trace + pss_trace)*weight**2

        self.tot_trace = tot_trace


    def save(self, title):
        """
        Method to save the `ccpimage` object to file.

        Parameters
        ----------
        title : str
            File name for the saved object

        """
        title = str(title)

        if title is None:
            title = "CCP_image.pkl"

        if not ".pkl" in title:
            file = open(title+".pkl", "wb")
        else:
            file = open(title, "wb")
        pickle.dump(self, file)
        file.close()


    def plot_ccp(self, vmin=-0.05, vmax=0.05, 
        save=False, title='', fmt='png'):
        """
        Method to plot the final CCP stacks along the line.

        Parameters
        ----------
        vmin : float
            Maximum negative value for the color scale
        vmax : float
            Maximum positive value for the color scale
        save : boolean
            Whether or not to save the figure.
        fmt : str
            Format of saved figure
        xlen : float
            Total length of profile

        """

        xlen = self.nx*self.dx
        xm = 1.5 + 0.0065*xlen
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(
            4, 1, figsize=(xm, 8))

        # plt.pcolormesh(xarray,zarray,xs_ps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im1 = ax1.pcolormesh(self.xarray, self.zarray,
                             self.xs_ps_avg*self.weights[0], 
                             cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im1, ax=ax1)
        ax1.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax1.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax1.set_xlabel('Lateral Distance (km)', size=10)
        ax1.set_ylabel('Depth (km)', size=10)
        ax1.set_title('Ps CCP image', size=10)
        ax1.invert_yaxis()

        # plt.pcolormesh(xarray,zarray,xs_pps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im2 = ax2.pcolormesh(self.xarray, self.zarray,
                             self.xs_pps_avg*self.weights[1], 
                             cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im2, ax=ax2)
        ax2.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax2.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax2.set_xlabel('Lateral Distance (km)', size=10)
        ax2.set_ylabel('Depth (km)', size=10)
        ax2.set_title('Pps CCP image', size=10)
        ax2.invert_yaxis()

        im3 = ax3.pcolormesh(self.xarray, self.zarray,
                             self.xs_pss_avg*self.weights[2], 
                             cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im3, ax=ax3)
        ax3.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax3.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax3.set_xlabel('Lateral Distance (km)', size=10)
        ax3.set_ylabel('Depth (km)', size=10)
        ax3.set_title('Pss CCP image', size=10)
        ax3.invert_yaxis()

        im4 = ax4.pcolormesh(self.xarray, self.zarray,
                             self.tot_trace, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im4, ax=ax4)
        ax4.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax4.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        ax4.set_xlabel('Lateral Distance (km)', size=10)
        ax4.set_ylabel('Depth (km)', size=10)
        ax4.set_title('Weighted CCP image', size=10)
        ax4.invert_yaxis()
        plt.tight_layout()

        if save:
            if not os.path.isdir("FIGURES"):
                os.makedirs("FIGURES")
            plt.savefig('FIGURES/ccp.' + title + '.' + fmt)

        plt.show()

    def plot_gccp(self, vmin=-0.015, vmax=0.015, 
        save=False, title='', fmt='png'):
        """
        Method to plot the final GCCP stacks along the line.

        Parameters
        ----------
        vmin : float
            Maximum negative value for the color scale
        vmax : float
            Maximum positive value for the color scale
        save : boolean
            Whether or not to save the figure.
        fmt : str
            Format of saved figure
        xlen : float
            Total length of profile

        """

        xlen = self.nx*self.dx
        xm = 1.5 + 0.0065*xlen
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(
            4, 1, figsize=(xm, 8))

        # plt.pcolormesh(xarray,zarray,xs_ps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im1 = ax1.pcolormesh(self.xarray, self.zarray,
                             self.xs_gauss_ps*self.weights[0], 
                             cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im1, ax=ax1)
        ax1.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax1.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax1.set_xlabel('Lateral Distance (km)', size=10)
        ax1.set_ylabel('Depth (km)', size=10)
        ax1.set_title('Ps GCCP image', size=10)
        ax1.invert_yaxis()

        # plt.pcolormesh(xarray,zarray,xs_pps_avg,cmap=cm.coolwarm,vmin=vmin,vmax=vmax)
        im2 = ax2.pcolormesh(self.xarray, self.zarray,
                             self.xs_gauss_pps*self.weights[1], 
                             cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im2, ax=ax2)
        ax2.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax2.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax2.set_xlabel('Lateral Distance (km)', size=10)
        ax2.set_ylabel('Depth (km)', size=10)
        ax2.set_title('Pps GCCP image', size=10)
        ax2.invert_yaxis()

        im3 = ax3.pcolormesh(self.xarray, self.zarray,
                             self.xs_gauss_pss*self.weights[2], 
                             cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im3, ax=ax3)
        ax3.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax3.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        # ax3.set_xlabel('Lateral Distance (km)', size=10)
        ax3.set_ylabel('Depth (km)', size=10)
        ax3.set_title('Pss GCCP image', size=10)
        ax3.invert_yaxis()

        import scipy.ndimage as ndimage
        self.tot_trace = ndimage.filters.gaussian_filter(
            self.tot_trace, sigma=(3, 0), order=0)

        im4 = ax4.pcolormesh(self.xarray, self.zarray,
                             self.tot_trace, cmap=cm.RdBu_r,
                             vmin=vmin, vmax=vmax)
        bar = plt.colorbar(im4, ax=ax4)
        ax4.set_xlim((min(self.xarray)),
                     (max(self.xarray)))
        ax4.set_ylim((min(self.zarray)), (max(self.zarray)))
        bar.ax.set_ylabel('Amplitude', size=10)
        ax4.set_xlabel('Lateral Distance (km)', size=10)
        ax4.set_ylabel('Depth (km)', size=10)
        ax4.set_title('Phase-weighted GCCP image', size=10)
        ax4.invert_yaxis()
        plt.tight_layout()

        if save:
            if not os.path.isdir("FIGURES"):
                os.makedirs("FIGURES")
            plt.savefig('FIGURES/gccp.' + title + '.' + fmt)

        plt.show()


def ppoint_distance(tr, delta_z, vs):
    """
    Calculate horizontal distance for interval delta_z and velocity vs

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Single trace object to migrate to depth
    delta_z : float
        Vertical sampling distance
    vs : float
        S-wave velocity (km/s)

    """

    # Calculate distance
    delta_x = delta_z*np.tan(np.arcsin(tr.stats.slow*vs))

    return delta_x


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


def ttime(tr, delta_z, vp, vs, phase=None):
    """
    Calculate travel time for interval delta_z and velocities vp and vs

    Parameters
    ----------
    tr : :class:`~obspy.core.Trace`
        Single trace object to migrate to depth
    delta_z : float
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
        tt = delta_z*(np.sqrt((1./vs)**2 - slow**2) -
                 np.sqrt((1./vp)**2 - slow**2))
    elif phase == 'Pps':
        tt = delta_z*(np.sqrt((1./vs)**2 - slow**2) +
                 np.sqrt((1./vp)**2 - slow**2))
    elif phase == 'Pss':
        tt = 2.*delta_z*(np.sqrt((1./vs)**2 - slow**2))
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


def raypath(tr, dep=None, vp=None, vs=None):
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

    # Get exact depth parameters
    delta_z = dep[1] - dep[0]
    nz = len(dep)

    # Define arrays with zeros
    plat = np.zeros(nz)
    plon = np.zeros(nz)
    ttps = np.zeros(nz)
    ttpps = np.zeros(nz)
    ttpss = np.zeros(nz)

    # Now loop through all depths
    for iz in range(nz):

        # Initialize travel time and distance counters
        dtps = 0.
        dtpps = 0.
        dtpss = 0.
        delta_x = 0.

        # Sum over depths from 0 to iz
        for i in range(iz):
            dtps += ttime(tr, delta_z, vp[i], vs[i], 'Ps')
            dtpps += ttime(tr, delta_z, vp[i], vs[i], 'Pps')
            dtpss += ttime(tr, delta_z, vp[i], vs[i], 'Pss')
            delta_x += ppoint_distance(tr, delta_z, vs[i])

        # Get piercing point from distance
        plo, pla = ppoint(tr, delta_x)

        # Assign values to arrays
        ttps[iz] = dtps
        ttpps[iz] = dtpps
        ttpss[iz] = dtpss
        plon[iz] = plo
        plat[iz] = pla

    return ttps, ttpps, ttpss, plon, plat


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
