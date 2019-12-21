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
Functions to calculate thickness (H) and Vp/Vs ratio (R) of the crust based on
moveout times of direct Ps and reverberated Pps and Pss phases.

The stacks are obtained from the median weighted by the phase of individual signals.

"""

import numpy as np
from obspy.core import Stream, Trace, AttribDict
from matplotlib import pyplot as plt
from scipy.signal import hilbert
from scipy import stats
import sys


class HkStack(object):
    """
    A HkStack object contains attributes and methods to stack radial
    receiver functions along moveout curves for point measurements 
    of the depth to Moho (H) and P-to-S velocity ratio (k) beneath
    a seismic stations. The object is initialized with at least one
    :class:`~obspy.core.Stream` containing observed (or synthetised)
    radial receiver functions. The methods available can produce linear
    weighted stacks or product stacks, and can be used in the presence
    of flat or dipping Moho (with known strike and dip). 

    Note
    ----
    The object is initialized with the ``rfV1`` field only, and
    other attributes are added to the object as the analysis proceeds.
    A second ``rfV2`` can be included, which is typically a copy of ``rfV1``
    filtered at different corner frequencies and is used to stack along the
    Pps and Pss moveout curves.

    Parameters
    ----------
    rfV1 : :class:`~obspy.core.Stream`
        Stream object containing the radial-component receiver function
        seismograms 
    rfV2 : :class:`~obspy.core.Stream`
        Stream object containing the radial-component receiver function
        seismograms (typically filtered at lower frequencies)
    strike : float
        Strike angle of dipping Moho (has to be known or estimated a priori)
    dip : float
        Dip angle of Moho (has to be known or estimated a priori)
    vp : float
        Mean P-wave velocity of the crust (km/s)

    Default Parameters
    ------------------
    kbound : list
        List of 2 floats that determine the range of Vp/Vs values to search
    dk : float
        Spacing between adjacent Vp/Vs search values
    hbound : list
        List of 2 floats that determine the range of Moho depth values to search
    dh : float
        Spacing between adjacent Moho depth search values
    weights : list
        List of 3 floats that determine the relative weights of the individual
        phase stacks to be used in the final stack. The third weight is negative
        since Pss amplitudes are opposite to those of the Ps and Pps phases.
    phases : list
        List of 3 strings ('ps', 'pps', 'pss') corresponding to the thre phases
        of interest (`do not modify this attribute`)

    Examples
    --------


    """

    def __init__(self, rfV1, rfV2=None, strike=None, dip=None, vp=6.0):

        self.rfV1 = rfV1
        self.rfV2 = rfV2
        self.strike = strike
        self.dip = dip
        self.vp = vp
        self.kbound = [1.56, 2.1]
        self.dk = 0.02
        self.hbound = [20., 50.]
        self.dh = 0.5
        self.weights = [0.5, 2., -1.]
        self.phases = ['ps', 'pps', 'pss']

    def stack(self, vp=None):
        """
        Method to calculate Hk stacks from radial receiver functions.
        The stacks are calculated using phase-weighted stacking for
        individual phases and take the median of the weighted stack 
        to avoid bias by outliers.
        
        Note
        ----
        If two streams are available as attributes, the method will assume
        that the second stream should be used for stacking along the Pps
        and Pss move out curves (e.g., if the second stream contains
        lower frequency signals). 

        Parameters
        ----------
        vp : float
            Mean crust P-wave velocity. If this parameter is not specified,
            the method will search the available stream for this attribute. 
            If not found, it will use the value set during initialization
            (default of 6.0 km/s)

        Attributes
        ----------
        pws : :class:`~numpy.ndarray`
            Array of phase stacks, where the outer dimension corresponds
            to the phase index (shape ``nH, nk, nph``)
        sig : :class:`~numpy.ndarray`
            Variance of phase stacks, where the outer dimension corresponds
            to the phase index (shape ``nH, nk, nph``)

        Example
        -------
        
        """

        # Mean crustal P-wave velocity
        if not vp:
            try:
                vp = self.rfV1[0].stats.vp
            except:
                vp = self.vp

        # Station name
        sta = self.rfV1[0].stats.station

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        # Initialize arrays
        weight = np.complex(0.)
        amp = np.zeros(len(self.rfV1))
        sig = np.zeros((len(H), len(k), len(self.phases)))
        pws = np.zeros((len(H), len(k), len(self.phases)))

        for ih in _progressbar(range(len(H)), 'Computing: ', 15):
            for ik, kk in enumerate(k):
                for ip, ph in enumerate(self.phases):
                    for i in range(len(self.rfV1)):

                        if self.rfV2 and (ph=='pps' or ph=='pss'):
                            rfV = self.rfV2[i].copy()
                        else:
                            rfV = self.rfV1[i].copy()

                        # Calculate move out for each phase and get
                        # median value, weighted by instantaneous phase (pws)
                        tt = _dtime_(rfV, H[ih], kk, vp, ph)
                        trace = _timeshift_(rfV, tt)
                        thilb = hilbert(trace)
                        tphase = np.arctan2(thilb.imag, thilb.real)
                        weight += np.exp(1j*tphase[0])
                        amp[i] = trace[0]

                    weight = abs(weight/len(self.rfV1))**4
                    sig[ih, ik, ip] = np.var(amp)*np.real(weight)
                    pws[ih, ik, ip] = np.median(amp)*np.real(weight)

        self.pws = pws
        self.sig = sig

    def stack_dip(self, vp=None):
        """
        Function to calculate HK stacks from radial receiver functions
        """

        # P-wave velocity
        if not vp:
            try:
                vp = self.rfV1[0].stats.vp
            except:
                vp = 6.0

        sta = self.rfV1[0].stats.station

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        # Initialize arrays
        weight = np.complex(0.)
        amp = np.zeros(len(self.rfV1))
        sig = np.zeros((len(H), len(k), len(self.phases)))
        pws = np.zeros((len(H), len(k), len(self.phases)))

        for ih in _progressbar(range(len(H)), 'Computing: ', 15):
            for ik, kk in enumerate(k):
                for ip, ph in enumerate(self.phases):
                    for i in range(len(self.rfV1)):

                        if self.rfV2 and (ph=='pps' or ph=='pss'):
                            rfV = self.rfV2[i].copy()
                        else:
                            rfV = self.rfV1[i].copy()

                        # Calculate move out for each phase and get
                        # median value, weighted by instantaneous phase (pws)
                        tt = _dtime_dip_(
                            rfV, H[ih], kk, vp, ph, self.strike, self.dip)
                        trace = _timeshift_(rfV, tt)
                        thilb = hilbert(trace)
                        tphase = np.arctan2(thilb.imag, thilb.real)
                        weight += np.exp(1j*tphase[0])
                        amp[i] = trace[0]

                    weight = abs(weight/np.float(len(self.rfV1)))**4
                    sig[ih, ik, ip] = np.var(amp)*np.real(weight)
                    pws[ih, ik, ip] = np.median(amp)*np.real(weight)

        self.pws = pws
        self.sig = sig

    def average(self, typ='sum'):
        """
        # Function to average stacks from weights
        """

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        # Multiply pws by weights
        ps = self.pws[:, :, 0]*self.weights[0]
        try:
            pps = self.pws[:, :, 1]*self.weights[1]
        except:
            pps = None
        try:
            pss = self.pws[:, :, 2]*self.weights[2]
        except:
            pss = None

        # Get stacks
        if typ == 'sum':
            stack = (ps + pps + pss)
        elif typ == 'mult':
            # Zero out negative values
            ps[ps < 0] = 0.
            try:
                pps[pps < 0] = 0.
            except:
                pps = 1.
            try:
                pss[pss < 0] = 0.
            except:
                pss = 1.
            stack = ps*pps*pss
        else:
            raise(Exception("'typ' must be either 'sum' or 'mult'"))

        self.typ = typ

        # Find maximum within stacks
        ind = np.where(stack == stack.max())

        self.h0 = H[ind[0]][0]
        self.k0 = k[ind[1]][0]
        self.stack = stack


    def error(self, q=0.05, method='amp'):
        """
        Function to determine the q(%) confidence interval of the misfit
        function.

        From Walsh, JGR, 2013

        """

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        msf = self.stack/self.stack.max()

        # Method 1 - based on stats
        if method == 'stats':

            # Get degrees of freedom
            dof = _dof(self._residuals())
            print(dof)
            if dof < 3:
                dof = 3
                print(
                    "Degrees of freedom < 3. Fixing to DOF = 3, which may " +
                    "result in accurate errors")

            n_par = 2
            msf = 1. - msf

            # Error contour
            vmin = msf.min()
            vmax = msf.max()
            self.err_contour = vmin*(1. + n_par/(dof - n_par)*
                           stats.f.ppf(1. - q, n_par, dof - n_par))
            print(vmin*(1. + n_par/(dof - n_par)*
                           stats.f.ppf(1. - q, n_par, dof - n_par)))
            # self.err_contour = (n_par/(dof - n_par) *
            err = np.where(msf < self.err_contour)

        # Method 2 - based on amplitude
        elif method == 'amp':

            self.err_contour = 0.5
            err = np.where(msf > self.err_contour)

        self.method = method

        # Estimate uncertainty (q confidence interval)
        self.err_k = max(0.25*(k[max(err[1])] - k[min(err[1])]), self.dk)
        self.err_H = max(0.25*(H[max(err[0])] - H[min(err[0])]), self.dh)


    def plot(self, save=False, title=None):
        """
        Function to plot H-K stacks
        """

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        # Extent of plots
        extent = (H.min(), H.max(), k.min(), k.max())

        # Multiply pws by weights
        ps = self.pws[:, :, 0]*self.weights[0]
        try:
            pps = self.pws[:, :, 1]*self.weights[1]
        except:
            pps = None
        try:
            pss = self.pws[:, :, 2]*self.weights[2]
        except:
            pss = None

        if self.typ == 'mult':
            # Zero out negative values
            ps[ps < 0] = 0.
            try:
                pps[pps < 0] = 0.
            except:
                pass
            try:
                pss[pss < 0] = 0.
            except:
                pass

        # Set up figure
        # plt.clf()
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
            2, 2, sharex=True, sharey=True)

        cmap = 'RdBu_r'

        # First subplot: Ps
        #vmin = -np.abs(pws[:,:,0].max())
        # vmin=0.
        vmax = np.abs(max(ps.max(), ps.min(), key=abs))
        im = ax1.imshow(np.rot90(ps), cmap=cmap,
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax1.set_ylabel('Vp/Vs')
        ax1.set_title('Ps')

        # Second subplot: Pps
        #vmin = -np.abs(pws[:,:,1].max())
        # vmin=0.
        #vmax = np.abs(pps.max())
        vmax = np.abs(max(pps.max(), pps.min(), key=abs))
        im = ax2.imshow(np.rot90(pps), cmap=cmap,
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax2.set_title('Pps')

        # Third subplot: Pss
        #vmin = -np.abs(pws[:,:,2].max())
        # vmin=0.
        #vmax = np.abs(pss.max())
        vmax = np.abs(max(pss.max(), pss.min(), key=abs))
        im = ax3.imshow(np.rot90(pss), cmap=cmap,
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax3.set_title('Pss')
        ax3.set_ylabel('Vp/Vs')
        ax3.set_xlabel('Thickness (km)')

        # Fourth subplot: Average
        # vmin=0.
        #vmax = np.abs(stack.max())
        vmax = np.abs(max(self.stack.max(), self.stack.min(), key=abs))
        im = ax4.imshow(np.rot90(self.stack), cmap=cmap,
                        # extent=extent, aspect='auto')
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax4.set_title('Stack')
        ax4.set_xlabel('Thickness (km)')

        #cbar = fig.colorbar(im, ticks=[-vmax, 0, vmax])
        #cbar.ax.set_yticklabels(['min', '0', 'max'])

        # Get confidence intervals
        if hasattr(self, 'err_contour'):
            # ax.contour(np.rot90(vmax-msf), (vmax-err_cont,),
            if self.method == 'stats':
                ax4.contour(
                    np.rot90(1.-self.stack/self.stack.max()), 
                    (self.err_contour,),
                    hold='on', colors='yellow', linewidths=1, origin='upper',
                    extent=extent)
            elif self.method == 'amp':
                ax4.contour(
                    np.rot90(self.stack/self.stack.max()), 
                    (self.err_contour,),
                    hold='on', colors='yellow', linewidths=1, origin='upper',
                    extent=extent)

        # Add star showing best fit
        try:
            ax4.scatter(self.h0, self.k0, 60, marker='*', color='white')
        except:
            print("'h0' and 'k0' are not available")

        plt.suptitle('H-k stacks, station: ' + self.rfV1[0].stats.station)

        if save:
            plt.savefig('RF_PLOTS/' + self.rfV1[0].stats.station +
                        '.' + title+'.eps', format='eps')
            plt.close()
        else:
            plt.show()

    def _residuals(self):
        """ 
        Internal method to obtain residuals between observed and predicted
        receiver functions given the Moho depth and Vp/Vs obtainedw from
        the Hk stack.
        """
        from telewavesim import utils

        # Simple 1-layer model over half-space
        model = utils.Model(
            [self.h0, 0.], 
            [2800., 3300.], 
            [self.vp, 8.0],
            [self.vp/self.k0, 4.5],
            ['iso', 'iso'])

        # Parameters for run
        slow = [tr.stats.slow for tr in self.rfV1]
        npts = self.rfV1[0].stats.npts
        dt = self.rfV1[0].stats.delta

        trR = Stream()

        for sl in slow:
            trxyz = utils.run_plane(model, sl, npts, dt)
            tfs = utils.tf_from_xyz(trxyz, pvh=True, vp=self.vp, vs=self.vp/self.k0)
            tfs[0].data = np.fft.fftshift(tfs[0].data)
            trR.append(tfs[0])

        trR.filter('bandpass', freqmin=0.05, freqmax=0.5, corners=2,
            zerophase=True)

        # Get stream of residuals
        res = trR.copy()
        for i in range(len(res)):
            res[i].data = self.rfV1[i].data - trR[i].data
        return res


def _dof(st):
    """
    Function to determine the degrees of freedom to calculate
    the confidence region of the misfit function.

    From Walsh, JGR, 2013

    """

    dof = []
    for tr in st:

        F = np.abs(np.fft.fft(tr.data)[0:int(len(tr.data)/2)+1])

        E2 = np.sum(F**2)
        E2 -= (F[0]**2 + F[-1]**2)/2.
        E4 = (1./3.)*(F[0]**4 + F[-1]**4)
        for i in range(1, len(F) - 1):
            E4 += (4./3.)*F[i]**4

        dof.append(int(4.*E2**2/E4 - 2.))

    dof_max = min(dof)

    return dof_max


# Function to calculate travel time for different scattered phases
def _dtime_(trace, z, r, vp, ph):

    # Horizontal slowness
    slow = trace.stats.slow

    # Vertical slownesses
    c1 = np.sqrt((r/vp)**2 - slow**2)
    c2 = np.sqrt((1./vp)**2 - slow**2)

    if ph == 'ps':
        tt = z*(c1 - c2)
    elif ph == 'pps':
        tt = z*(c1 + c2)
    elif ph == 'pss':
        tt = 2.*z*c1

    return tt

# Function to calculate travel time for different scattered phases


def _dtime_dip_(trace, z, r, ph, strike, dip):

    vp = trace.stats.vp
    n = np.zeros(3)
    pinc = np.zeros(3)
    ai = 8.1
    br = vp/r

    dip = dip*np.pi/180.
    strike = strike*np.pi/180.
    n[0] = np.sin(strike)*np.sin(dip)
    n[1] = -np.cos(strike)*np.sin(dip)
    n[2] = np.cos(dip)

    # Horizontal slowness
    slow = trace.stats.slow
    baz = trace.stats.baz

    # Assemble constants
    c1 = n[2]
    theta = baz*np.pi/180.+np.pi
    pinc[0] = slow*np.cos(theta)
    pinc[1] = slow*np.sin(theta)
    pinc[2] = -np.sqrt(1./(ai*ai) - slow*slow)

    ndotpi = n[0]*pinc[0] + n[1]*pinc[1] + n[2]*pinc[2]
    c2 = 1./(ai*ai) - ndotpi**2

    if ph == 'ps':
        a = vp*vp
        tt = z*(c1*(np.sqrt(r*r/a - c2) - np.sqrt(1./a - c2)))

    elif ph == 'pps':
        pref = pinc - ndotpi*n - np.sqrt(1./(vp*vp)-c2)*n
        pref[2] = -pref[2]
        ndotpr = n[0]*pref[0]+n[1]*pref[1]+n[2]*pref[2]
        c4 = 1./(vp*vp) - ndotpr**2
        c3 = 2.*pref[2] - 2.*n[2]*ndotpr
        a = vp*vp
        tt = z*(c1*(np.sqrt(r*r/a - c4) + np.sqrt(1./a - c4)))

    elif ph == 'pss':
        tt = 2.*z*c1
        pref = pinc - ndotpi*n - np.sqrt(1./(vp*vp) - c2)*n
        pref[2] = np.sqrt(1./(br*br) - pref[0]*pref[0] - pref[1]*pref[1])
        ndotpr = n[0]*pref[0] + n[1]*pref[1] + n[2]*pref[2]
        c6 = 1./(br*br) - ndotpr**2
        a = vp*vp
        tt = z*(2.*c1*np.sqrt(r*r/a - c6))

    return tt

# Function to shift traces in time given travel time


def _timeshift_(trace, tt):

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
    rtrace = np.real(np.fft.ifft(ftrace))

    return rtrace


def _progressbar(it, prefix="", size=60, file=sys.stdout):
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
