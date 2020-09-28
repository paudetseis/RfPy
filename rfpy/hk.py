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
from scipy.signal import hilbert
from scipy import stats
import sys
from matplotlib import pyplot as plt


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

    Other Parameters
    ----------------
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

    """

    def __init__(self, rfV1, rfV2=None, strike=None, dip=None, vp=6.0):

        # Load example data if initializing empty object
        if rfV1 == 'demo' or rfV1 == 'Demo':
            print("Uploading demo data - station NY.MMPY")
            import os
            import pickle
            file = open(os.path.join(
                os.path.dirname(__file__),
                "examples/data", "demo_streams.pkl"), 'rb')
            rfV1 = pickle.load(file)
            file.close()

        # fftshift if the time axis starts at negative lags
        if rfV1[0].stats.taxis[0] < 0.:
            nn = rfV1[0].stats.npts
            for tr in rfV1:
                tr.data = np.fft.fftshift(tr.data)[0:int(nn/2)]
                tr.stats.taxis = np.fft.fftshift(tr.data)[0:int(nn/2)]
            try:
                for tr in rfV2:
                    tr.data = np.fft.fftshift(tr.data)[0:int(nn/2)]
                    tr.stats.taxis = np.fft.fftshift(tr.data)[0:int(nn/2)]
            except:
                pass

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
        lower frequency signals). Furthermore, If the ``vp`` argument is 
        not specified, the method will use the 
        value set during initialization (``vp=6.0`` km/s)

        Parameters
        ----------
        vp : float
            Mean crust P-wave velocity (km/s). 

        Attributes
        ----------
        pws : :class:`~numpy.ndarray`
            Array of phase stacks, where the outer dimension corresponds
            to the phase index (shape ``nH, nk, nph``)
        sig : :class:`~numpy.ndarray`
            Variance of phase stacks, where the outer dimension corresponds
            to the phase index (shape ``nH, nk, nph``)

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

                        if self.rfV2 and (ph == 'pps' or ph == 'pss'):
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

                        # ### Attempt at speeding things up
                        # ind = (np.abs(rfV.stats.taxis - tt)).argmin()
                        # trace = rfV.copy()
                        # thilb = hilbert(trace.data)
                        # tphase = np.arctan2(thilb.imag, thilb.real)
                        # weight += np.exp(1j*tphase[ind])
                        # amp[i] = trace.data[ind]

                    weight = abs(weight/len(self.rfV1))**4
                    sig[ih, ik, ip] = np.var(amp)*np.real(weight)
                    pws[ih, ik, ip] = np.median(amp)*np.real(weight)

        self.pws = pws
        self.sig = sig

    def stack_dip(self, vp=None, strike=None, dip=None):
        """
        Method to calculate Hk stacks from radial receiver functions
        using known stike and dip angles of the Moho.
        The stacks are calculated using phase-weighted stacking for
        individual phases and take the median of the weighted stack 
        to avoid bias by outliers.

        Note
        ----
        If two streams are available as attributes, the method will assume
        that the second stream should be used for stacking along the Pps
        and Pss move out curves (e.g., if the second stream contains
        lower frequency signals). Furthermore, 
        If the arguments are not specified, the method will use the 
        values set during initialization (``vp=6.0`` km/s, 
        ``strike=0.``, ``dip=0.``)

        Parameters
        ----------
        vp : float
            Mean crust P-wave velocity (km/s). 
        strike : float
            Strike angle of dipping Moho (has to be known or estimated a priori)
        dip : float
            Dip angle of Moho (has to be known or estimated a priori)

        Attributes
        ----------
        pws : :class:`~numpy.ndarray`
            Array of phase stacks, where the outer dimension corresponds
            to the phase index (shape ``nH, nk, nph``)
        sig : :class:`~numpy.ndarray`
            Variance of phase stacks, where the outer dimension corresponds
            to the phase index (shape ``nH, nk, nph``)

        """

        if not strike:
            strike = self.strike
        else:
            self.strike = strike
        if not dip:
            dip = self.dip
        else:
            self.dip = dip

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

                        if self.rfV2 and (ph == 'pps' or ph == 'pss'):
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

    def average(self, typ='sum', q=0.05, err_method='amp'):
        """
        Method to combine the phase-weighted stacks to produce a final
        stack, from which to estimate the H and k parameters and their 
        associated errors.

        Parameters
        ----------
        typ : str
            How the phase-weigthed stacks should be combined to produce
            a final stack. Available options are: weighted sum (``typ=sum``) 
            or product (``typ=product``).
        q : float
            Confidence level for the error estimate
        err_method : str
            How errors should be estimated. Options are ``err_method='amp'``
            to estimate errors from amplitude, or ``err_method='stats'`` to 
            use a statistical F test from the residuals.

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
        elif typ == 'product':
            # Zero out negative values
            ps[ps < 0] = 0.
            if self.weights[1] != 0.:
                pps[pps < 0] = 0.
            else:
                pps = 1.
            if self.weights[2] != 0.:
                pss[pss < 0] = 0.
            else:
                pss = 1.
            stack = ps*pps*pss
        else:
            raise(Exception("'typ' must be either 'sum' or 'product'"))

        self.typ = typ

        # Find maximum within stacks
        ind = np.where(stack == stack.max())

        self.h0 = H[ind[0]][0]
        self.k0 = k[ind[1]][0]
        self.stack = stack

        try:
            self.error()
        except:
            self.err_k0 = 0.
            self.err_h0 = 0.

    def error(self, q=0.05, err_method='amp'):
        """
        Method to determine the error on H and k estimates.

        From Walsh, JGR, 2013

        Parameters
        ----------
        q : float
            Confidence level for the error estimate
        err_method : str
            How errors should be estimated. Options are ``err_method='amp'``
            to estimate errors from amplitude, or ``err_method='stats'`` to 
            use a statistical F test from the residuals.

        """

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        msf = self.stack/self.stack.max()

        # Method 1 - based on stats
        if err_method == 'stats':

            # Get degrees of freedom
            dof = _dof(self._residuals())
            # print(dof)
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
            self.err_contour = vmin*(1. + n_par/(dof - n_par) *
                                     stats.f.ppf(1. - q, n_par, dof - n_par))
            # print(vmin*(1. + n_par/(dof - n_par)*
            #                stats.f.ppf(1. - q, n_par, dof - n_par)))
            # self.err_contour = (n_par/(dof - n_par) *
            err = np.where(msf < self.err_contour)

        # Method 2 - based on amplitude
        elif err_method == 'amp':

            self.err_contour = 0.5
            err = np.where(msf > self.err_contour)

        else:
            raise(Exception("'err_method' must be either 'stats' or 'amp'"))
        self.err_method = err_method

        # Estimate uncertainty (q confidence interval)
        self.err_k0 = max(0.25*(k[max(err[1])] - k[min(err[1])]), self.dk)
        self.err_h0 = max(0.25*(H[max(err[0])] - H[min(err[0])]), self.dh)

    def plot(self, save=False, title=None, form='png'):
        """
        Method to plot H-K stacks. By default all 4 panels
        are plotted: The ``ps``, ``pps`` and ``pss`` stacks, and the
        final (averaged) stack. Error contours are also plotted,
        as well as the position of the maximum stack values.

        Parameters
        ----------
        save : bool
            Whether or not to save the Figure
        title : str
            Title of plot
        """

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        # Extent of plots
        extent = (H.min(), H.max(), k.min(), k.max())

        # Extract phase stacks
        ps = self.pws[:, :, 0]
        pps = self.pws[:, :, 1]
        pss = self.pws[:, :, 2]

        if self.typ == 'product':
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
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
            2, 2, sharex=True, sharey=True)

        cmap = 'RdBu_r'

        # First subplot: Ps
        vmax = np.abs(max(ps.max(), ps.min(), key=abs))
        im = ax1.imshow(np.rot90(ps), cmap=cmap,
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax1.set_ylabel('Vp/Vs')
        ax1.set_title('Ps - weight: {0:.1f}'.format(
            self.weights[0]), fontsize=10)

        # Second subplot: Pps
        vmax = np.abs(max(pps.max(), pps.min(), key=abs))
        im = ax2.imshow(np.rot90(pps), cmap=cmap,
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax2.set_title('Pps - weight: {0:.1f}'.format(
            self.weights[1]), fontsize=10)

        # Third subplot: Pss
        vmax = np.abs(max(pss.max(), pss.min(), key=abs))
        im = ax3.imshow(np.rot90(pss), cmap=cmap,
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax3.set_title('Pss - weight: {0:.1f}'.format(
            self.weights[2]), fontsize=10)
        ax3.set_ylabel('Vp/Vs')
        ax3.set_xlabel('Thickness (km)')

        # Fourth subplot: Average
        vmax = np.abs(max(self.stack.max(), self.stack.min(), key=abs))
        im = ax4.imshow(np.rot90(self.stack), cmap=cmap,
                        extent=extent, vmin=-vmax, vmax=vmax, aspect='auto')
        ax4.set_title('Stack')
        ax4.set_xlabel('Thickness (km)')

        #cbar = fig.colorbar(im, ticks=[-vmax, 0, vmax])
        #cbar.ax.set_yticklabels(['min', '0', 'max'])

        # Get confidence intervals
        if hasattr(self, 'err_contour'):
            # ax.contour(np.rot90(vmax-msf), (vmax-err_cont,),
            if self.err_method == 'stats':
                ax4.contour(
                    np.rot90(1.-self.stack/self.stack.max()),
                    (self.err_contour,),
                    hold='on', colors='yellow', linewidths=1, origin='upper',
                    extent=extent)
            elif self.err_method == 'amp':
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

        if title:
            plt.suptitle(title)
        else:
            plt.suptitle('H-k stacks, station: ' + self.rfV1[0].stats.station)

        if save:
            plt.savefig('HK_PLOTS/hk.' + self.rfV1[0].stats.station +
                        '.' + title+'.'+self.typ+'.'+form, format=form)
        else:
            plt.show()

        plt.close()

## JMG ##
    def save(self, file):
        ## JMG ##
        """
        Saves HkStack object to file

        Parameters
        ----------
        file : str
            File name for HkStack object

        """

        import pickle
        output = open(file, 'wb')
        pickle.dump(self, output)
        output.close()

    def _residuals(self):
        """ 
        Internal method to obtain residuals between observed and predicted
        receiver functions given the Moho depth and Vp/Vs obtained from
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
            tfs = utils.tf_from_xyz(
                trxyz, pvh=True, vp=self.vp, vs=self.vp/self.k0)
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
    Method to determine the degrees of freedom to calculate
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


def _dtime_(trace, z, r, vp, ph):
    """
    Method to calculate travel time for different scattered phases

    """

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


def _dtime_dip_(trace, z, r, vp, ph, strike, dip):
    """
    Method to calculate travel time for different scattered phases
    using strike and dip angles

    """

    # Initialize some parameters
    n = np.zeros(3)
    pinc = np.zeros(3)
    ai = 8.1
    br = vp/r

    # Get vector normal to dipping interface
    dip = dip*np.pi/180.
    strike = strike*np.pi/180.
    n[0] = np.sin(strike)*np.sin(dip)
    n[1] = -np.cos(strike)*np.sin(dip)
    n[2] = np.cos(dip)

    # Horizontal slowness
    slow = trace.stats.slow
    # Back-azimuth
    baz = trace.stats.baz

    # Assemble constants of incident wave
    c1 = n[2]
    theta = baz*np.pi/180.+np.pi
    pinc[0] = slow*np.cos(theta)
    pinc[1] = slow*np.sin(theta)
    pinc[2] = -np.sqrt(1./(ai*ai) - slow*slow)

    # Calculate scalar product n * pinc
    ndotpi = n[0]*pinc[0] + n[1]*pinc[1] + n[2]*pinc[2]
    # Incident normal slowness
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


def _timeshift_(trace, tt):
    """
    Function to shift traces in time given travel time

    """

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
