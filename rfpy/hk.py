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


class HkStack(object):

    def __init__(self, rfV1, rfV2=None, strike=None, dip=None):

        self.rfV1 = rfV1
        self.rfV2 = rfV2
        self.kbound = [1.56, 2.1]
        self.dk = 0.01
        self.hbound = [20., 50.]
        self.dh = 0.25
        self.weights = [0.5, 2., -1.]
        self.phases = ['ps', 'pps', 'pss']
        self.strike = strike
        self.dip = dip

    def stack(self, vp=None):
        """
        Method to calculate HK stacks from radial receiver functions.

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

        for ih in range(len(H)):

            print()
            print('Thickness(ih)', H[ih], ', ih = ', ih, '/', len(H)-1)

            for ik in range(len(k)):

                for ip in range(len(self.phases)):

                    for i in range(len(self.rfV1)):

                        if self.rfV2 and (self.phases[ip]=='pps' or
                                          self.phases[ip]=='pss'):
                            rfV = self.rfV2[i].copy()
                        else:
                            rfV = self.rfV1[i].copy()

                        # Calculate move out for each phase and get
                        # median value, weighted by instantaneous phase (pws)
                        tt = _dtime_(rfV, H[ih], k[ik], vp, self.phases[ip])
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

        for ih in range(len(H)):

            print()
            print('Thickness(ih)', H[ih], ', ih = ', ih, '/', len(H)-1)

            for ik in range(len(k)):

                for ip in range(len(self.phases)):

                    for i in range(len(self.rfV1)):

                        if self.rfV2 and (self.phases[ip]=='pps' or
                                          self.phases[ip]=='pss'):
                            rfV = self.rfV2[i].copy()
                        else:
                            rfV = self.rfV1[i].copy()

                        # Calculate move out for each phase and get
                        # median value, weighted by instantaneous phase (pws)
                        tt = _dtime_dip_(
                            rfV, H[ih], k[ik], vp, 
                            self.phases[ip], self.strike, self.dip)
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
        # return H[ind[0]][0], k[ind[1]][0], stack  # , stack.max()

    def error(res, q, method='amp'):
        """
        Function to determine the q(%) confidence interval of the misfit
        function.

        From Walsh, JGR, 2013

        """

        # Initialize arrays based on bounds
        H = np.arange(self.hbound[0], self.hbound[1] + self.dh, self.dh)
        k = np.arange(self.kbound[0], self.kbound[1] + self.dk, self.dk)

        # Get degrees of freedom
        dof = _dof(res)
        n_par = 2

        msf = self.stack/self.stack.max()

        # Method 1 - based on stats
        if method == 'stats':

            msf = 1. - msf

            # Error contour
            vmin = msf.min()
            vmax = msf.max()
            # err_contour = vmin*(1. + n_par/(dof - n_par)*
            err_contour = (n_par/(dof - n_par) *
                           stats.f.ppf(1.-q, n_par, dof-n_par))
            err = np.where(msf < err_contour)

        # Method 2 - based on amplitude
        elif method == 'amp':

            err_contour = 0.5
            err = np.where(msf > err_contour)

        self.method = method

        # Estimate uncertainty (q confidence interval)
        err_k = 0.25*(k[max(err[1])] - k[min(err[1])])
        err_H = 0.25*(H[max(err[0])] - H[min(err[0])])

        return err_k, err_H, err_contour

    def plot(self, err_c=None, save=False, title=None):
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
        if err_c:
            # ax.contour(np.rot90(vmax-msf), (vmax-err_cont,),
            if self.method == 'stats':
                ax4.contour(
                    np.rot90(1.-self.stack/self.stack.max()), (err_c,),
                    hold='on', colors='yellow', linewidths=1, origin='upper',
                    extent=extent)
            elif self.method == 'amp':
                ax4.contour(
                    np.rot90(self.stack/self.stack.max()), (err_c,),
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


def _dof(st):
    """
    Function to determine the degrees of freedom to calculate
    the confidence region of the misfit function.

    From Walsh, JGR, 2013

    """

    dof = []
    for tr in st:

        ft = np.fft.fft(tr.data)[0:len(tr.data)/2+1]
        ai = np.ones(len(ft))
        ai[0] = 0.5
        ai[-1] = 0.5

        E2 = np.sum(np.dot(ai, np.abs(ft)**2))
        E4 = np.sum(np.dot(ai, np.abs(ft)**4))

        dof.append(2.*(2.*E2**2/E4 - 1.))

    dof_max = max(dof)

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
