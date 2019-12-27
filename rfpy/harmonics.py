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
Harmonic decomposition module.

"""

# Import modules and functions
import numpy as np
from obspy.core import Stream, Trace
from matplotlib import pyplot as plt


class Harmonics(object):

    def __init__(self, strV, strH, azim=0, xmin=0., xmax=40.):
        self.strV = strV
        self.strH = strH
        self.azim = azim
        self.xmin = xmin
        self.xmax = xmax

    def dcomp_find_azim(self, xmin=None, xmax=None):
        """
        Method to decompose radial and transverse receiver function 
        streams into back-azimuth harmonics and determine the main 
        orientation ``azim``, obtained by minimizing the B1 component
        between ``xmin`` and ``xmax`` (i.e., time or depth).

        Parameters
        ----------
        xmin : float
            Minimum x axis value over which to calculate ``azim``
        xmax : float
            Maximum x axis value over which to calculate ``azim``

        Attributes
        ----------
        hstream : :class:`~obspy.core.Stream`
            Stream containing the 5 harmonics, oriented in direction ``azim``
        azim : float
            Direction (azimuth) along which the B1 component of the stream
            is minimized (between ``xmin`` and ``xmax``)
        var : :class:`~numpy.ndarray`
            Variance of the 5 harmonics between ``xmin`` and ``xmax``

        """

        if not xmin:
            xmin = self.xmin
        if not xmax:
            xmax = self.xmax

        print()
        print('Decomposing receiver functions into baz harmonics')

        # Some integers
        nbin = len(self.strV)
        nz = len(self.strV[0].data)
        naz = 180
        daz = np.float(360/naz)
        deg2rad = np.pi/180.

        # Define depth range over which to calculate azimuth
        indmin = int(xmin/self.strV[0].stats.delta)
        indmax = int(xmax/self.strV[0].stats.delta)

        # Copy stream stats
        str_stats = self.strV[0].stats

        # Initialize work arrays
        C0 = np.zeros((nz, naz))
        C1 = np.zeros((nz, naz))
        C2 = np.zeros((nz, naz))
        C3 = np.zeros((nz, naz))
        C4 = np.zeros((nz, naz))

        # Loop over each depth step
        for iz in range(nz):

            # Build matrices OBS and H for each azimuth
            for iaz in range(naz):

                # Initialize work arrays
                OBS = np.zeros(2*nbin)
                H = np.zeros((2*nbin, 5))

                azim = iaz*daz

                # Radial component
                for irow, trace in enumerate(self.strV):

                    baz = trace.stats.baz
                    OBS[irow] = trace.data[iz]
                    H[irow, 0] = 1.0
                    H[irow, 1] = np.cos(deg2rad*(baz-azim))
                    H[irow, 2] = np.sin(deg2rad*(baz-azim))
                    H[irow, 3] = np.cos(2.*deg2rad*(baz-azim))
                    H[irow, 4] = np.sin(2.*deg2rad*(baz-azim))

                shift = 90.

                # Transverse component
                for irow, trace in enumerate(self.strH):

                    baz = trace.stats.baz
                    OBS[irow+nbin] = trace.data[iz]
                    H[irow+nbin, 0] = 0.0
                    H[irow+nbin, 1] = np.cos(deg2rad*(baz+shift-azim))
                    H[irow+nbin, 2] = np.sin(deg2rad*(baz+shift-azim))
                    H[irow+nbin, 3] = np.cos(2.*deg2rad*(baz+shift/2.0-azim))
                    H[irow+nbin, 4] = np.sin(2.*deg2rad*(baz+shift/2.0-azim))

                # Solve system of equations with truncated SVD
                u, s, v = np.linalg.svd(H)
                s[s < 0.001] = 0.
                CC = np.linalg.solve(s[:, None] * v, u.T.dot(OBS)[:5])

                # Fill up arrays
                C0[iz, iaz] = np.float(CC[0])
                C1[iz, iaz] = np.float(CC[1])
                C2[iz, iaz] = np.float(CC[2])
                C3[iz, iaz] = np.float(CC[3])
                C4[iz, iaz] = np.float(CC[4])

        # Minimize variance of third component over specific depth range to
        # find azim
        C1var = np.zeros(naz)
        for iaz in range(naz):
            C1var[iaz] = np.sqrt(np.mean(np.square(C1[indmin:indmax, iaz])))
        indaz = np.argmin(C1var)

        C0var = np.sqrt(np.mean(np.square(C0[indmin:indmax, indaz])))
        C1var = np.sqrt(np.mean(np.square(C1[indmin:indmax, indaz])))
        C2var = np.sqrt(np.mean(np.square(C2[indmin:indmax, indaz])))
        C3var = np.sqrt(np.mean(np.square(C3[indmin:indmax, indaz])))
        C4var = np.sqrt(np.mean(np.square(C4[indmin:indmax, indaz])))

        # Put back into traces
        A = Trace(data=C0[:, indaz], header=str_stats)
        B1 = Trace(data=C1[:, indaz], header=str_stats)
        B2 = Trace(data=C2[:, indaz], header=str_stats)
        C1 = Trace(data=C3[:, indaz], header=str_stats)
        C2 = Trace(data=C4[:, indaz], header=str_stats)

        # Put all treaces into stream
        self.hstream = Stream(traces=[A, B1, B2, C1, C2])
        self.azim = indaz*daz
        self.var = [C0var, C1var, C2var, C3var, C4var]

    def dcomp_fix_azim(self, azim=None):
        """
        Method to decompose radial and transverse receiver function 
        streams into back-azimuth harmonics along direction ``azim``.

        Parameters
        ----------
        azim : float
            Direction (azimuth) along which the B1 component of the stream
            is minimized (between ``xmin`` and ``xmax``)

        Attributes
        ----------
        hstream : :class:`~obspy.core.Stream`
            Stream containing the 5 harmonics, oriented in direction ``azim``

        """

        if azim is None:
            azim = self.azim
        else:
            self.azim = azim

        print('Decomposing receiver functions into baz harmonics for az = ',
              azim)

        # Some integers
        nbin = len(self.strV)
        nz = len(self.strV[0].data)
        deg2rad = np.pi/180.

        # Copy stream stats
        str_stats = self.strV[0].stats

        # Initialize work arrays
        C0 = np.zeros(nz)
        C1 = np.zeros(nz)
        C2 = np.zeros(nz)
        C3 = np.zeros(nz)
        C4 = np.zeros(nz)

        # Loop over each depth step
        for iz in range(nz):

            # Initialize working arrays
            OBS = np.zeros(2*nbin)
            H = np.zeros((2*nbin, 5))

            # Radial component
            for irow, trace in enumerate(self.strV):

                baz = trace.stats.baz
                OBS[irow] = trace.data[iz]
                H[irow, 0] = 1.0
                H[irow, 1] = np.cos(deg2rad*(baz-azim))
                H[irow, 2] = np.sin(deg2rad*(baz-azim))
                H[irow, 3] = np.cos(2.*deg2rad*(baz-azim))
                H[irow, 4] = np.sin(2.*deg2rad*(baz-azim))

            shift = 90.

            # Transverse component
            for irow, trace in enumerate(self.strH):

                baz = trace.stats.baz
                OBS[irow+nbin] = trace.data[iz]
                H[irow+nbin, 0] = 0.0
                H[irow+nbin, 1] = np.cos(deg2rad*(baz+shift-azim))
                H[irow+nbin, 2] = np.sin(deg2rad*(baz+shift-azim))
                H[irow+nbin, 3] = np.cos(2.*deg2rad*(baz+shift/2.0-azim))
                H[irow+nbin, 4] = np.sin(2.*deg2rad*(baz+shift/2.0-azim))

            # Solve system of equations with truncated SVD
            u, s, v = np.linalg.svd(H)
            s[s < 0.001] = 0.
            CC = np.linalg.solve(s[:, None] * v, u.T.dot(OBS)[:5])

            # Fill up arrays
            C0[iz] = np.float(CC[0])
            C1[iz] = np.float(CC[1])
            C2[iz] = np.float(CC[2])
            C3[iz] = np.float(CC[3])
            C4[iz] = np.float(CC[4])

        # Put back into traces
        A = Trace(data=C0, header=str_stats)
        B1 = Trace(data=C1, header=str_stats)
        B2 = Trace(data=C2, header=str_stats)
        C1 = Trace(data=C3, header=str_stats)
        C2 = Trace(data=C4, header=str_stats)

        # Put all traces into stream
        self.hstream = Stream(traces=[A, B1, B2, C1, C2])

    def forward(self, baz_list=None):
        """
        Method to forward calculate radial and transverse component
        receiver functions given the 5 harmonics and a list of back-azimuth
        values. The receiver function signal parameters (length, sampling
        rate, etc.) will be identical to those in the stream of
        harmonic components.

        Parameters
        ----------
        baz_list : list
            List of back-azimuth directions over which to calculate
            the receiver functions. If no list is specified, the method
            will use the same back-azimuths as those in the original
            receiver function streams

        Attributes
        ----------
        forwardV : :class:`~obspy.core.Stream`
            Stream containing the radial receiver functions
        forwardH : :class:`~obspy.core.Stream`
            Stream containing the transverse receiver functions


        """

        if not hasattr(self, 'hstream'):
            raise(Exception("Decomposition has not been performed yet"))

        if not baz_list:
            print("Warning: no BAZ specified - using all baz from " +
                  "stored streams")
            baz_list = [tr.stats.baz for tr in self.strV]
        if not isinstance(baz_list, list):
            baz_list = [baz_list]

        # Some constants
        nz = len(self.hstream[0].data)
        deg2rad = np.pi/180.

        # Copy traces
        forwardV = Stream()
        forwardH = Stream()

        for baz in baz_list:
            trV = Trace(header=self.hstream[0].stats)
            trH = Trace(header=self.hstream[0].stats)

            # Loop over each time/depth step
            for iz in range(nz):

                # Initialize working arrays
                X = np.zeros(5)
                H = np.zeros((2, 5))

                # Fill up X array
                X[0] = hstream[0].data[iz]
                X[1] = hstream[1].data[iz]
                X[2] = hstream[2].data[iz]
                X[3] = hstream[3].data[iz]
                X[4] = hstream[4].data[iz]

                # Fill up H arrays (for V and H)
                H[0, 0] = 1.0
                H[0, 1] = np.cos(deg2rad*(baz-self.azim))
                H[0, 2] = np.sin(deg2rad*(baz-self.azim))
                H[0, 3] = np.cos(2.*deg2rad*(baz-self.azim))
                H[0, 4] = np.sin(2.*deg2rad*(baz-self.azim))

                shift = 90.

                H[1, 0] = 0.0
                H[1, 1] = np.cos(deg2rad*(baz+shift-self.azim))
                H[1, 2] = np.sin(deg2rad*(baz+shift-self.azim))
                H[1, 3] = np.cos(2.*deg2rad*(baz+shift/2.0-self.azim))
                H[1, 4] = np.sin(2.*deg2rad*(baz+shift/2.0-self.azim))

                # Calculate dot product B = H*X
                B = np.dot(H, X)

                # Extract receiver functions
                trV.data[iz] = B[0]
                trH.data[iz] = -B[1]

        self.forwardR = forwardR
        self.forwardT = forwardT

    def plot(self, ymax=30., maxval=10., save=False, title=None):
        """
        Method to plot the 5 harmonic components.

        Parameters
        ----------
        ymax : float
            Maximum y axis value (time or depth) over which to 
            plot the harmonic components
        maxval : float
            Scaling factor for the amplitudes (typically > 1)
        save : bool
            Whether or not to save the plot
        title : str
            Title of plot, to be used in the Figure and the 
            file name (if ``save==True``)

        """

        # Y axis
        y = np.arange(self.hstream[0].stats.npts) /\
            self.hstream[0].stats.sampling_rate

        # Station name
        sta = self.hstream[0].stats.station

        # Initialize count
        i = 0

        # Initialize figure
        fig = plt.figure()
        plt.clf()

        # Get more control on subplots
        # ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.5])
        ax1 = fig.add_subplot(111)

        for i, trace in enumerate(self.hstream):
            i += 1
            ax1.fill_betweenx(
                y, i, i+trace.data*maxval,
                where=trace.data+1e-6 <= 0.,
                facecolor='blue',
                linewidth=0)
            ax1.fill_betweenx(
                y, i, i+trace.data*maxval,
                where=trace.data+1e-6 >= 0.,
                facecolor='red',
                linewidth=0)

        ax1.set_ylim(ymax, 0)
        ax1.set_ylabel('Depth (km)')
        ax1.set_xlabel('Harmonic components')
        if title:
            ax1.set_title('Station '+sta)
        labels = [item.get_text() for item in ax1.get_xticklabels()]
        labels[1] = '$A$'
        labels[2] = '$B_1$'
        labels[3] = '$B_2$'
        labels[4] = '$C_1$'
        labels[5] = '$C_2$'
        ax1.set_xticklabels(labels)
        off = ax1.xaxis.get_offset_text()
        ax1.tick_params(axis=u'x', pad=10)
        ax1.grid()

        if save:
            plt.savefig('RF_PLOTS/'+sta+'.'+title+'.eps', dpi=300,
                        bbox_inches='tight', format='eps')
        else:
            plt.show()

    # Plot wiggles for harmonic decomposition
