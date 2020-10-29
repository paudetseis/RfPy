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
import matplotlib.pyplot as plt


class Harmonics(object):
    """
    A Harmonics object contains attributes and methods to decompose
    radial and transverse component receiver functions into
    back-azimuth harmonics. The object is initialized with two
    :class:`~obspy.core.Stream` objects containing observed (or synthetised)
    radial and transverse receiver functions. The methods available
    can decompose the receiver functions along a fixed azimuth, or
    search for the optimal azimuth within a time range by minimizing
    one component.

    Note
    ----
    The object is initialized with the ``rfV1`` field only, and
    other attributes are added to the object as the analysis proceeds.
    A second ``rfV2`` can be included, which is typically a copy of ``rfV1``
    filtered at different corner frequencies and is used to stack along the
    Pps and Pss moveout curves.

    Parameters
    ----------
    radialRF : :class:`~obspy.core.Stream`
        Stream object containing the radial-component receiver function
        seismograms
    transvRF : :class:`~obspy.core.Stream`
        Stream object containing the transverse-component receiver function
        seismograms
    azim : float
        Direction (azimuth) along which the B1 component of the stream
        is minimized (between ``xmin`` and ``xmax``)
    xmin : float
        Minimum x axis value over which to calculate ``azim``
    xmax : float
        Maximum x axis value over which to calculate ``azim``

    Other Parameters
    ----------------
    hstream : :class:`~obspy.core.Stream`
        Stream containing the 5 harmonics, oriented in direction ``azim``
    radial_forward : :class:`~obspy.core.Stream`, optional
        Stream containing the radial receiver functions
    transv_forward : :class:`~obspy.core.Stream`, optional
        Stream containing the transverse receiver functions

    """

    def __init__(self, radialRF, transvRF=None, azim=0, xmin=0., xmax=10.):

        # Load example data if initializing empty object
        if radialRF == 'demo' or radialRF == 'Demo':
            print("Uploading demo data - station NY.MMPY")
            import os
            import pickle
            file = open(os.path.join(
                    os.path.dirname(__file__),
                    "examples/data", "demo_streams.pkl"), 'rb')
            radialRF = pickle.load(file)
            transvRF = pickle.load(file)
            file.close()

        if not transvRF:
            raise TypeError("__init__() missing 1 required positional argument: 'transvRF'")

        # fftshift if the time axis starts at negative lags
        if radialRF[0].stats.taxis[0]<0.:
            for tr in radialRF:
                tr.data = np.fft.fftshift(tr.data)
            for tr in transvRF:
                tr.data = np.fft.fftshift(tr.data)

        self.radialRF = radialRF
        self.transvRF = transvRF
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
        nbin = len(self.radialRF)
        nz = len(self.radialRF[0].data)
        naz = 180
        daz = np.float(360/naz)
        deg2rad = np.pi/180.

        # Define depth range over which to calculate azimuth
        indmin = int(xmin/self.radialRF[0].stats.delta)
        indmax = int(xmax/self.radialRF[0].stats.delta)

        # Copy stream stats
        str_stats = self.radialRF[0].stats

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
                for irow, trace in enumerate(self.radialRF):

                    baz = trace.stats.baz
                    OBS[irow] = trace.data[iz]
                    H[irow, 0] = 1.0
                    H[irow, 1] = np.cos(deg2rad*(baz-azim))
                    H[irow, 2] = np.sin(deg2rad*(baz-azim))
                    H[irow, 3] = np.cos(2.*deg2rad*(baz-azim))
                    H[irow, 4] = np.sin(2.*deg2rad*(baz-azim))

                shift = 90.

                # Transverse component
                for irow, trace in enumerate(self.transvRF):

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

        print('Decomposing receiver functions into baz harmonics for azimuth = ',
              azim)

        # Some integers
        nbin = len(self.radialRF)
        nz = len(self.radialRF[0].data)
        deg2rad = np.pi/180.

        # Copy stream stats
        str_stats = self.radialRF[0].stats

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
            for irow, trace in enumerate(self.radialRF):

                baz = trace.stats.baz
                OBS[irow] = trace.data[iz]
                H[irow, 0] = 1.0
                H[irow, 1] = np.cos(deg2rad*(baz-azim))
                H[irow, 2] = np.sin(deg2rad*(baz-azim))
                H[irow, 3] = np.cos(2.*deg2rad*(baz-azim))
                H[irow, 4] = np.sin(2.*deg2rad*(baz-azim))

            shift = 90.

            # Transverse component
            for irow, trace in enumerate(self.transvRF):

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
        receiver functions given the 5 pre-determined harmonics and
        a list of back-azimuth values. The receiver function signal
        parameters (length, sampling rate, etc.) will be identical
        to those in the stream of harmonic components.

        Parameters
        ----------
        baz_list : list
            List of back-azimuth directions over which to calculate
            the receiver functions. If no list is specified, the method
            will use the same back-azimuths as those in the original
            receiver function streams

        Attributes
        ----------
        radial_forward : :class:`~obspy.core.Stream`
            Stream containing the radial receiver functions
        transv_forward : :class:`~obspy.core.Stream`
            Stream containing the transverse receiver functions


        """

        if not hasattr(self, 'hstream'):
            raise(Exception("Decomposition has not been performed yet"))

        if not baz_list:
            print("Warning: no BAZ specified - using all baz from " +
                  "stored streams")
            baz_list = [tr.stats.baz for tr in self.radialRF]
        if not isinstance(baz_list, list):
            baz_list = [baz_list]

        # Some constants
        nz = len(self.hstream[0].data)
        deg2rad = np.pi/180.

        # Copy traces
        self.radial_forward = Stream()
        self.transv_forward = Stream()

        for baz in baz_list:
            trR = Trace(header=self.hstream[0].stats)
            trT = Trace(header=self.hstream[0].stats)

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
                trR.data[iz] = B[0]
                trT.data[iz] = -B[1]

            self.radial_forward.append(trR)
            self.transv_forward.append(trT)


    def plot(self, ymax=30., scale=10., save=False, title=None, form='png'):
        """
        Method to plot the 5 harmonic components.

        Parameters
        ----------
        ymax : float
            Maximum y axis value (time or depth) over which to
            plot the harmonic components
        scale : float
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

        # Initialize figure
        fig = plt.figure()
        plt.clf()

        # Get more control on subplots
        # ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.5])
        ax1 = fig.add_subplot(111)

        for i, trace in enumerate(self.hstream):
            # i += 1
            ax1.fill_betweenx(
                y, i, i+trace.data*scale,
                where=trace.data+1e-6 <= 0.,
                facecolor='blue',
                linewidth=0)
            ax1.fill_betweenx(
                y, i, i+trace.data*scale,
                where=trace.data+1e-6 >= 0.,
                facecolor='red',
                linewidth=0)

        ax1.set_ylim(ymax, 0)
        ax1.set_xlabel('Harmonic components')
        if title:
            ax1.set_title(title)

        ax1.set_xticks([0, 1, 2, 3, 4], minor=False)
        ax1.set_xticklabels(['$A$', '$B_1$', '$B_2$', '$C_1$', '$C_2$'], minor=False)
        off = ax1.xaxis.get_offset_text()
        ax1.tick_params(axis=u'x', pad=10)
        ax1.grid()

        if save:
            plt.savefig('FIGURES/'+sta+'.'+title+'.'+form, dpi=300,
                        bbox_inches='tight', format=form)
        plt.show()

    def save(self, file):
        """
        Saves harmonics object to file

        Parameters
        ----------
        file : str
            File name for Harmonics object

        """

        import pickle
        output = open(file, 'wb')
        pickle.dump(self, output)
        output.close()
