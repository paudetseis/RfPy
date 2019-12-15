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


def depth_azim(strV, strH, mindep=0., mdep=40.):
    """
    Function to decompose receiver function streams into back-azimuth harmonics
    and determine the main orientation alpha, obtained fro the surface to mdep.

    """
    print()
    print('Decomposing receiver functions into baz harmonics')

    # Some integers
    nbin = len(strV)
    nz = len(strV[0].data)
    naz = 180
    daz = np.float(360/naz)
    deg2rad = np.pi/180.

    # Define depth range over which to calculate azimuth
    indmin = int(mindep/strV[0].stats.delta)
    indmax = int(mdep/strV[0].stats.delta)

    # Copy stream stats
    str_stats = strV[0].stats

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
            irow = 0

            # Radial component
            for ibin in range(nbin):

                baz = strV[ibin].stats.sac.baz
                OBS[irow] = strV[ibin].data[iz]
                H[irow, 0] = 1.0
                H[irow, 1] = np.cos(deg2rad*(baz-azim))
                H[irow, 2] = np.sin(deg2rad*(baz-azim))
                H[irow, 3] = np.cos(2.*deg2rad*(baz-azim))
                H[irow, 4] = np.sin(2.*deg2rad*(baz-azim))
                irow += 1

            shift = 90.

            # Transverse component
            for ibin in range(nbin):

                baz = strV[ibin].stats.sac.baz
                OBS[irow] = strH[ibin].data[iz]
                H[irow, 0] = 0.0
                H[irow, 1] = np.cos(deg2rad*(baz+shift-azim))
                H[irow, 2] = np.sin(deg2rad*(baz+shift-azim))
                H[irow, 3] = np.cos(2.*deg2rad*(baz+shift/2.0-azim))
                H[irow, 4] = np.sin(2.*deg2rad*(baz+shift/2.0-azim))
                irow += 1

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
    stream = Stream(traces=[A, B1, B2, C1, C2])

    return stream, indaz*daz, [C0var, C1var, C2var, C3var, C4var]


def depth_fix_azim(strV, strH, azim=0.):
    """
    Function to decompose receiver function streams into back-azimuth harmonics
    at orientation given by azim

    """

    print('Decomposing receiver functions into baz harmonics for az = ', azim)

    # Some integers
    nbin = len(strV)
    nz = len(strV[0].data)
    deg2rad = np.pi/180.

    # Copy stream stats
    str_stats = strV[0].stats

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

        irow = 0

        # Radial component
        for ibin in range(nbin):

            baz = strV[ibin].stats.sac.baz
            OBS[irow] = strV[ibin].data[iz]
            H[irow, 0] = 1.0
            H[irow, 1] = np.cos(deg2rad*(baz-azim))
            H[irow, 2] = np.sin(deg2rad*(baz-azim))
            H[irow, 3] = np.cos(2.*deg2rad*(baz-azim))
            H[irow, 4] = np.sin(2.*deg2rad*(baz-azim))
            irow += 1

        shift = 90.

        # Transverse component
        for ibin in range(nbin):

            baz = strV[ibin].stats.sac.baz
            OBS[irow] = strH[ibin].data[iz]
            H[irow, 0] = 0.0
            H[irow, 1] = np.cos(deg2rad*(baz+shift-azim))
            H[irow, 2] = np.sin(deg2rad*(baz+shift-azim))
            H[irow, 3] = np.cos(2.*deg2rad*(baz+shift/2.0-azim))
            H[irow, 4] = np.sin(2.*deg2rad*(baz+shift/2.0-azim))
            irow += 1

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
    stream = Stream(traces=[A, B1, B2, C1, C2])

    return stream


def forward(trV, trH, str_harm, azim):

    # Some integers
    nz = len(trV.data)
    deg2rad = np.pi/180.

    # Copy traces
    rfV = trV.copy()
    rfH = trH.copy()

    # Loop over each depth step
    for iz in range(nz):

        # Initialize working arrays
        X = np.zeros(5)
        H = np.zeros((2, 5))

        baz = trV.stats.sac.baz

        # Fill up X array
        X[0] = str_harm[0].data[iz]
        X[1] = str_harm[1].data[iz]
        X[2] = str_harm[2].data[iz]
        X[3] = str_harm[3].data[iz]
        X[4] = str_harm[4].data[iz]

        # Fill up H arrays (for V and H)
        H[0, 0] = 1.0
        H[0, 1] = np.cos(deg2rad*(baz-azim))
        H[0, 2] = np.sin(deg2rad*(baz-azim))
        H[0, 3] = np.cos(2.*deg2rad*(baz-azim))
        H[0, 4] = np.sin(2.*deg2rad*(baz-azim))

        shift = 90.

        H[1, 0] = 0.0
        H[1, 1] = np.cos(deg2rad*(baz+shift-azim))
        H[1, 2] = np.sin(deg2rad*(baz+shift-azim))
        H[1, 3] = np.cos(2.*deg2rad*(baz+shift/2.0-azim))
        H[1, 4] = np.sin(2.*deg2rad*(baz+shift/2.0-azim))

        # Calculate dot product B = H*X
        B = np.dot(H, X)

        # Extract receiver functions
        rfV.data[iz] = B[0]
        rfH.data[iz] = -B[1]

    return rfV, rfH
