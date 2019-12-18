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
Functions to plot single station P-receiver functions as wiggle plots.

Options to plot by receiver functions # vs time/depth, by back-azimuth vs 
time/depth or slowness vs time. 

"""

# Import modules and functions
import os
import fnmatch
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from obspy.core import Stream, Trace, AttribDict
from scipy.interpolate import griddata

# PLot wiggles


def wiggle(str1, str2=None, tmax=30, normalize=True, save=False, title=None):

    print()
    print('Plotting Wiggles by trace number')
    print()

    # Time axis
    t = np.arange(str1[0].stats.npts)/str1[0].stats.sampling_rate

    # Normalize?
    if normalize:
        ntr = len(str1)
        maxamp = np.max([np.max(np.abs(tr.data[t<tmax])) for tr in str1])

    f = plt.figure()

    # Clear figure
    plt.clf()

    if str2 is not None:
        ax1 = f.add_subplot(121)
    else:
        ax1 = f.add_subplot(111)

    # loop in stream
    for itr, tr in enumerate(str1):

        # Fill positive in red, negative in blue
        ax1.fill_between(t,itr+1,itr+1+tr.data/maxamp/2.,\
               where=tr.data+1e-10<=0.,facecolor='blue',linewidth=0)
        ax1.fill_between(t,itr+1,itr+1+tr.data/maxamp/2.,\
               where=tr.data+1e-10>=0.,facecolor='red',linewidth=0)
        # plt.plot(t, y+tr.data*maxamp, c='k')

    ax1.set_xlim(0, tmax)
    ax1.set_ylabel('Radial RF')
    ax1.grid()

    if str2 is not None:
        ax2 = f.add_subplot(122)


        # loop in stream
        for itr, tr in enumerate(str2):

            # Fill positive in red, negative in blue
            ax2.fill_between(t,itr+1,itr+1+tr.data/maxamp/2.,\
                   where=tr.data+1e-10<=0.,facecolor='blue',linewidth=0)
            ax2.fill_between(t,itr+1,itr+1+tr.data/maxamp/2.,\
                   where=tr.data+1e-10>=0.,facecolor='red',linewidth=0)
            # plt.plot(t, y+tr.data*maxamp, c='k')

        ax2.set_xlim(0, tmax)
        ax2.set_ylabel('Transverse RF')
        ax2.grid()

    plt.suptitle('Station '+str1[0].stats.station)

    if save:
        plt.savefig('RF_PLOTS/'+str1[0].stats.station+title+'.png', 
            dpi=300, bbox_inches='tight')
    else:
        plt.show()


# PLot wiggles along line
def wiggle_line(str1, str2, zmax=30, maxval=1000, ll=[], profile=[],
                save=False, title=None, suptitle=None):

    print()
    print('Plotting Wiggles by distance along line')

    # Time axis
    t = np.arange(str1[0].stats.npts)/str1[0].stats.sampling_rate

    # Initialize figure
    #fig = plt.figure(figsize=(15,6))
    #fig = plt.figure(figsize=(10,4))
    fig = plt.figure(0, figsize=(6.5, 2.75))
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.85, 0.85, 0.1])
    ax2 = fig.add_axes([0.1, 0.475, 0.85, 0.35])
    ax3 = fig.add_axes([0.1, 0.1, 0.85, 0.35])
    #ax1 = fig.add_axes([0.1, 0.85, 0.425, 0.1])
    #ax2 = fig.add_axes([0.1, 0.475, 0.425, 0.35])
    #ax3 = fig.add_axes([0.1, 0.1, 0.425, 0.35])

    # First elevation plot
    ax1.plot(profile[0], profile[1]/1000., 'k')
    #ax1.set_ylabel('Elevation (km)')
    ax1.set_ylim(-1, 1.)
    ax1.set_xticks(())
    ax1.set_yticks((-1, 0, 1))
    # ax1.set_xlim((-25,300))
    # ax1.set_xlim((-85,150))

    # Initialize index
    i = 0

    # loop in stream
    for tr in str1:

        # Get distance from list
        x = ll[i]
        #print('str1, x= ',x)

        # Update count
        i += 1

        # Fill positive in red, negative in blue
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    ax2.set_ylim(zmax, 0)
    ax2.set_ylabel('Time (s)', size=10)
    # ax2.grid()
    ax2.set_xticks(())
    # ax2.set_xlim((-25,300))
    # ax2.set_xlim((-85,150))

    # Re-initialize index
    i = 0

    # loop in stream
    for tr in str2:

        # Get distance from list
        x = ll[i]
        #print('str2, x= ',x)

        # Update count
        i += 1

        # Fill positive in red, negative in blue
        ax3.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax3.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    ax3.set_ylim(zmax, 0)
    ax3.set_ylabel('Time (s)', size=10)
    ax3.set_xlabel('Distanceh (km)', size=10)
    #ax3.set_xlabel('Distance from trench (km)', size=10)
    # ax3.set_xlim((-25,300))
    # ax3.set_xlim((-85,150))
    # ax3.grid()

    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax3.tick_params(axis='both', which='major', labelsize=10)

    if suptitle:
        plt.suptitle(suptitle)

    if save:
        plt.savefig('PLOTS/'+title+'.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()


# PLot wiggles along line
def wiggle_cascadia(str1, zmax=30, maxval=1000, ll=[], profile=[],
                    save=False, title=None, suptitle=None):

    print()
    print('Plotting Wiggles by distance along line')

    # Time axis
    t = np.arange(str1[0].stats.npts)/str1[0].stats.sampling_rate

    # Initialize figure
    #fig = plt.figure(figsize=(15,6))
    #fig = plt.figure(figsize=(10,4))
    fig = plt.figure(0, figsize=(6.5, 2.75))
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.75, 0.85, 0.1])
    ax2 = fig.add_axes([0.1, 0.15, 0.85, 0.6])

    # First elevation plot
    #ax1.fill_between(profile[0], 0., profile[1]/1000., facecolor='lightblue')
    #ax1.fill_between(profile[0], -2., profile[1]/1000., facecolor='lightgrey')
    ax1.plot(profile[0], profile[1]/1000., 'k')
    #ax1.set_ylabel('Elevation (km)')
    # ax1.set_ylim(-2,0.8)
    ax1.set_xticks(())
    #ax1.set_yticks((-1, 0, 1))

    # Initialize index
    i = 0

    # loop in stream
    for tr in str1:

        # Get distance from list
        x = ll[i]

        # Update count
        i += 1

        # Fill positive in red, negative in blue
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    ax2.set_ylim(zmax, 0)
    ax2.set_ylabel('Time (s)', size=10)
    ax2.set_xlabel('Distance along line (km)', size=10)
    # ax2.grid()
    # ax2.set_xticks(())
    # ax2.set_xlim((-25,300))
    ax2.set_xlim((680, 1150))

    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax2.tick_params(axis='both', which='major', labelsize=10)

    if suptitle:
        plt.suptitle(suptitle)

    if save:
        plt.savefig('PLOTS/'+title+'.eps', dpi=300,
                    bbox_inches='tight', format='eps')
    else:
        plt.show()

# PLot wiggles along line


def wiggle_line2(str1, zmax=30, maxval=1000, ll=[], profile=[],
                 save=False, title=None, suptitle=None):

    print()
    print('Plotting Wiggles by distance along line')

    # Time axis
    t = np.arange(str1[0].stats.npts)/str1[0].stats.sampling_rate

    # Initialize figure
    #fig = plt.figure(figsize=(15,6))
    #fig = plt.figure(figsize=(10,4))
    fig = plt.figure(0, figsize=(8.5, 3.75))
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.8, 0.85, 0.2])
    axs = fig.add_axes([0.1, 0.6, 0.85, 0.2])
    ax2 = fig.add_axes([0.1, 0.1, 0.85, 0.5])

    # First elevation plot
    ax1.plot(profile[0], profile[1]/1000., 'k')
    ax1.set_xticks(())
    ax1.set_yticks((1, 2))
    ax1.set_ylabel('Elev. (km)', size=10)
    axs.set_ylabel('Sed. (km)', size=10)
    axs.set_xticks(())
    axs.set_ylim((-5., 0.))
    axs.set_yticks((-4, -2, 0))

    sed = np.loadtxt('line2_sed.xyz').transpose()[2]

    # Initialize index
    i = 0

    # loop in stream
    for tr in str1:

        # Get distance from list
        x = ll[i]
        ss = -sed[i]

        # Update count
        i += 1

        axs.scatter(x, ss)

        # Fill positive in red, negative in blue
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    ax2.set_ylim(zmax, 0)
    ax2.set_ylabel('Time (s)', size=10)
    ax2.set_xlabel('Distance along line (km)', size=10)
    ax2.set_xlim((680, 1150))
    axs.set_xlim((680, 1150))

    ax1.tick_params(axis='both', which='major', labelsize=10)
    axs.tick_params(axis='both', which='major', labelsize=10)
    ax2.tick_params(axis='both', which='major', labelsize=10)

    if suptitle:
        plt.suptitle(suptitle)

    if save:
        plt.savefig('PLOTS/'+title+'.eps', dpi=300,
                    bbox_inches='tight', format='eps')
    else:
        plt.show()


# PLot wiggles according to either baz or slowness
def wiggle_bins(str1, str2, tr1, tr2, sta, btyp='baz', xmax=30,
                xtyp='time', scale=None, save=False, title=None, wtype='P'):

    if not (btyp == 'baz' or btyp == 'slow' or btyp == 'dist'):
        print('type has to be "baz" or "slow" or "dist"')
        return
    if not (xtyp == 'time' or xtyp == 'depth'):
        print('type has to be "time" or "depth"')
        return
    if btyp == 'slow' and xtyp == 'depth':
        print('Cannot plot by slowness if data is migrated')
        return

    print()
    print('Plotting Wiggles by '+btyp)

    # X axis
    nn = str1[0].stats.npts
    sr = str1[0].stats.sampling_rate
    if wtype == 'P':
        x = np.arange(nn)/sr
    elif wtype == 'S' or wtype == 'SKS':
        x = np.arange(-nn/2, nn/2)/sr
        tr1.data = np.fft.fftshift(tr1.data)
        tr2.data = np.fft.fftshift(tr2.data)

    # Initialize figure
    fig = plt.figure()
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.825, 0.3, 0.05])
    ax2 = fig.add_axes([0.1, 0.1, 0.3, 0.7])
    ax3 = fig.add_axes([0.45, 0.825, 0.3, 0.05])
    ax4 = fig.add_axes([0.45, 0.1, 0.3, 0.7])

    # Plot stack of all SV traces on top left
    ax1.fill_between(
        x, 0., tr1.data,
        where=tr1.data+1e-6 <= 0.,
        facecolor='blue',
        linewidth=0)
    ax1.fill_between(
        x, 0., tr1.data,
        where=tr1.data+1e-6 >= 0.,
        facecolor='red',
        linewidth=0)
    ax1.set_ylim(-np.max(np.abs(tr1.data)), np.max(np.abs(tr1.data)))
    ax1.set_yticks(())
    ax1.set_xticks(())
    ax1.set_title('Radial')
    if wtype == 'P':
        ax1.set_xlim(0, xmax)
    elif wtype == 'S' or wtype == 'SKS':
        ax1.set_xlim(-20., xmax)

    # Plot binned SV traces in back-azimuth on bottom left
    for tr in str1:

        if wtype == 'S' or wtype == 'SKS':
            tr.data = np.fft.fftshift(tr.data)

        if scale:
            maxval = scale
            # Define y axis
            if btyp == 'baz':
                y = tr.stats.baz
            elif btyp == 'slow':
                y = tr.stats.slow
            elif btyp == 'dist':
                y = tr.stats.sac.user0
        else:
            # Define y axis
            if btyp == 'baz':
                y = tr.stats.baz
                maxval = 100
            elif btyp == 'slow':
                y = tr.stats.user0
                maxval = 0.02
            elif btyp == 'dist':
                y = tr.stats.sac.user0
                maxval = 20

        # Fill positive in red, negative in blue
        ax2.fill_between(
            x, y, y+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax2.fill_between(
            x, y, y+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    if wtype == 'P':
        ax2.set_xlim(0, xmax)
    elif wtype == 'S' or wtype == 'SKS':
        ax2.set_xlim(-20., xmax)

    if btyp == 'baz':
        ax2.set_ylim(-5, 370)
        ax2.set_ylabel('Back-azimuth (deg)')

    elif btyp == 'slow':
        if wtype == 'P':
            ax2.set_ylim(0.038, 0.082)
        elif wtype == 'S':
            ax2.set_ylim(0.07, 0.125)
        elif wtype == 'SKS':
            ax2.set_ylim(0.03, 0.06)
        ax2.set_ylabel('Slowness (s/km)')
    elif btyp == 'dist':
        if wtype == 'P':
            ax2.set_ylim(28., 92.)
        elif wtype == 'S':
            ax2.set_ylim(53., 107.)
        elif wtype == 'SKS':
            ax2.set_ylim(83., 117.)
        ax2.set_ylabel('Distance (deg)')

    if xtyp == 'time':
        ax2.set_xlabel('Time (sec)')
    elif xtyp == 'depth':
        ax2.set_xlabel('Depth (km)')
    ax2.grid()

    # Plot stack of all SH traces on top right
    ax3.fill_between(
        x, 0., tr2.data,
        where=tr2.data+1e-6 <= 0.,
        facecolor='blue',
        linewidth=0)
    ax3.fill_between(
        x, 0., tr2.data,
        where=tr2.data+1e-6 >= 0.,
        facecolor='red',
        linewidth=0)
    if wtype == 'P':
        ax3.set_xlim(0, xmax)
    elif wtype == 'S' or wtype == 'SKS':
        ax3.set_xlim(-20., xmax)
    ax3.set_ylim(-np.max(np.abs(tr1.data)), np.max(np.abs(tr1.data)))
    ax3.set_yticks(())
    ax3.set_xticks(())
    ax3.set_title('Transverse')
    # ax3.set_title('Synthetic')
    # ax3.set_title('SH')

    # Plot binned SH traces in back-azimuth on bottom right
    for tr in str2:

        if wtype == 'S' or wtype == 'SKS':
            tr.data = np.fft.fftshift(tr.data)

        if scale:
            maxval = scale
            # Define y axis
            if btyp == 'baz':
                y = tr.stats.baz
            elif btyp == 'slow':
                y = tr.stats.slow
            elif btyp == 'dist':
                y = tr.stats.sac.user0
        else:
            # Define y axis
            if btyp == 'baz':
                y = tr.stats.baz
                maxval = 150
            elif btyp == 'slow':
                y = tr.stats.slow
                maxval = 0.02
            elif btyp == 'dist':
                y = tr.stats.sac.user0
                maxval = 20

        # Fill positive in red, negative in blue
        ax4.fill_between(
            x, y, y+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax4.fill_between(
            x, y, y+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    if wtype == 'P':
        ax4.set_xlim(0, xmax)
    elif wtype == 'S' or wtype == 'SKS':
        ax4.set_xlim(-20., xmax)

    if btyp == 'baz':
        ax4.set_ylim(-5, 370)

    elif btyp == 'slow':
        if wtype == 'P':
            ax4.set_ylim(0.038, 0.082)
        elif wtype == 'S':
            ax4.set_ylim(0.074, 0.125)
        elif wtype == 'SKS':
            ax4.set_ylim(0.03, 0.06)
    elif btyp == 'dist':
        if wtype == 'P':
            ax4.set_ylim(28., 92.)
        elif wtype == 'S':
            ax4.set_ylim(53., 107.)
        elif wtype == 'SKS':
            ax4.set_ylim(83., 117.)

    if xtyp == 'time':
        ax4.set_xlabel('Time (sec)')
    elif xtyp == 'depth':
        ax4.set_xlabel('Depth (km)')
    ax4.set_yticklabels([])
    ax4.grid()

    #plt.suptitle('Station '+sta)

    if save:
        # plt.savefig('RF_PLOTS/'+sta+title+'.png',dpi=300,bbox_inches='tight')
        plt.savefig('RF_PLOTS/'+sta+title+'.eps', format='eps')
    else:
        plt.show()

# Plot wiggles for harmonic decomposition


def wiggle_harmonics(stream, xmax=30, maxval=10, save=False,
                     title=None, nz=50, dep=None, vp=None):

    # Y axis
    y = np.arange(stream[0].stats.npts)/stream[0].stats.sampling_rate

    # Station name
    sta = stream[0].stats.station
    print(sta)

    # Initialize count
    i = 0

    # Initialize figure
    fig = plt.figure()
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.5])
    #ax2 = fig.add_axes([0.85, 0.1, 0.1, 0.5])

    #ax1 = fig.add_axes([0.1, 0.825, 0.3, 0.05])
    #ax2 = fig.add_axes([0.1, 0.1, 0.3, 0.7])
    #ax3 = fig.add_axes([0.45, 0.825, 0.3, 0.05])
    #ax4 = fig.add_axes([0.45, 0.1, 0.3, 0.7])

    for trace in stream:
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

    ax1.set_ylim(xmax, 0)
    ax1.set_ylabel('Depth (km)')
    #ax1.set_xlabel('Harmonic component')
    # ax1.set_title('Station '+sta)
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels[1] = '$A$'
    labels[1] = '$A$'
    labels[2] = '$B_1$'
    labels[3] = '$B_2$'
    labels[4] = '$C_1$'
    labels[5] = '$C_2$'
    # labels[2] = '$B_{\parallel}$'
    # labels[3] = '$B_{\perp}$'
    # labels[4] = '$C_{\parallel}$'
    # labels[5] = '$C_{\perp}$'
    ax1.set_xticklabels(labels)
    off = ax1.xaxis.get_offset_text()
    #tcex = ax1.xaxis.get_ticklabel_extents()
    ax1.tick_params(axis=u'x', pad=10)
    # ax1.xaxis.set_label_position('top')
    #ax1.set_xticks(['A', 'Bper', 'Bpar', 'Cpar', 'Cpar'])
    ax1.grid()

    if (dep is None) and (vp is None):
        dep = np.array([0., 4., 8., 14., 25.9, 35.7, 45., 110.])
        vp = np.array([4.0, 5.9, 6.2, 6.3, 6.8, 7.2, 8.0, 8.1])

    # Get interpolated velocity profile
    #idep = np.linspace(dep.min(), dep.max(), nz)
    #ivp = sp.interpolate.interp1d(dep, vp, kind='slinear')(idep)
    #ax2.plot(ivp, idep, 'purple',linewidth=1.5)
    # ax2.set_ylim(xmax,0)
    # ax2.set_xlim(4,8.5)
    # ax2.yaxis.tick_right()
    #xt = [4, 5, 6, 7, 8]
    # ax2.set_xticks(xt)
    #ax2.set_xlabel('Vp (km/s)')

    if save:
        plt.savefig('RF_PLOTS/'+sta+title+'.eps', dpi=300,
                    bbox_inches='tight', format='eps')
    else:
        plt.show()

# Plot wiggles for harmonic decomposition


def wiggle_harmonics_comp(
        stream1, stream2, xmax=30, maxval=10, save=False,
        title=None, nz=50, dep=None, vp=None):

    print('wiggle harmonics comp')

    # Y axis
    y = np.arange(stream1[0].stats.npts)/stream1[0].stats.sampling_rate

    # Station name
    sta = stream1[0].stats.station

    # Initialize count
    i = 0

    # Initialize figure
    fig = plt.figure()
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.1, 0.7, 0.5])

    for trace in stream1:
        i += 1
        ax1.plot(i+trace.data*maxval, y, c='r')

    i = 0

    for trace in stream2:
        i += 1
        ax1.plot(i+trace.data*maxval, y, c='k')

    ax1.set_ylim(xmax, 0)
    ax1.set_ylabel('Depth (km)')
    #ax1.set_xlabel('Harmonic component')
    ax1.set_title('Station '+sta)
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels[1] = '$A$'
    labels[2] = '$B_{\parallel}$'
    labels[3] = '$B_{\perp}$'
    labels[4] = '$C_{\parallel}$'
    labels[5] = '$C_{\perp}$'
    ax1.set_xticklabels(labels)
    off = ax1.xaxis.get_offset_text()
    #tcex = ax1.xaxis.get_ticklabel_extents()
    ax1.tick_params(axis=u'x', pad=10)
    # ax1.xaxis.set_label_position('top')
    #ax1.set_xticks(['A', 'Bper', 'Bpar', 'Cpar', 'Cpar'])
    ax1.grid()

    if (dep is None) and (vp is None):
        dep = np.array([0., 4., 8., 14., 25.9, 35.7, 45., 110.])
        vp = np.array([4.0, 5.9, 6.2, 6.3, 6.8, 7.2, 8.0, 8.1])

    if save:
        plt.savefig('RF_PLOTS/'+sta+title+'.eps', dpi=300,
                    bbox_inches='tight', format='eps')
    else:
        plt.show()


# Plot wiggles for harmonic decomposition
def wiggle_tohoku(stream1, stream2, xmax=30, maxval=10, save=False,
                  title=None):

    # Y axis
    y = np.arange(stream1[0].stats.npts)/stream1[0].stats.sampling_rate

    # Station name
    sta = stream1[0].stats.station

    # Initialize count
    i = 0

    # Clear figure
    plt.clf()

    for trace in stream1:
        i += 1
        plt.plot(i+trace.data*maxval, y, 'b-', label='Pre Tohoku')

    i = 0
    for trace in stream2:
        i += 1
        plt.plot(i+trace.data*maxval, y, 'r--', label='Post Tohoku')

    plt.ylim(xmax, 0)
    plt.ylabel('Depth (km)')
    plt.xlabel('Harmonic component')
    plt.title('Station '+sta)
    #plt.xticks(('A', 'Bper', 'Bpar', 'Cpar', 'Cpar'))
    plt.grid()
    plt.legend()

    if save:
        plt.savefig('RF_PLOTS/'+sta+title+'.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()

# PLot dt wiggles


def wiggle_dt(str1, str2, xmax=30, maxval=1, save=False, title=None):

    print()
    print('Plotting Wiggles by start time')

    # Station name
    sta = str1[0].stats.station

    # Time axis
    t = np.arange(str1[0].stats.npts)/str1[0].stats.sampling_rate

    # Clear figure
    plt.clf()

    # First SV subplot
    plt.subplot(121)

    # loop in stream
    for tr in str1:

        # Timestamp
        y = (tr.stats.starttime - str1[0].stats.starttime)/3600./24./30.

        # Fill positive in red, negative in blue
        plt.fill_between(
            t, y, y+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        plt.fill_between(
            t, y, y+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    plt.xlim(0, xmax)
    #plt.ylabel('Radial RF')
    plt.grid()

    # Second SH subplot
    plt.subplot(122)

    # loop in stream
    for tr in str2:

        # Timestamp
        y = (tr.stats.starttime - str1[0].stats.starttime)/3600./24./30.

        # Fill positive in red, negative in blue
        plt.fill_between(
            t, y, y+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        plt.fill_between(
            t, y, y+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    plt.xlim(0, xmax)
    #plt.ylabel('Transverse RF')
    plt.suptitle('Station '+sta)
    plt.grid()

    if save:
        plt.savefig('RF_PLOTS/'+sta+title+'.png', dpi=300, bbox_inches='tight')
    else:
        plt.show()


# PLot interpolated image along line
def interp_line(str1, str2, zmax=30, maxval=1000, ll=[], profile=[],
                save=False, title=None):

    print()
    print('Plotting interpolated image by distance along line')

    # Time axis
    t = np.arange(str1[0].stats.npts)/str1[0].stats.sampling_rate

    # Initialize figure
    #fig = plt.figure(figsize=(15,6))
    #fig = plt.figure(figsize=(10,4))
    fig = plt.figure(0, figsize=(6.5, 2.75))
    plt.clf()

    # Get more control on subplots
    ax1 = fig.add_axes([0.1, 0.85, 0.85, 0.1])
    ax2 = fig.add_axes([0.1, 0.475, 0.85, 0.35])
    ax3 = fig.add_axes([0.1, 0.1, 0.85, 0.35])
    #ax1 = fig.add_axes([0.1, 0.85, 0.425, 0.1])
    #ax2 = fig.add_axes([0.1, 0.475, 0.425, 0.35])
    #ax3 = fig.add_axes([0.1, 0.1, 0.425, 0.35])

    # First elevation plot
    ax1.plot(profile[0], profile[1]/1000., 'k')
    #ax1.set_ylabel('Elevation (km)')
    ax1.set_ylim(-1, 1.)
    ax1.set_xticks(())
    ax1.set_yticks((-1, 0, 1))
    # ax1.set_xlim((-25,300))
    ax1.set_xlim((-100, 175))

    # Initialize index
    i = 0

#    x = ll
#    gridxi, gridt = np.arange(-100,175,5
#    grixt = t
#
#    lst = []
#    gridR = [lst.append(st.data) for st in str1]
#    grid_int = griddata(

    # loop in stream
    for tr in str1:

        # Get distance from list
        x = ll[i]
        print('str1, x= ', x)

        # Update count
        i += 1

        # Fill positive in red, negative in blue
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax2.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    ax2.set_ylim(zmax, 0)
    ax2.set_ylabel('Time (s)', size=10)
    ax2.grid()
    ax2.set_xticks(())
    # ax2.set_xlim((-25,300))
    ax2.set_xlim((-100, 175))
    ax2.grid()

    # Re-initialize index
    i = 0

    # loop in stream
    for tr in str2:

        # Get distance from list
        x = ll[i]
        print('str2, x= ', x)

        # Update count
        i += 1

        # Fill positive in red, negative in blue
        ax3.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            facecolor='blue',
            linewidth=0)
        ax3.fill_betweenx(
            t, x, x+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='red',
            linewidth=0)

    ax3.set_ylim(zmax, 0)
    ax3.set_ylabel('Time (s)', size=10)
    ax3.set_xlabel('Distance (km)', size=10)
    # ax3.set_xlim((-25,300))
    ax3.set_xlim((-100, 175))
    # ax3.grid()
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax3.tick_params(axis='both', which='major', labelsize=10)

    if save:
        plt.savefig('PLOTS/'+title+'.eps', dpi=300,
                    bbox_inches='tight', format='eps')
        plt.savefig('PLOTS/'+title+'.png')
    else:
        plt.show()
