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
import numpy as np
import scipy as sp
from obspy.core import Stream, Trace, AttribDict
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


def wiggle(stream1, stream2=None, sort=None, tmin=0., tmax=30, normalize=True,
           save=False, title=None, form='png'):
    """
    Function to plot receiver function traces by index in stream. By default, 
    only one stream is required, which produces a Figure with a single panel.
    Optionally, if a second stream is present, the Figure will contain two panels.
    The stream(s) can be sorted by stats attribute ``sort``, normalized, and
    the Figure can be saved as .eps.

    Parameters
    ----------
    stream1 : :class:`~obspy.core.Stream`
        Stream of receiver function traces
    stream2 : :class:`~obspy.core.Stream`
        Stream of receiver function traces, for a two-panel Figure
    sort : str
        Name of attribute in the ``stats`` attribute of indivudial traces, used
        to sort the traces in the Stream(s).
    xmax : float
        Maximum x-axis value displayed in the Figure.
    normalize : bool
        Whether or not to normalize the traces according to the max amplitude 
        in ``stream1``
    save : bool
        Whether or not to save the Figure 
    title : str
        Title of plot

    """

    if sort:
        try:
            stream1.traces.sort(key=lambda x: x.stats[sort], reverse=False)
        except:
            print("Warning: stats attribute " + sort +
                  " is not available in stats")
            pass

    # Time axis
    nn = stream1[0].stats.npts
    sr = stream1[0].stats.sampling_rate
    time = np.arange(-nn/2+1, nn/2+1)/sr

    # Normalize?
    if normalize:
        ntr = len(stream1)
        maxamp = np.median(
            [np.max(np.abs(
                tr.data[(time > tmin) & (time < tmax)])) for tr in stream1])

    f = plt.figure()

    # Clear figure
    plt.clf()

    if stream2 is not None:
        ax1 = f.add_subplot(121)
    else:
        ax1 = f.add_subplot(111)

    # loop in stream
    for itr, tr in enumerate(stream1):

        # Fill positive in red, negative in blue
        ax1.fill_between(
            time, itr+1, itr+1+tr.data/maxamp/2.,
            where=tr.data+1e-10 <= 0., facecolor='blue', linewidth=0)
        ax1.fill_between(
            time, itr+1, itr+1+tr.data/maxamp/2.,
            where=tr.data+1e-10 >= 0., facecolor='red', linewidth=0)
        # plt.plot(t, y+tr.data*maxamp, c='k')

    ax1.set_xlim(tmin, tmax)
    ax1.set_ylabel('Radial RF')
    ax1.grid()

    if stream2 is not None:
        ax2 = f.add_subplot(122)

        # loop in stream
        for itr, tr in enumerate(stream2):

            # Fill positive in red, negative in blue
            ax2.fill_between(
                time, itr+1, itr+1+tr.data/maxamp/2.,
                where=tr.data+1e-10 <= 0., facecolor='blue', linewidth=0)
            ax2.fill_between(
                time, itr+1, itr+1+tr.data/maxamp/2.,
                where=tr.data+1e-10 >= 0., facecolor='red', linewidth=0)

        ax2.set_xlim(tmin, tmax)
        ax2.set_ylabel('Transverse RF')
        ax2.grid()

    plt.suptitle('Station ' + stream1[0].stats.station)

    if save:
        plt.savefig('RF_PLOTS/' + stream1[0].stats.station +
                    '.' + title + '.' + form,
                    dpi=300, bbox_inches='tight', format=form)
    else:
        plt.show()

    plt.close()


# PLot wiggles according to either baz or slowness
def wiggle_bins(stream1, stream2=None, tr1=None, tr2=None,
                btyp='baz', tmin=0., tmax=30., xtyp='time', scale=None,
                norm=None, save=False, title=None, form='png'):
    """
    Function to plot receiver function according to either baz or
    slowness bins. By default, 
    only one stream is required, which produces a Figure with a single panel.
    Optionally, if a second stream is present, the Figure will contain two panels.
    If the single trace arguments are present, they will be plotted at the top.

    Parameters
    ----------
    stream1 : :class:`~obspy.core.Stream`
        Stream of receiver function traces
    stream2 : :class:`~obspy.core.Stream`
        Stream of receiver function traces, for a two-panel Figure
    tr1 : :class:`~obspy.core.Trace`
        Trace to plot at the top of ``stream1``
    tr2 : :class:`~obspy.core.Trace`
        Trace to plot at the top of ``stream2``
    btyp : str
        Type of plot to produce (either 'baz', 'slow', or 'dist')
    tmax : float
        Maximum x-axis value displayed in the Figure.
    xtyp : str
        Type of x-axis label (either 'time' or 'depth')
    scale : float
        Scale factor applied to trace amplitudes for plotting
    save : bool
        Whether or not to save the Figure 
    title : str
        Title of plot

    """

    if not (btyp == 'baz' or btyp == 'slow' or btyp == 'dist'):
        raise(Exception("Type has to be 'baz' or 'slow' or 'dist'"))
    if not (xtyp == 'time' or xtyp == 'depth'):
        raise(Exception("Type has to be 'time' or 'depth'"))
    if btyp == 'slow' and xtyp == 'depth':
        raise(Exception("Cannot plot by slowness if data is migrated"))

    # Figure out scaling here
    if scale:
        maxval = scale
        maxvalT = maxval
    else:
        if norm:
            for tr in stream1:
                tr.data /= norm
            if stream2:
                for tr in stream2:
                    tr.data /= norm
            if btyp == 'baz':
                maxval = 10.
                maxvalT = maxval
            elif btyp == 'slow':
                maxval = 0.001
                maxvalT = 2.*maxval
            elif btyp == 'dist':
                maxval = 1
                maxvalT = maxval
        else:
            if btyp == 'baz':
                maxval = 100
                maxvalT = maxval
            elif btyp == 'slow':
                maxval = 0.02
                maxvalT = maxval
            elif btyp == 'dist':
                maxval = 20
                maxvalT = maxval

    # Time axis
    nn = stream1[0].stats.npts
    sr = stream1[0].stats.sampling_rate
    time = np.arange(-nn/2, nn/2)/sr

    # Initialize figure
    fig = plt.figure()
    plt.clf()

    if stream2 and tr1 and tr2:
        # Get more control on subplots
        ax1 = fig.add_axes([0.1, 0.825, 0.3, 0.05])
        ax2 = fig.add_axes([0.1, 0.1, 0.3, 0.7])
        ax3 = fig.add_axes([0.45, 0.825, 0.3, 0.05])
        ax4 = fig.add_axes([0.45, 0.1, 0.3, 0.7])
    elif stream2 and not tr1 and not tr2:
        ax1 = None
        ax2 = fig.add_subplot(121)
        ax3 = None
        ax4 = fig.add_subplot(122)
    else:
        ax1 = None
        ax2 = fig.add_subplot(111)
        ax3 = None
        ax4 = None

    if ax1:

        # Plot stack of all SV traces on top left
        ax1.fill_between(
            time, 0., tr1.data,
            where=tr1.data+1e-6 <= 0.,
            # facecolor='grey',
            facecolor='red',
            linewidth=0)
        ax1.fill_between(
            time, 0., tr1.data,
            where=tr1.data+1e-6 >= 0.,
            # facecolor='red',
            facecolor='k',
            linewidth=0)
        ax1.plot(time, tr1.data,
            linewidth=0.25, c='k')
        ylim = np.max([np.max(np.abs(tr1.data)), np.max(np.abs(tr2.data))])
        if btyp=='slow':
            ylimT = ylim/2.
        else:
            ylimT = ylim
        ax1.set_ylim(-1.*ylim, ylim)
        ax1.set_yticks(())
        ax1.set_xticks(())
        ax1.set_title('Radial')
        ax1.set_xlim(tmin, tmax)

    # Plot binned SV traces in back-azimuth on bottom left
    for tr in stream1:

        # Define y axis
        if btyp == 'baz':
            y = tr.stats.baz
        elif btyp == 'slow':
            y = tr.stats.slow
        elif btyp == 'dist':
            y = tr.stats.dist

        # Fill positive in red, negative in blue
        ax2.fill_between(
            time, y, y+tr.data*maxval,
            where=tr.data+1e-6 <= 0.,
            # facecolor='grey',
            facecolor='red',
            linewidth=0)
        ax2.fill_between(
            time, y, y+tr.data*maxval,
            where=tr.data+1e-6 >= 0.,
            facecolor='k',
            # facecolor='red',
            linewidth=0)
        ax2.plot(time, y+tr.data*maxval,
            linewidth=0.25, c='k')

    ax2.set_xlim(tmin, tmax)

    if btyp == 'baz':
        ax2.set_ylim(-5, 370)
        ax2.set_ylabel('Back-azimuth (deg)')
    elif btyp == 'slow':
        # ax2.set_ylim(0.03, 0.13)
        ax2.set_ylim(0.038, 0.082)
        ax2.set_ylabel('Slowness (s/km)')
    elif btyp == 'dist':
        ax2.set_ylim(28., 92.)
        ax2.set_ylabel('Distance (deg)')

    if xtyp == 'time':
        ax2.set_xlabel('Time (sec)')
    elif xtyp == 'depth':
        ax2.set_xlabel('Depth (km)')
    ax2.grid(ls=':')
    if not ax1:
        ax2.set_title('Radial')

    if ax3:

        # Plot stack of all SH traces on top right
        ax3.fill_between(
            time, 0., tr2.data,
            where=tr2.data+1e-6 <= 0.,
            # facecolor='grey',
            facecolor='red',
            linewidth=0)
        ax3.fill_between(
            time, 0., tr2.data,
            where=tr2.data+1e-6 >= 0.,
            facecolor='k',
            # facecolor='red',
            linewidth=0)
        ax3.plot(time, tr2.data,
            linewidth=0.25, c='k')
        ax3.set_xlim(tmin, tmax)
        ax3.set_ylim(-1.*ylimT, ylimT)
        ax3.set_yticks(())
        ax3.set_xticks(())
        ax3.set_title('Transverse')

    if ax4:
        # Plot binned SH traces in back-azimuth on bottom right
        for tr in stream2:

            # Define y axis
            if btyp == 'baz':
                y = tr.stats.baz
            elif btyp == 'slow':
                y = tr.stats.slow
            elif btyp == 'dist':
                y = tr.stats.dist

            # Fill positive in red, negative in blue
            ax4.fill_between(
                time, y, y+tr.data*maxvalT,
                where=tr.data+1e-6 <= 0.,
                facecolor='red',
                # facecolor='grey',
                linewidth=0)
            ax4.fill_between(
                time, y, y+tr.data*maxvalT,
                where=tr.data+1e-6 >= 0.,
                # facecolor='red',
                facecolor='k',
                linewidth=0)
            ax4.plot(time, y+tr.data*maxvalT,
                linewidth=0.25, c='k')

        ax4.set_xlim(tmin, tmax)

        if btyp == 'baz':
            ax4.set_ylim(-5, 370)
        elif btyp == 'slow':
            # ax4.set_ylim(0.08, 0.13)
            ax4.set_ylim(0.038, 0.082)
        elif btyp == 'dist':
            ax4.set_ylim(28., 92.)

        if xtyp == 'time':
            ax4.set_xlabel('Time (sec)')
        elif xtyp == 'depth':
            ax4.set_xlabel('Depth (km)')
        ax4.set_yticklabels([])
        ax4.grid(ls=':')
        if not ax3:
            ax4.set_title('Transverse')

    if title:
        plt.suptitle(title)

    if save:
        plt.savefig('RF_PLOTS/' + stream1[0].stats.station +
                    '.' + title + '.' + form, format=form)
    else:
        plt.show()

    plt.close()


def wiggle_single_event(rfdata, filt=None, pre_filt=None, trange=None):

    lqtcopy = rfdata.data.copy()
    rfcopy = rfdata.rf.copy()
    nn = lqtdata[0].stats.npts
    sr = lqtdata[0].stats.sampling_rate
    taxis = np.arange(-nn/2., nn/2.)/sr

    if pre_filt:
        lqtcopy.filter('bandpass', freqmin=pre_filt[0], 
            freqmax=pre_filt[1], corners=2, zerophase=True)

    if filt:
        rfcopy.filter('bandpass', freqmin=filt[0], 
            freqmax=filt[1], corners=2, zerophase=True)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, figsize=(7,5))

    ax1.plot(taxis, lqtcopy[0], label=lqtcopy[0].stats.channel, lw=1)
    ax2.plot(taxis, lqtcopy[1], label=lqtcopy[0].stats.channel, lw=1)
    ax3.plot(taxis, lqtcopy[2], label=lqtcopy[0].stats.channel, lw=1)

    ax1.legend()
    ax2.legend()
    ax3.legend()

    nn = rfcopy[0].stats.npts
    sr = rfcopy[0].stats.sampling_rate
    taxis = np.arange(-nn/2., nn/2.)/sr

    ax4.plot(taxis, rfcopy[0], label=rfcopy[0].stats.channel, lw=1)
    ax4.set_ylim(-1., 1.)
    ax4.plot(taxis, rfcopy[1], label=rfcopy[0].stats.channel, lw=1)
    ax4.plot(taxis, rfcopy[2], label=rfcopy[0].stats.channel, lw=1)
    if trange:
        ax4.set_xlim(trange[0], trange[1])
    ax4.legend()
    plt.suptitle(
        'AZ corr: {0:.1f}; BAZ: {1:.1f}\n SNR: {2:.1f}; CC: {3:.1f}'.format(
            rfdata.sta.azcorr, rfdata.meta.baz, 
            rfdata.meta.snr, rfdata.meta.cc))

    plt.show()


def event_dist(stream, phase='P', save=False, title=None, form='png'):

    import cartopy.crs as ccrs

    # Specify azimuthal equidistant projection
    aeqd = ccrs.AzimuthalEquidistant(central_longitude=stream[0].stats.stlo,
                                     central_latitude=stream[0].stats.stla,
                                     globe=ccrs.Globe())

    # Extract latitude and longitude
    evlat = []
    evlon = []
    phlist = []
    for tr in stream:
        evlat.append(tr.stats.evla)
        evlon.append(tr.stats.evlo)
        phlist.append(tr.stats.phase)

    # Turn into arrays
    evlat = np.array(evlat)
    evlon = np.array(evlon)
    phlist = np.array(phlist)

    # Select phase to plot
    phP = phlist == 'P'
    phPP = phlist == 'PP'

    # Now plot
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(1, 1, 1, projection=aeqd)
    if np.sum(phP) > 0:
        # ax.scatter(evlon[phP], evlat[phP], c='royalblue', label='P')
        ax.scatter(evlon[phP], evlat[phP], c='royalblue', label='P',
                   transform=ccrs.Geodetic())
    if np.sum(phPP) > 0:
        ind = phase == 'P'
        ax.scatter(evlon[phPP], evlat[phPP], c='coral', label='PP',
                   transform=ccrs.Geodetic())

    ax.scatter(stream[0].stats.stlo, stream[0].stats.stla, c='grey',
               marker='v', transform=ccrs.Geodetic())
    ax.coastlines()

    if title:
        plt.suptitle(title)

    if save:
        plt.savefig('RF_PLOTS/' + stream[0].stats.station +
                    '.' + title + '.event_dist.' + form, format=form)
    else:
        plt.show()

    plt.close()

    # ax.gridlines()
