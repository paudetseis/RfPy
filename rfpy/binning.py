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
Functions to bin receiver functions along back-azimuth or slowness
bins, or both, or for all data, regardless of the x axis (time or depth).

"""

# Import modules and functions
import numpy as np
from obspy.core import Stream, Trace
from scipy.signal import hilbert


def bin(stream1, stream2=None, typ='baz', nbin=36+1, pws=False):
    """ 
    Function to stack receiver functions into (baz or slow) bins
    This can be done using a linear stack (i.e., simple
    mean), or using phase-weighted stacking.

    Parameters
    ----------
    stream1 : :class:`~obspy.core.Stream`
        Stream of equal-length seismograms to be stacked into
        a single trace.
    stream2 : :class:`~obspy.core.Stream`
        Optionally stack a second stream in the same operation.
    dbaz : int
        Number of bazk-azimuth samples in bins
    dslow : int
        Number of slowness samples in bins
    pws : bool
        Whether or not to perform phase-weighted stacking

    Returns
    -------
    stack : :class:`~obspy.core.Stream`
        Stream containing one or two stacked traces,
        depending on the number of input streams

    """

    if not typ in ['baz', 'slow', 'dist']:
        raise(Exception("Type has to be 'baz' or 'slow' or 'dist'"))

    if typ == 'baz':
        bmin = 0
        bmax = 360
        stat = [stream1[i].stats.baz for i in range(len(stream1))]
    elif typ == 'slow':
        stat = [stream1[i].stats.slow for i in range(len(stream1))]
        bmin = np.min(np.array(stat))
        bmax = np.max(np.array(stat))
    elif typ == 'dist':
        stat = [stream1[i].stats.gac for i in range(len(stream1))]
        bmin = np.min(np.array(stat))
        bmax = np.max(np.array(stat))

    # Define bins
    bins = np.linspace(bmin, bmax, nbin)

    # Digitize stat
    ind = np.digitize(stat, bins)

    final_stream = []

    for stream in [stream1, stream2]:
        try:
            # Define empty streams
            binned_stream = Stream()

            # Loop through bins
            for i in range(nbin):

                nb = 0
                array = np.zeros(len(stream[0].data))
                weight = np.zeros(len(stream[0].data), dtype=complex)

                # Loop through stat
                for j, tr in enumerate(stream):

                    # If index of bins is equal to ind
                    if i == ind[j]:

                        nb += 1
                        array += tr.data
                        hilb = hilbert(tr.data)
                        phase = np.arctan2(hilb.imag, hilb.real)
                        weight += np.exp(1j*phase)
                        
                        continue

                if nb > 0:

                    # Average and update stats
                    array /= nb
                    weight = np.real(abs(weight/np.float(nb)))

                    trace = Trace(header=stream[0].stats)
                    trace.stats.nbin = nb
                    if typ == 'baz':
                        trace.stats.baz = bins[i]
                        trace.stats.slow = None
                        trace.stats.nbin = nb
                    elif typ == 'slow': 
                        trace.stats.slow = bins[i]
                        trace.stats.baz = None
                        trace.stats.nbin = nb
                    elif typ == 'dist':
                        trace.stats.dist = bins[i]
                        trace.stats.slow = None
                        trace.stats.baz = None
                        trace.stats.nbin = nb
                    if not pws:
                        weight = np.ones(len(stream[0].data))
                    trace.data = weight*array
                    binned_stream.append(trace)

            final_stream.append(binned_stream)

        except:
            continue

    return final_stream


def bin_baz_slow(stream1, stream2=None, nbaz=36+1, nslow=20+1, pws=False):
    """ 
    Function to stack receiver functions into back-azimuth and slowness bins.
    This can be done using a linear stack (i.e., simple
    mean), or using phase-weighted stacking.

    Parameters
    ----------
    stream1 : :class:`~obspy.core.Stream`
        Stream of equal-length seismograms to be stacked into
        a single trace.
    stream2 : :class:`~obspy.core.Stream`
        Optionally stack a second stream in the same operation.
    dbaz : int
        Number of bazk-azimuth samples in bins
    dslow : int
        Number of slowness samples in bins
    pws : bool
        Whether or not to perform phase-weighted stacking

    Returns
    -------
    stack : :class:`~obspy.core.Stream`
        Stream containing one or two stacked traces,
        depending on the number of input streams

    """

    # Define back-azimuth and slowness bins
    baz_bins = np.linspace(0, 360, nbaz)
    slow_bins = np.linspace(0.04, 0.08, nslow)

    # Extract baz and slowness
    baz = [stream1[i].stats.baz for i in range(len(stream1))]
    slow = [stream1[i].stats.slow for i in range(len(stream1))]

    # Digitize baz and slowness
    ibaz = np.digitize(baz, baz_bins)
    islow = np.digitize(slow, slow_bins)

    final_stream = []

    for stream in [stream1, stream2]:
        try:
            # Define empty streams
            binned_stream = Stream()

            # Loop through baz_bins
            for i in range(nbaz):
                for j in range(nslow):

                    nbin = 0
                    array = np.zeros(len(stream[0].data))
                    weight = np.zeros(len(stream[0].data), dtype=complex)

                    # Loop all traces
                    for k, tr in enumerate(stream):

                        # If index of baz_bins is equal to ibaz
                        if i == ibaz[k] and j == islow[k]:

                            nbin += 1
                            array += tr.data
                            hilb = hilbert(tr.data)
                            phase = np.arctan2(hilb.imag, hilb.real)
                            weight += np.exp(1j*phase)
                            
                            continue

                    if nbin > 0:

                        # Average and update stats
                        array /= nbin
                        weight = np.real(abs(weight/nbin))

                        trace = Trace(header=stream[0].stats)
                        trace.stats.baz = baz_bins[i]
                        trace.stats.slow = slow_bins[j]
                        trace.stats.nbin = nbin

                        if not pws:
                            weight = np.ones(len(stream[0].data))
                        trace.data = weight*array
                        binned_stream.append(trace)

            final_stream.append(binned_stream)

        except:
            continue

    return final_stream

def bin_all(stream1, stream2=None, pws=False):
    """ 
    Function to bin all streams into a single trace.
    This can be done using a linear stack (i.e., simple
    mean), or using phase-weighted stacking.

    Parameters
    ----------
    stream1 : :class:`~obspy.core.Stream`
        Stream of equal-length seismograms to be stacked into
        a single trace.
    stream2 : :class:`~obspy.core.Stream`
        Optionally stack a second stream in the same operation.
    pws : bool
        Whether or not to perform phase-weighted stacking

    Returns
    -------
    stack : :class:`~obspy.core.Stream`
        Stream containing one or two stacked traces,
        depending on the number of input streams

    """

    # Initialize empty stack stream
    stack = Stream()
    for stream in [stream1, stream2]:
        try:
            # Copy stats from stream1
            stats = stream[0].stats

            # Initialize arrays
            array = np.zeros(len(stream[0].data))
            pweight = np.zeros(len(stream[0].data), dtype=complex)

            # Get phase weights
            for tr in stream:
                array += tr.data
                hilb = hilbert(tr.data)
                phase = np.arctan2(hilb.imag, hilb.real)
                pweight += np.exp(1j*phase)

            # Normalize
            array = array/len(stream)
            weight = np.real(abs(pweight/len(stream)))
            # Regular linear stack
            if not pws:
                weight = np.ones(len(stream[0].data))

            # Put back into traces
            stack.append(Trace(data=weight*array, header=stats))

        except:
            continue

    return stack

