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


def bin(stream1, stream2=None, typ='baz', nbin=36+1, pws=False,
        include_empty=False):
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
    typ: str
        Attribute to bin
        'baz': backazimuth (degree)
        'slow': Horizontal slowness (s/km)
        'dist': epicentral distance (degree)
    nbin : int
        Number of equally sized bins within data range
    pws : bool
        Whether or not to perform phase-weighted stacking
    include_empty : bool
        Return empty bins as null traces (default omits them)

    Returns
    -------
    stack : :class:`~obspy.core.Stream`
        Stream containing one or two stacked traces,
        depending on the number of input streams

    Note
    ----
    Sets the following attributes of the stack:
        nbin: Number of bins
        index_list: Indices of constituent traces in source stream
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
                index_list = []

                # Loop through stat
                for j, tr in enumerate(stream):

                    # If index of bins is equal to ind
                    if i == ind[j]:

                        nb += 1
                        array += tr.data
                        index_list.append(j)

                        if pws:
                            hilb = hilbert(tr.data)
                            phase = np.arctan2(hilb.imag, hilb.real)
                            weight += np.exp(1j*phase)
                        
                        continue

                if nb > 0 or include_empty:

                    # Average and update stats
                    array /= nb
                    weight = np.real(abs(weight/np.float(nb)))

                    trace = Trace(header=stream[0].stats)
                    trace.stats.nbin = nb
                    trace.stats.index_list = index_list

                    if typ == 'baz':
                        trace.stats.baz = bins[i]
                        trace.stats.slow = None
                    elif typ == 'slow': 
                        trace.stats.slow = bins[i]
                        trace.stats.baz = None
                    elif typ == 'dist':
                        trace.stats.dist = bins[i]
                        trace.stats.slow = None
                        trace.stats.baz = None
                    if not pws:
                        weight = np.ones(len(stream[0].data))
                    trace.data = weight*array
                    binned_stream.append(trace)

            final_stream.append(binned_stream)

        except:
            continue

    return final_stream


def bin_baz_slow(stream1, stream2=None, nbaz=36+1, nslow=20+1, baz_range=(0, 360),
        slow_range=(0.04, 0.08), pws=False, include_empty=False):
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
    nbaz : int
        Number of back-azimuth samples in bins
    nslow : int
        Number of slowness samples in bins
    baz_range : 2-tuple
        Minimum and maximum of the back-azimuth bins
    slow_range : 2-tuple
        Minimum and maximum of the slowness bins
    pws : bool
        Whether or not to perform phase-weighted stacking
    include_empty : bool
        Return empty bins as null traces (default omits them)

    Returns
    -------
    stack : :class:`~obspy.core.Stream`
        Stream containing one or two stacked traces,
        depending on the number of input streams

    Note
    ----
    Sets the following attributes of the stack:
        nbin: Number of bins
        index_list: Indices of constituent traces in source stream
    """

    # Define back-azimuth and slowness bins
    baz_bins = np.linspace(*baz_range, nbaz)
    slow_bins = np.linspace(*slow_range, nslow)

    # Extract baz and slowness
    baz = [stream1[i].stats.baz for i in range(len(stream1))]
    slow = [stream1[i].stats.slow for i in range(len(stream1))]

    # Digitize baz and slowness
    ibaz = np.digitize(baz, baz_bins)
    islow = np.digitize(slow, slow_bins)

    final_stream = []
    if stream2:
        stream1 = [stream1, stream2]
    else:
        stream1 = [stream1]

    for stream in stream1:
        # try:
        # Define empty streams
        binned_stream = Stream()

        # Loop through baz_bins
        for i in range(nbaz):
            for j in range(nslow):

                nbin = 0
                index_list = []
                array = np.zeros(len(stream[0].data), dtype=type(stream[0].data[0]))
                weight = np.zeros(len(stream[0].data), dtype=complex)

                # Loop all traces
                for k, tr in enumerate(stream):

                    # If index of baz_bins is equal to ibaz
                    if i == ibaz[k] and j == islow[k]:

                        nbin += 1
                        array += tr.data
                        index_list.append(k)

                        if pws:
                            hilb = hilbert(np.real(tr.data))
                            phase = np.arctan2(hilb.imag, hilb.real)
                            weight += np.exp(1j*phase)
                        
                        continue

                if nbin > 0 or include_empty:

                    # Average and update stats
                    array /= nbin
                    weight = np.real(abs(weight/nbin))

                    trace = Trace(header=stream[0].stats)
                    trace.stats.baz = baz_bins[i]
                    trace.stats.slow = slow_bins[j]
                    trace.stats.nbin = nbin
                    trace.stats.index_list = index_list

                    if not pws:
                        weight = np.ones(len(stream[0].data))
                    trace.data = weight*array
                    binned_stream.append(trace)

        final_stream.append(binned_stream)

        # except:
        #     final_stream = [Stream()]
        #     continue

    return final_stream


def stack(rfdatas, which='rf', pws=False,
        attributes={'bin': None, 'slow': None, 'index_list': None}):
    """ 
    Function to stack receiver functions stored in a RFData list. This can
    be done using a linear stack (i.e., simple mean), or using phase-weighted
    stacking.

    Parameters
    ----------
    rfdatas : list of :class:`rfpy.rfdata.RFdata`
        List of RFData objects to be stacked
    which : str
        'rf' stack receiver functions
        'specs' stack spectra
    pws : bool
        Whether or not to perform phase-weighted stacking (which=rf only)
    attributes : dict
        Attributes passed to the traces of stack

    Returns
    -------
    stack : :class:`~obspy.core.Stream`
        Stream containing stacked traces
    """

    nbin = len(rfdatas)

    template = rfdatas[0].__dict__[which]
    typ = float
    if which == 'specs':
        typ = complex

    array = np.zeros((len(template), len(template[0].data)), dtype=typ)
    weight = np.zeros((len(template), len(template[0].data)), dtype=complex)

    for rfdata in rfdatas:
        stream = rfdata.__dict__[which]

        for n, tr in enumerate(stream):
            array[n, :] += tr.data

            if pws:
                hilb = hilbert(np.real(tr.data))
                phase = np.arctan2(hilb.imag, hilb.real)
                weight[n, :] += np.exp(1j*phase)

    if pws:
        weight = np.real(abs(weight/nbin))
        array *= weight
    array /= nbin

    stacks = Stream()
    for n in range(len(template)):
        stack = Trace(header=template[n].stats)
        stack.data = array[n, :]
        stack.stats.nbin = nbin
        for att in attributes:
            val = attributes[att]
            stack.stats.__dict__[att] = val
        stacks.append(stack)

    return stacks


def baz_slow_bin_indices(
        rfs, nbaz=36+1, nslow=20+1, baz_range=(0, 360), slow_range=(0.04, 0.08), include_empty=False):
    """ 
    Sort traces of streams into backazimuth (nbaz) and slowness (nslow) bins.

    Parameters
    ----------
    rfs : list of :class:`~rfpy.RFData`
        List of receiver functions to be sorted in bins
    bazs : int
        Number of back-azimuth samples in bins
    slows : int
        Number of slowness samples in bins
    baz_range : 2-tuple
        Minimum and maximum of the back-azimuth bins
    slow_range : 2-tuple
        Minimum and maximum of the slowness bins
    include_empty : bool
        Include empty bins (default omits them)

    Returns
    -------
    index_lists : list of lists
        Indices that sorts traces of stream into bins defined by nbaz and nslow
    bazslow_list : list of 2*tuple
        Backazimuth and slowness values of each element of index list
    """

    # Define back-azimuth and slowness bins
    baz_bins = np.linspace(*baz_range, nbaz)
    slow_bins = np.linspace(*slow_range, nslow)

    # Extract baz and slowness
    baz = [rf.meta.baz for rf in rfs]
    slow = [rf.meta.slow for rf in rfs]

    # Digitize baz and slowness
    ibaz = np.digitize(baz, baz_bins)
    islow = np.digitize(slow, slow_bins)

    index_lists = []
    bazslow_list = []

    for i in range(nbaz):
        for j in range(nslow):
            index_list = []

            for k in range(len(rfs)):
                if i == ibaz[k] and j == islow[k]:
                    index_list.append(k)

            if len(index_list) > 0 or include_empty:
                index_lists.append(index_list)
                bazslow_list.append((baz_bins[i], slow_bins[j]))

    return index_lists, bazslow_list


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
                if pws:
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

