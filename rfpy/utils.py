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

# -*- coding: utf-8 -*-
import math
from obspy import UTCDateTime, Stream, read
import numpy as np
import copy
import re
from stdb import StDbElement


def get_stkeys(inventory, keys=[]):

    allkeys = []
    station_list = inventory.get_contents()['stations']
    allkeys = [s.split(' ')[0] for s in station_list]

    if len(keys) > 0:
        # Extract key subset

        stkeys = []
        for key in keys:
            # Convert the pattern to a regex pattern
            # Replace '.' with '\.' to match literal dots
            # Replace '*' with '.*' to match any sequence of characters
            # Replace '?' with '.' to match any single character
            pattern = key.replace('.', r'\.').replace('*', '.*').replace('?', '.')

            # Compile the regex pattern
            regex = re.compile(f'^.*{pattern}.*$')

            # Filter allkeys based on the compiled regex
            stkeys.extend([key for key in allkeys if regex.match(key)])

    else:
        stkeys = allkeys

    return stkeys


def inv2stdb(inventory, keys=[]):

    stkeys = get_stkeys(inventory, keys)

    stations = {}
    for key in stkeys:
        net = key.split('.')[0]
        sta = key.split('.')[1]
        cha = '?H?'
        inv = inventory.select(network=net, station=sta, channel=cha)
        seed_id = inv.get_contents()['channels'][0]
        coords = inv.get_coordinates(seed_id)
        channel = seed_id.split('.')[3][0:2]
        location = list(seed_id.split('.')[2])
        if len(location) == 0:
            location.append('--')

        stdb_element = StDbElement(
            station=sta,
            network=net,
            channel=channel,
            location=list(set(location)),
            latitude=coords['latitude'],
            longitude=coords['longitude'],
            elevation=coords['elevation'],
            startdate=inv[0].stations[0].start_date,
            enddate=inv[0].stations[0].end_date
            )
        stations[key] = stdb_element

    return stations, stkeys


def floor_decimal(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier


def traceshift(trace, tt):
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
        ftrace[i] = ftrace[i]*np.exp(-2.*np.pi*1j*freq[i]*tt)

    # Back Fourier transform and return as trace
    rtrace = trace.copy()
    rtrace.data = np.real(np.fft.ifft(ftrace))

    # Update start time
    rtrace.stats.starttime -= tt

    return rtrace


def download_data(client=None, sta=None, start=UTCDateTime(),
                  end=UTCDateTime(), new_sr=0., verbose=False,
                  remove_response=False, zcomp='Z'):
    """
    Function to build a stream object for a seismogram in a given time window
    by getting data from a client object, either from a local SDS archive or
    from an FDSN web-service. The function performs sanity checks for
    the start times, sampling rates and window lengths.

    Parameters
    ----------
    client : :class:`~obspy.client.fdsn.Client`
        Client object
    sta : Dict
        Station metadata from :mod:`~StDb` data base
    start : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time for request
    end : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time for request
    new_sr : float
        New sampling rate (Hz)
    verbose : bool
        Whether or not to print messages to screen during run-time
    remove_response : bool
        Remove instrument response from seismogram and resitute to true ground
        velocity (m/s) using obspy.core.trace.Trace.remove_response()
    zcomp: str
        Vertical Component Identifier. Should be a single character.
        This is different then 'Z' only for fully unknown component
        orientation (i.e., components are 1, 2, 3)

    Returns
    -------
    err : bool
        Boolean for error handling (`False` is associated with success)
    trN : :class:`~obspy.core.Trace`
        Trace of North component of motion
    trE : :class:`~obspy.core.Trace`
        Trace of East component of motion
    trZ : :class:`~obspy.core.Trace`
        Trace of Vertical component of motion

    """

    for loc in sta.location:

        # Construct location name
        if loc == "--":
            tloc = ""
        else:
            tloc = copy.copy(loc)

        # Construct Channel List
        cha = sta.channel.upper() + '?'
        msg = "*     {0:s}.{1:2s}?.{2:2s} - Checking Network".format(
            sta.station, sta.channel.upper(), loc)
        print(msg)

        # Get waveforms, with extra 1 second to avoid
        # traces cropped too short - traces are trimmed later
        try:
            st = client.get_waveforms(
                network=sta.network,
                station=sta.station,
                location=tloc,
                channel=cha,
                starttime=start,
                endtime=end+1.)
        except Exception as e:
            if verbose:
                print("* Met exception:")
                print("* " + e.__repr__())
            st = None
        else:
            if len(st) == 3:
                # It's possible if len(st)==1 that data is Z12
                print("*                 - Data Downloaded")
            elif len(st) < 3:
                print("* Error retrieving waveforms")
                print("**************************************************")
                return True, None

    # Check the correct 3 components exist
    if st is None:
        print("* Error retrieving waveforms")
        print("**************************************************")
        return True, None

    # Three components successfully retrieved
    if remove_response:
        try:
            st.remove_response()
            if verbose:
                print("*")
                print("* Restituted stream to true ground velocity.")
                print("*")
        except Exception as e:
            print("*")
            print('* Cannot remove response, moving on.')
            print("*")

    # Detrend and apply taper
    st.detrend('demean').detrend('linear').taper(
        max_percentage=0.05, max_length=5.)

    # Check start times
    if not np.all([tr.stats.starttime == start for tr in st]):
        if verbose:
            print("* Start times are not all close to true start: ")
            [print("*   "+tr.stats.channel+" " +
                   str(tr.stats.starttime)+" " +
                   str(tr.stats.endtime)) for tr in st]
            print("*   True start: "+str(start))
            print("*   -> Shifting traces to true start")
        delay = [tr.stats.starttime - start for tr in st]
        st_shifted = Stream(
            traces=[traceshift(tr, dt) for tr, dt in zip(st, delay)])
        st = st_shifted.copy()

    # Check sampling rate
    sr = st[0].stats.sampling_rate
    sr_round = float(floor_decimal(sr, 0))
    if not sr == sr_round:
        if verbose:
            print("* Sampling rate is not an integer value: ", sr)
            print("* -> Resampling")
        st.resample(sr_round, no_filter=False)

    # Try trimming
    try:
        st.trim(start, end)
    except Exception as e:
        print("* Unable to trim")
        print("* -> Skipping")
        print("**************************************************")

        return True, None

    # Check final lengths - they should all be equal if start times
    # and sampling rates are all equal and traces have been trimmed
    if not np.allclose([tr.stats.npts for tr in st[1:]], st[0].stats.npts):
        print("* Lengths are incompatible: ")
        [print("*     "+str(tr.stats.npts)) for tr in st]
        print("* -> Skipping")
        print("**************************************************")

        return True, None

    elif not np.allclose([st[0].stats.npts], int((end - start)*sr),
                         atol=1):
        print("* Length is too short: ")
        print("*    "+str(st[0].stats.npts) +
              " ~= "+str(int((end - start)*sr)))
        print("* -> Skipping")
        print("**************************************************")

        return True, None

    else:
        print("* Waveforms Retrieved...")
        return False, st
