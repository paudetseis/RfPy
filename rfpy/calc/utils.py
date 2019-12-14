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
"""

from obspy.core import UTCDateTime
from numpy import nan


def parse_localdata_for_comp(comp='Z', stdata=list, sta=None,
                             start=UTCDateTime, end=UTCDateTime, ndval=nan):
    """
    Function to determine the path to data for a given component and alternate network

    Parameters
    ----------
    comp : str
        Channel for seismogram (one letter only)
    stdata : List
        Station list
    sta : Dict
        Station metadata from :mod:`~StDb` data base
    start : :class:`~obspy.core.utcdatetime.UTCDateTime`
        Start time for request
    end : :class:`~obspy.core.utcdatetime.UTCDateTime`
        End time for request
    ndval : float or nan
        Default value for missing data

    Returns
    -------
    err : bool
        Boolean for error handling (`False` is associated with success)
    st : :class:`~obspy.core.Stream`
        Stream containing North, East and Vertical components of motion

    """

    from fnmatch import filter
    from obspy import read

    # Get start and end parameters
    styr = start.strftime("%Y")
    stjd = start.strftime("%j")
    edyr = end.strftime("%Y")
    edjd = end.strftime("%j")

    # Intialize to default positive error
    erd = True

    print(
        ("*          {0:2s}{1:1s} - Checking Disk".format(sta.channel.upper(),
                                                          comp.upper())))

    # Time Window Spans Single Day
    if stjd == edjd:
        # Format 1
        lclfiles = list(filter(
            stdata,
            '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.{4:2s}{5:1s}.SAC'.format(
                styr, stjd, sta.network.upper(
                ), sta.station.upper(), sta.channel.upper()[0:2],
                comp.upper())))
        # Format 2
        if len(lclfiles) == 0:
            lclfiles = list(filter(
                stdata,
                '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.*{4:1s}.SAC'.format(
                    styr, stjd, sta.network.upper(), sta.station.upper(),
                    comp.upper())))
        # Alternate Nets (for CN/PO issues) Format 1
        if len(lclfiles) == 0:
            lclfiles = []
            for anet in sta.altnet:
                lclfiles.extend(
                    list(
                        filter(
                            stdata,
                            '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.{4:2s}{5:1s}.SAC'.format(
                                styr, stjd, anet.upper(), sta.station.upper(),
                                sta.channel.upper()[0:2], comp.upper()))))

        # Alternate Nets (for CN/PO issues) Format 2
        if len(lclfiles) == 0:
            # Check Alternate Networks
            lclfiles = []
            for anet in sta.altnet:
                lclfiles.extend(
                    list(
                        filter(
                            stdata,
                            '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.*{4:1s}.SAC'.format(
                                styr, stjd, sta.network.upper(),
                                sta.station.upper(), comp.upper()))))

        # If still no Local files stop
        if len(lclfiles) == 0:
            print("*              - Data Unavailable")
            return erd, None

        # Process the local Files
        for sacfile in lclfiles:
            # Read File
            st = read(sacfile, format="SAC")

            # Should only be one component, otherwise keep reading If more
            # than 1 component, error
            if len(st) != 1:
                pass

            else:
                # Check for NoData and convert to NaN
                stnd = st[0].stats.sac['user9']
                eddt = False
                if (not stnd == 0.0) and (not stnd == -12345.0):
                    st[0].data[st[0].data == stnd] = ndval
                    eddt = True

                # Check start/end times in range
                if (st[0].stats.starttime <= start and
                        st[0].stats.endtime >= end):
                    st.trim(starttime=start, endtime=end)

                    # Check for Nan in stream
                    if True in isnan(st[0].data):
                        print(
                            "*          !!! Missing Data Present !!! " +
                            "Skipping (NaNs)")
                    else:
                        if eddt and (ndval == 0.0):
                            if any(st[0].data == 0.0):
                                print(
                                    "*          !!! Missing Data Present " +
                                    "!!! (Set to Zero)")

                        st[0].stats.update()
                        tloc = st[0].stats.location
                        if len(tloc) == 0:
                            tloc = "--"

                        # Processed succesfully...Finish
                        print(("*          {1:3s}.{2:2s}  - From Disk".format(
                            st[0].stats.station, st[0].stats.channel.upper(),
                            tloc)))
                        return False, st

    # Time Window spans Multiple days
    else:
        # Day 1 Format 1
        lclfiles1 = list(
            filter(stdata,
                   '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.{4:2s}{5:1s}.SAC'.format(
                       styr, stjd, sta.network.upper(), sta.station.upper(),
                       sta.channel.upper()[0:2], comp.upper())))
        # Day 1 Format 2
        if len(lclfiles1) == 0:
            lclfiles1 = list(
                filter(stdata,
                       '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.*{4:1s}.SAC'.format(
                           styr, stjd, sta.network.upper(),
                           sta.station.upper(), comp.upper())))
        # Day 1 Alternate Nets (for CN/PO issues) Format 1
        if len(lclfiles1) == 0:
            lclfiles1 = []
            for anet in sta.altnet:
                lclfiles1.extend(
                    list(
                        filter(
                            stdata,
                            '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.{4:2s}{5:1s}.SAC'.format(
                                styr, stjd, anet.upper(), sta.station.upper(
                                ), sta.channel.upper()[0:2],
                                comp.upper()))))
        # Day 1 Alternate Nets (for CN/PO issues) Format 2
        if len(lclfiles1) == 0:
            lclfiles1 = []
            for anet in sta.altnet:
                lclfiles1.extend(
                    list(
                        filter(
                            stdata,
                            '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.*{4:1s}.SAC'.format(
                                styr, stjd, anet.upper(),
                                sta.station.upper(), comp.upper()))))

        # Day 2 Format 1
        lclfiles2 = list(
            filter(stdata,
                   '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.{4:2s}{5:1s}.SAC'.format(
                       edyr, edjd, sta.network.upper(
                       ), sta.station.upper(), sta.channel.upper()[0:2],
                       comp.upper())))
        # Day 2 Format 2
        if len(lclfiles2) == 0:
            lclfiles2 = list(
                filter(stdata,
                       '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.*{4:1s}.SAC'.format(
                           edyr, edjd, sta.network.upper(), sta.station.upper(),
                           comp.upper())))
        # Day 2 Alternate Nets (for CN/PO issues) Format 1
        if len(lclfiles2) == 0:
            lclfiles2 = []
            for anet in sta.altnet:
                lclfiles2.extend(
                    list(
                        filter(
                            stdata,
                            '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.{4:2s}{5:1s}.SAC'.format(
                                edyr, edjd, anet.upper(), sta.station.upper(),
                                sta.channel.upper()[0:2], comp.upper()))))
        # Day 2 Alternate Nets (for CN/PO issues) Format 2
        if len(lclfiles2) == 0:
            lclfiles2 = []
            for anet in sta.altnet:
                lclfiles2.extend(
                    list(
                        filter(
                            stdata,
                            '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.*{4:1s}.SAC'.format(
                                edyr, edjd, anet.upper(), sta.station.upper(),
                                comp.upper()))))

        # If still no Local files stop
        if len(lclfiles1) == 0 and len(lclfiles2) == 0:
            print("*              - Data Unavailable")
            return erd, None

        # Now try to merge the two separate day files
        if len(lclfiles1) > 0 and len(lclfiles2) > 0:
            # Loop over first day file options
            for sacf1 in lclfiles1:
                st1 = read(sacf1, format='SAC')
                # Loop over second day file options
                for sacf2 in lclfiles2:
                    st2 = read(sacf2, format='SAC')

                    # Check time overlap of the two files.
                    if st1[0].stats.endtime >= \
                            st2[0].stats.starttime-st2[0].stats.delta:
                        # Check for NoData and convert to NaN
                        st1nd = st1[0].stats.sac['user9']
                        st2nd = st2[0].stats.sac['user9']
                        eddt1 = False
                        eddt2 = False
                        if (not st1nd == 0.0) and (not st1nd == -12345.0):
                            st1[0].data[st1[0].data == st1nd] = ndval
                            eddt1 = True
                        if (not st2nd == 0.0) and (not st2nd == -12345.0):
                            st2[0].data[st2[0].data == st2nd] = ndval
                            eddt2 = True

                        st = st1 + st2
                        # Need to work on this HERE (AJS OCT 2015).
                        # If Calibration factors are different,
                        try:
                                # then the traces cannot be merged.
                            st.merge()

                            # Should only be one component, otherwise keep
                            # reading If more than 1 component, error
                            if len(st) != 1:
                                print(st)
                                print("merge failed?")

                            else:
                                if (st[0].stats.starttime <= start and
                                        st[0].stats.endtime >= end):
                                    st.trim(starttime=start, endtime=end)

                                    # Check for Nan in stream
                                    if True in isnan(st[0].data):
                                        print(
                                            "*          !!! Missing Data " +
                                            "Present !!! Skipping (NaNs)")
                                    else:
                                        if (eddt1 or eddt2) and (ndval == 0.0):
                                            if any(st[0].data == 0.0):
                                                print(
                                                    "*          !!! Missing " +
                                                    "Data Present !!! (Set " +
                                                    "to Zero)")

                                        st[0].stats.update()
                                        tloc = st[0].stats.location
                                        if len(tloc) == 0:
                                            tloc = "--"

                                        # Processed succesfully...Finish
                                        print(("*          {1:3s}.{2:2s}  - " +
                                               "From Disk".format(
                                                   st[0].stats.station,
                                                   st[0].stats.channel.upper(),
                                                   tloc)))
                                        return False, st

                        except:
                            pass
                    else:
                        print(("*                 - Merge Failed: No " +
                               "Overlap {0:s} - {1:s}".format(
                                   st1[0].stats.endtime,
                                   st2[0].stats.starttime -
                                   st2[0].stats.delta)))

    # If we got here, we did not get the data.
    print("*              - Data Unavailable")
    return erd, None


def get_data_NEZ(client=None, sta=None, start=UTCDateTime, end=UTCDateTime,
                 stdata=list, ndval=nan, new_sr=0.):
    """
    Function to build a stream object for a seismogram in a given time window either
    by downloading data from the client object or alternatively first checking if the
    given data is already available locally.

    Note 
    ----
    Currently only supports NEZ Components!

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
    stdata : List
        Station list
    ndval : float or nan
        Default value for missing data

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

    from fnmatch import filter
    from obspy import read
    from os.path import dirname, join, exists
    from numpy import any

    # Output
    print(("*     {0:s}.{1:2s} - NEZ:".format(sta.station,
                                              sta.channel.upper())))

    # Set Error Default to True
    erd = True

    # Check if there is local data
    if len(stdata) > 0:
        # Only a single day: Search for local data
        # Get Z localdata
        errZ, stZ = parse_localdata_for_comp(
            comp='Z', stdata=stdata, sta=sta, start=start, end=end,
            ndval=ndval)
        # Get N localdata
        errN, stN = parse_localdata_for_comp(
            comp='N', stdata=stdata, sta=sta, start=start, end=end,
            ndval=ndval)
        # Get E localdata
        errE, stE = parse_localdata_for_comp(
            comp='E', stdata=stdata, sta=sta, start=start, end=end,
            ndval=ndval)
        # Retreived Succesfully?
        erd = errZ or errN or errE
        if not erd:
            # Combine Data
            st = stZ + stN + stE

    # No local data? Request using client
    if erd:
        erd = False

        for loc in sta.location:
            tloc = loc
            # Constract location name
            if len(tloc) == 0:
                tloc = "--"
            # Construct Channel List
            channels = sta.channel.upper() + 'E,' + sta.channel.upper() + \
                'N,' + sta.channel.upper() + 'Z'
            print(("*          {1:2s}[ENZ].{2:2s} - Checking Network".format(
                sta.station, sta.channel.upper(), tloc)))
            try:
                st = client.get_waveforms(network=sta.network,
                                          station=sta.station, location=loc,
                                          channel=channels, starttime=start,
                                          endtime=end, attach_response=False)
                if len(st) == 3:
                    erd = False
                    print("*             - Data Downloaded")
                else:
                    erd = True
                    st = None
                    print("*             - Data Unavailable")
            except:
                erd = True
                st = None
                print("*             - Data Unavailable")

            # Break if we successfully obtained 3 components in st
            if not erd:
                break

    # Check the correct 3 components exist
    if st is None:
        print("* Error retrieving waveforms")
        print("****************************************************")
        return True, None, None, None
    elif (not st.select(component='Z')[0] or not st.select(component='E'[0])
          or not st.select(component='N')[0]):
        print("* Error retrieving waveforms")
        if not st.select(component='Z')[0]:
            print("*   --Z Missing")
        if not st.select(component='E')[0]:
            print("*   --E Missing")
        if not st.select(component='N')[0]:
            print("*   --N Missing")
        print("****************************************************")
        return True, None, None, None

    # Three components successfully retrieved
    print("* Waveforms Retrieved...")

    # Check trace lengths
    ll0 = len(st.select(component='Z')[0].data)
    ll1 = len(st.select(component='E')[0].data)
    ll2 = len(st.select(component='Z')[0].data)
    if not (ll0 == ll1 and ll0 == ll2):
        print(("* Error:  Trace lengths (Z,E,N): ", ll0, ll1, ll2))
        print("****************************************************")
        return True, None, None, None

    # Detrend data
    st.detrend('demean')
    st.detrend('linear')

    # Filter Traces
    st.filter('lowpass', freq=0.5*new_sr, corners=2, zerophase=True)
    st.resample(new_sr)

    # Extract traces
    trE = st.select(component='E')[0]
    trN = st.select(component='N')[0]
    trZ = st.select(component='Z')[0]

    # Return Flag and Data
    return False, trN, trE, trZ
