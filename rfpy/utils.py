import math
from obspy import UTCDateTime
from numpy import nan, isnan, abs
import numpy as np
from obspy import Stream, Inventory
from obspy import read, read_inventory


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


def list_local_data_stn(lcldrs=list, sta=None, net=None, dtype='SAC', altnet=[]):
    """
    Function to take the list of local directories and recursively
    find all data that matches the station name

    Parameters
    ----------
    lcldrs : List
        List of local directories
    sta : Dict
        Station metadata from :mod:`~StDb`
    net : str
        Network name
    altnet : List
        List of alternative networks

    Returns
    -------
    fpathmatch : List
        Sorted list of matched directories

    """
    from fnmatch import filter
    from os import walk
    from os.path import join

    if sta is None:
        return []
    else:
        if net is None:
            sstrings = ['*.{0:s}.*.{1:s}'.format(sta, dtype)]
        else:
            sstrings = ['*.{0:s}.{1:s}.*.{2:s}'.format(net, sta, dtype)]
            if len(altnet) > 0:
                for anet in altnet:
                    sstrings.append(
                        '*.{0:s}.{1:s}.*.{2:s}'.format(anet, sta, dtype))

    fpathmatch = []
    # Loop over all local data directories
    for lcldr in lcldrs:
        if not lcldr:
            continue
        # Recursiely walk through directory
        for root, dirnames, filenames in walk(lcldr):
            # Keep paths only for those matching the station
            for sstring in sstrings:
                for filename in filter(filenames, sstring):
                    # No hidden directories
                    if not '/.' in root:
                        fpathmatch.append(join(root, filename))

    fpathmatch.sort()

    return fpathmatch


def parse_localdata_for_comp(comp='Z', stdata=[], dtype='SAC', sta=None,
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

    # Possible naming conventions
    f1 = '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.{4:2s}{5:1s}.{6:s}'
    f2 = '*/{0:4s}.{1:3s}.{2:s}.{3:s}.*.*{4:1s}.{5:s}'
    f3 = '*/{0:4s}.{1:3s}.*.{2:s}.{3:s}.*.{4:2s}{5:1s}*.{6:s}'
    f4 = '*/{0:4s}.{1:3s}.*.{2:s}.{3:s}.*.*{4:1s}.D.{5:s}'
    f5 = '*/{0:4s}.{1:3s}.*.??.{2:s}.*.*{3:1s}.D.{4:s}'
    f6 = '*/{0:}/{1:}/{2:}/{3:}{4:}.D/{1:}.{2:}.*.*.D.{0:}.{5:}'

    # Time Window Spans Single Day
    if stjd == edjd:

        lclfiles = []
        nets = [sta.network] + sta.altnet
        for net in nets:
            # Start day
            s1 = f1.format(styr, stjd, net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), dtype)
            s2 = f2.format(styr, stjd, net.upper(), sta.station.upper(),
                           comp.upper(), dtype)
            s3 = f3.format(styr, stjd, net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), dtype)
            s4 = f4.format(styr, stjd, net.upper(), sta.station.upper(),
                           comp.upper(), dtype)
            s5 = f5.format(styr, stjd, sta.station.upper(),
                           comp.upper(), dtype)
            s6 = f6.format(styr,net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), stjd)

            print("*          Trying formats:")
            print("*          " + s1)
            print("*          " + s2)
            print("*          " + s3)
            print("*          " + s4)
            print("*          " + s5)
            print("*          " + s6)
            print("*          ")

            lclfiles.extend(list(filter(stdata, s1)))
            lclfiles.extend(list(filter(stdata, s2)))
            lclfiles.extend(list(filter(stdata, s3)))
            lclfiles.extend(list(filter(stdata, s4)))
            lclfiles.extend(list(filter(stdata, s5)))
            lclfiles.extend(list(filter(stdata, s6)))

        # If still no Local files stop
        if len(lclfiles) == 0:
            print("*              - Data Unavailable")
            return erd, None

        # Process the local Files
        for sacfile in lclfiles:
            # Read File
            try:
                st = read(sacfile)
            except OSError:
                print("*              - Met OSError.")
                print(f"*              - Possibly corrupt file: {sacfile}.")
                return False, Stream()

            if dtype.upper() == 'MSEED':
                if len(st) > 1:
                    st.merge(method=1, interpolation_samples=-
                             1, fill_value=-123456789)

            # Should only be one component, otherwise keep reading If more
            # than 1 component, error
            if len(st) != 1:
                pass

            else:
                # Check start/end times in range
                if (st[0].stats.starttime <= start and
                        st[0].stats.endtime >= end):
                    st.trim(starttime=start, endtime=end)

                    eddt = False
                    # Check for NoData and convert to NaN if a SAC file
                    if dtype.upper() == 'SAC':
                        try:
                            stnd = st[0].stats.sac['user9']
                        except KeyError:
                            stnd = 0.0
                        if (not stnd == 0.0) and (not stnd == -12345.0):
                            st[0].data[st[0].data == stnd] = ndval
                            eddt = True

                    # Check for Nan in stream for SAC
                    if True in isnan(st[0].data):
                        print(
                            "*          !!! Missing Data Present !!! " +
                            "Skipping (NaNs)")
                    # Check for ND Val in stream for MSEED
                    elif -123456789 in st[0].data:
                        print(
                            "*          !!! Missing Data Present !!! " +
                            "Skipping (MSEED fill)")
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
        lclfiles1 = []
        lclfiles2 = []
        nets = [sta.network] + sta.altnet
        for net in nets:
            # Start day
            s1 = f1.format(styr, stjd, net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), dtype)
            s2 = f2.format(styr, stjd, net.upper(), sta.station.upper(),
                           comp.upper(), dtype)
            s3 = f3.format(styr, stjd, net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), dtype)
            s4 = f4.format(styr, stjd, net.upper(), sta.station.upper(),
                           comp.upper(), dtype)
            s5 = f5.format(styr, stjd, sta.station.upper(),
                           comp.upper(), dtype)
            s6 = f6.format(styr,net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), stjd)

            lclfiles1.extend(list(filter(stdata, s1)))
            lclfiles1.extend(list(filter(stdata, s2)))
            lclfiles1.extend(list(filter(stdata, s3)))
            lclfiles1.extend(list(filter(stdata, s4)))
            lclfiles1.extend(list(filter(stdata, s5)))
            lclfiles1.extend(list(filter(stdata, s6)))

            # End day
            s1 = f1.format(edyr, edjd, net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), dtype)
            s2 = f2.format(edyr, edjd, net.upper(), sta.station.upper(),
                           comp.upper(), dtype)
            s3 = f3.format(edyr, edjd, net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), dtype)
            s4 = f4.format(edyr, edjd, net.upper(), sta.station.upper(),
                           comp.upper(), dtype)
            s5 = f5.format(edyr, edjd, sta.station.upper(),
                           comp.upper(), dtype)
            s6 = f6.format(edyr, net.upper(), sta.station.upper(),
                           sta.channel.upper()[0:2], comp.upper(), edjd)

            lclfiles2.extend(list(filter(stdata, s1)))
            lclfiles2.extend(list(filter(stdata, s2)))
            lclfiles2.extend(list(filter(stdata, s3)))
            lclfiles2.extend(list(filter(stdata, s4)))
            lclfiles2.extend(list(filter(stdata, s5)))
            lclfiles2.extend(list(filter(stdata, s6)))

        # If still no Local files stop
        if len(lclfiles1) == 0 and len(lclfiles2) == 0:
            print("*              - Data Unavailable")
            return erd, None

        # Now try to merge the two separate day files
        if len(lclfiles1) > 0 and len(lclfiles2) > 0:
            # Loop over first day file options
            for sacf1 in lclfiles1:
                st1 = read(sacf1)
                if dtype.upper() == 'MSEED':
                    if len(st1) > 1:
                        st1.merge(method=1, interpolation_samples=-1,
                                  fill_value=-123456789)

                # Loop over second day file options
                for sacf2 in lclfiles2:
                    st2 = read(sacf2)
                    if dtype.upper() == 'MSEED':
                        if len(st2) > 1:
                            st2.merge(
                                method=1, interpolation_samples=-1,
                                fill_value=-123456789)

                    # Check time overlap of the two files.
                    if st1[0].stats.endtime >= \
                            st2[0].stats.starttime-st2[0].stats.delta:
                        # eddt1 = False
                        # eddt2 = False
                        # if dtype.upper() == 'SAC':
                        #     # Check for NoData and convert to NaN
                        #     st1nd = st1[0].stats.sac['user9']
                        #     st2nd = st2[0].stats.sac['user9']
                        #     if (not st1nd == 0.0) and (not st1nd == -12345.0):
                        #         st1[0].data[st1[0].data == st1nd] = ndval
                        #         eddt1 = True
                        #     if (not st2nd == 0.0) and (not st2nd == -12345.0):
                        #         st2[0].data[st2[0].data == st2nd] = ndval
                        #         eddt2 = True

                        st = st1 + st2
                        # Need to work on this HERE (AJS OCT 2015).
                        # If Calibration factors are different,
                        # then the traces cannot be merged.
                        try:
                            st.merge(method=1, interpolation_samples=-
                                     1, fill_value=-123456789)

                            # Should only be one component, otherwise keep
                            # reading If more than 1 component, error
                            if len(st) != 1:
                                print(st)
                                print("merge failed?")

                            else:
                                if (st[0].stats.starttime <= start and
                                        st[0].stats.endtime >= end):
                                    st.trim(starttime=start, endtime=end)

                                    eddt = False
                                    # Check for NoData and convert to NaN if a SAC file
                                    if dtype.upper() == 'SAC':
                                        try:
                                            stnd = st[0].stats.sac['user9']
                                        except KeyError:
                                            stnd = 0.0
                                        if (not stnd == 0.0) and (not stnd == -12345.0):
                                            st[0].data[st[0].data == stnd] = ndval
                                            eddt = True

                                    # Check for Nan in stream for SAC
                                    if True in isnan(st[0].data):
                                        print(
                                            "*          !!! Missing Data " +
                                            "Present !!! Skipping (NaNs)")
                                    # Check for ND Val in stream for MSEED
                                    elif -123456789 in st[0].data:
                                        print(
                                            "*          !!! Missing Data Present !!! " +
                                            "Skipping (MSEED fill)")
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
                        st2ot = st2[0].stats.endtime-st2[0].stats.delta
                        print("*                 - Merge Failed: No " +
                              "Overlap {0:s} - {1:s}".format(
                                  st1[0].stats.endtime.strftime(
                                      "%Y-%m-%d %H:%M:%S"),
                                  st2ot.strftime("%Y-%m-%d %H:%M:%S")))

    # If we got here, we did not get the data.
    print("*              - Data Unavailable")
    return erd, None


def attach_local_response_data(stream, local_response_dir):
    """
    Function to restrieve response data from locally stored dataless seed

    Parameters
    ----------
    stream : obspy.Stream
        Event seismogram
    local_response_dir: str
        Directory holding response information. All files containing the station
        name are read. If that fails, all files beginning with the network code are
        read.
    """
    stations = set([t.stats.station for t in stream.traces])
    networks = set([t.stats.network for t in stream.traces])
    inventory = Inventory()
    try:
        for station in stations:
            inventory += read_inventory(
                '{:}/*{:}*'.format(local_response_dir, station))
    except Exception:
        for net in networks:
            inventory += read_inventory('{:}/{:}*'.format(local_response_dir, net))

    stream.attach_response(inventories=inventory)


def download_data(client=None, sta=None, start=UTCDateTime, end=UTCDateTime,
                  stdata=[], dtype='SAC', ndval=nan, new_sr=0., verbose=False,
                  remove_response=False, local_response_dir=''):
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
        List of directories holding local waveform data
    ndval : float or nan
        Default value for missing data
    remove_response : bool
        Remove instrument response from seismogram and resitute to true ground
        velocity (m/s) using obspy.core.trace.Trace.remove_response()
    local_response_dir: str
        Directory holding response files to be read by obspy.read_inventory().
        Required when ``remove_response`` and using locally stored waveform data
        via ``stdata``.

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
    from obspy import read, Stream
    from os.path import dirname, join, exists
    from numpy import any
    from math import floor

    # Output
    print(("*     {0:s}.{1:2s} - ZNE:".format(sta.station,
                                              sta.channel.upper())))

    # Set Error Default to True
    erd = True

    # Check if there is local data
    if len(stdata) > 0:
        # Only a single day: Search for local data
        # Get Z localdata
        errZ, stZ = parse_localdata_for_comp(
            comp='Z', stdata=stdata, dtype=dtype, sta=sta, start=start, end=end,
            ndval=ndval)
        # Get N localdata
        errN, stN = parse_localdata_for_comp(
            comp='N', stdata=stdata, dtype=dtype, sta=sta, start=start, end=end,
            ndval=ndval)
        # Get E localdata
        errE, stE = parse_localdata_for_comp(
            comp='E', stdata=stdata, dtype=dtype, sta=sta, start=start, end=end,
            ndval=ndval)
        # Retreived Succesfully?
        erd = errZ or errN or errE
        if not erd:
            # Combine Data
            st = stZ + stN + stE
            if remove_response:
                if local_response_dir:
                    attach_local_response_data(st, local_response_dir)
                else:
                    print('*')
                    print('* Warning: No local_response_dir given.')
                    print('* Warning: Continuing without.')
                    print('*')


    # No local data? Request using client
    if erd:
        erd = False

        for loc in sta.location:
            tloc = loc
            # Construct location name
            if len(tloc) == 0:
                tloc = "--"

            # Construct Channel List
            chaZNE = (sta.channel.upper() + 'Z,' +
                      sta.channel.upper() + 'N,' +
                      sta.channel.upper() + 'E')
            msgZNE = "*          {1:2s}[ZNE].{2:2s} - Checking Network".format(
                sta.station, sta.channel.upper(), tloc)

            chaZ12 = (sta.channel.upper() + 'Z,' +
                      sta.channel.upper() + '1,' +
                      sta.channel.upper() + '2')
            msgZ12 = "*          {1:2s}[Z12].{2:2s} - Checking Network".format(
                sta.station, sta.channel.upper(), tloc)

            # Loop over possible channels
            for channel, msg in zip([chaZNE, chaZ12], [msgZNE, msgZ12]):
                print(msg)

                kwargs = {}
                if remove_response:
                    kwargs = dict(attach_response=remove_response)

                # Get waveforms, with extra 1 second to avoid
                # traces cropped too short - traces are trimmed later
                try:
                    st = client.get_waveforms(
                        network=sta.network,
                        station=sta.station, location=loc,
                        channel=channel, starttime=start,
                        endtime=end+1., **kwargs)
                except Exception as e:
                    if verbose:
                        print("* Met exception:")
                        print("* " + e.__repr__())
                    st = None
                else:
                    if len(st) == 3:
                        # It's possible if len(st)==1 that data is Z12
                        print("*              - ZNE Data Downloaded")
                        break

            # Break if we successfully obtained 3 components in st
            if not erd:
                break

    # Check the correct 3 components exist
    if st is None or len(st) < 3:
        print("* Error retrieving waveforms")
        print("**************************************************")
        return True, None

    # Three components successfully retrieved
    if remove_response:
        st.remove_response()
        print("*")
        print("* Restituted stream to true ground velocity.")
        print("*")

    # Detrend and apply taper
    st.detrend('demean').detrend('linear').taper(
        max_percentage=0.05, max_length=5.)

    # Check start times
    if not np.all([tr.stats.starttime == start for tr in st]):
        print("* Start times are not all close to true start: ")
        [print("*   "+tr.stats.channel+" " +
               str(tr.stats.starttime)+" " +
               str(tr.stats.endtime)) for tr in st]
        print("*   True start: "+str(start))
        print("* -> Shifting traces to true start")
        delay = [tr.stats.starttime - start for tr in st]
        st_shifted = Stream(
            traces=[traceshift(tr, dt) for tr, dt in zip(st, delay)])
        st = st_shifted.copy()

    # Check sampling rate
    sr = st[0].stats.sampling_rate
    sr_round = float(floor_decimal(sr, 0))
    if not sr == sr_round:
        print("* Sampling rate is not an integer value: ", sr)
        print("* -> Resampling")
        st.resample(sr_round, no_filter=False)

    # Try trimming
    try:
        st.trim(start, end)
    except:
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


def optimal_wlevel(rfdatas, wlevels=10**(np.arange(-5, 0.25, 0.25)), iplot=False):
    """
    Optimal water level for 'water2' deconvolution using Generalized Cross
    Validation

    Accepts
    -------

    rfdatas : list of rfpy.RFData
        rfdata elements hold rotated seismograms in rfdata.data attribute
    wlevels : iterable
        Range of water levels (fraction of maximum stacked amplitude) to try.
        Default: 1e-5 to 10

    Returns:
    --------
    wlevel : float
        Optimal water level (fraction of maximum stacked amplitude)
    """

    nt = len(rfdatas)
    nf = len(rfdatas[0].specs[0].data)

    wft = np.zeros((nt, nf), dtype=complex)
    vft = np.zeros((nt, nf), dtype=complex)

    for it, rfdata in enumerate(rfdatas):

        if not rfdata.meta.rotated:
            msg = 'Element {:} of rfdatas must be rotated.'.format(it)
            raise ValueError(msg)

        wft[it, :] = np.fft.fft(rfdata.data[0].data, n=nf)
        vft[it, :] = np.fft.fft(rfdata.data[1].data, n=nf)

    # Auto and cross spectra
    wwft = wft*np.conjugate(wft)
    vwft = vft*np.conjugate(wft)

    if nt > 0:
        wwft = np.sum(wwft)
        vwft = np.sum(vwft)

    betas = wlevels*np.amax(wwft)
    misfits = np.zeros_like(betas)
    modnorm = np.zeros_like(betas)
    gcvf = np.zeros_like(betas)

    for ib, beta in enumerate(betas):
        # Define operator W W* / (W W* + B) and deconvolve to get
        # impulse response in frequency domain.
        wwft2 = wwft + beta
        rft = vwft / wwft2
        xft = wwft / wwft2

        # Compute model norm.
        modnorm[ib] = np.linalg.norm(rft)**2

        # Compute data misfit. Note misfit is numerator of GCV function
        # Note also earlier mistake where norm(nft)^2 was norm(nft).
        if nt == 1:
            nft = vft - wft*rft  # this should be SP
            misfits[ib] = np.linalg.norm(nft)**2
        else:
            for it in range(nt):
                nft = vft[it, :] - wft[it, :]*rft
                misfits[ib] = misfits[ib] + np.linalg.norm(nft)**2

        # Compute denominator and GCV function.
        den = (nf*nt - np.real(np.sum(xft)))**2
        gcvf[ib] = misfits[ib]/den

    # Compute best beta.
    ibest = np.argmin(gcvf)

    if iplot:
        import matplotlib.pyplot as mp

        fig, axs = mp.subplots(ncols=2)
        ax = axs[0]
        ax.plot(modnorm, misfits, color='black')
        ax.plot(modnorm, misfits, marker='+', color='red')
        ax.plot(modnorm[ibest], misfits[ibest], marker='o', color='green')
        ax.set_xlabel('Model Norm')
        ax.set_ylabel('Data Misfit')
        ax = axs[1]
        ax.plot(betas, gcvf)
        ax.plot(betas, gcvf, marker='+', color='red')
        ax.plot(betas[ibest], gcvf[ibest], color='green', marker='o')
        ax.set_xscale('log')
        ax.set_xlabel('Regularization Parameter')
        ax.set_ylabel('GCV Function')
        fig.show()
        input('Press key to continue')
        mp.close(fig)

    if ibest == 0 or ibest == len(betas):
        # Minimum not found
        msg = 'No minimum found. Try extending wlevels.'
        raise ValueError(msg)

    wlevel = float(betas[ibest]/np.amax(wwft))
    return wlevel
